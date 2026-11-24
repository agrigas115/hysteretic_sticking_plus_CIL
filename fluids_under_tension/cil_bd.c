/* Run tensioned fluid: hysteretic sticky cells with CIL-C or CIL-P activity */
/* Alex Grigas */
/*
 * Overdamped dynamics at T = 0 with drag 1: dr_i/dt = F_i, forward Euler.
 * Radii are uniform in [0.5, 0.7) and sigma_ij = R_i + R_j.
 *
 * Bonds: an unbonded pair is repulsive while it overlaps and binds on contact. A bond is
 * a double-sided linear spring, F_ij = (r0_ij - r) r_hat_ij, that unbinds at rate koff
 * (tau_b = 1/koff). Contact bonds have r0 = sigma_ij.
 *
 * --model C (CIL-C, crawling): each bonded cell crawls with force f0 directly away from
 *   the sum of the unit vectors to its bonded neighbours.
 * --model P (CIL-P, pulling): no crawling. A cell binds a radical-Voronoi neighbour, at any
 *   distance, at rate kon per pair, if the new bond points into an open angular gap wider
 *   than theta_m on both cells. The new bond starts at r0 = r and ratchets in,
 *   dr0/dt = -kr sigma_ij, until r0 = sigma_ij.
 *
 * Stages: bulk (random start, passive FIRE at phi, t_equil, then t_run of measurement),
 * then optionally a strain step (--shear, --dilation) or a body-force probe (--probe 1,
 * --vF0 transverse or --V0 longitudinal), measured for t_strain. --start runs only the
 * second stage from a saved end state; --strain_only 1 repeats the bulk equilibration
 * and goes straight to it; --extend 1 continues a finished bulk stage to a longer t_run.
 *
 * Outputs per stage: isf.txt, chi4.txt, msd.txt (q = 2 pi sqrt(rho), centre of mass
 * removed), stress_acf.txt, response.txt (pair stress), traj.bin (snapshots with bonds;
 * read with trajectory.py), and vf.txt for a probe. Runs checkpoint and resume exactly.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include "multitau_isf.h"
#include "multitau_scalar.h"

#define pi 3.14159265358979323846
#define MAXBURST 16 // overlapping snapshot bursts
#define MAXDEG 64   // bonds per cell in the angular-gap test
#define MAXV 128    // vertices of a Voronoi cell while clipping

static void fail(const char *msg){
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

static void *xmalloc(size_t n){
	if (n == 0) n = 1;
	void *p = malloc(n);
	if (!p) fail("out of memory");
	return p;
}

static void *xrealloc(void *p, size_t n){
	if (n == 0) n = 1;
	p = realloc(p, n);
	if (!p) fail("out of memory");
	return p;
}

// RNG (splitmix64), uniform in [0, 1)
static double splitmix_unif(uint64_t *state){
	uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return ((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
}

// Main stream: initial positions, bond lifetimes and binding attempts
static uint64_t rng_state;
static double rand_unif(void){
	return splitmix_unif(&rng_state);
}

// Cells, forces and the Verlet list
typedef struct {
	int N;
	double Lx, Ly, inv_Lx, inv_Ly; // periodic box
	double r_l;          // list cutoff / sigma_ij = 1 + skin
	double sig_min, sig_max;
	double *x, *y;       // unwrapped positions
	double *R;           // radii
	double *fx, *fy;     // total force (= velocity, drag 1)
	double *mx, *my;     // CIL-C crawl direction
	int *nb;             // CIL-C bonds per cell
	double *x0, *y0;     // positions at the last list build

	// Verlet list: each pair stores its periodic image (in box lengths), sigma_ij,
	// and its bond (-1 if none). Every bond is in the list.
	int npairs, cap;
	int *pair_i, *pair_j, *pair_bond;
	double *pair_kx, *pair_ky, *pair_sig;

	// Cell list (used when the box holds >= 3 cells per side)
	int use_cells, ncx, ncy;
	double csx, csy;
	int *head, *next;
} sys_t;

// Bonds (i < j), and a min-heap of their expiry times
typedef struct {
	int n, cap;
	int *i, *j, *pair;   // cells and Verlet list entry
	double *r0;          // rest length
	double *deadline;    // expiry time
	int *heap, *hpos, nheap;
	int timed;           // lifetimes are drawn (off during FIRE)
	int ratchet;         // rest lengths can differ from sigma_ij (CIL-P)
	double now, koff;
	long formed, broken; // counts
	long version;        // changes whenever the bonds do
	int *first, *nextb;  // list-build workspace: bonds by their lower cell
	char *found;
} bonds_t;

// CIL-P binding: Voronoi neighbour pairs and the attempt clock
typedef struct {
	double kon, theta_m;
	double clock;        // time (in units of 1/rate) left before the next attempt
	long attempts, binds;
	int n, cap;          // Voronoi pairs (i < j)
	int *i, *j;
	long version;        // bonds version of the partner lists
	int *start, *partner, partner_cap; // bonded partners of each cell
	int nvx, nvy, *vhead, *vnext;      // grid for the Voronoi search
	double *xw, *yw;
} bind_t;

// Pair virial sums, sum over pairs of r_ij,a F_ij,b
typedef struct {
	double Wxx, Wyy, Wxy;
} obs_t;

// Grow the Verlet list arrays
static void grow_pairs(sys_t *s, int n){
	if (n <= s->cap) return;
	while (s->cap < n) s->cap *= 2;
	s->pair_i = xrealloc(s->pair_i, s->cap*sizeof(int));
	s->pair_j = xrealloc(s->pair_j, s->cap*sizeof(int));
	s->pair_bond = xrealloc(s->pair_bond, s->cap*sizeof(int));
	s->pair_kx = xrealloc(s->pair_kx, s->cap*sizeof(double));
	s->pair_ky = xrealloc(s->pair_ky, s->cap*sizeof(double));
	s->pair_sig = xrealloc(s->pair_sig, s->cap*sizeof(double));
}

// Grow the bond arrays
static void grow_bonds(bonds_t *b, int n){
	if (n <= b->cap) return;
	while (b->cap < n) b->cap *= 2;
	b->i = xrealloc(b->i, b->cap*sizeof(int));
	b->j = xrealloc(b->j, b->cap*sizeof(int));
	b->pair = xrealloc(b->pair, b->cap*sizeof(int));
	b->r0 = xrealloc(b->r0, b->cap*sizeof(double));
	b->deadline = xrealloc(b->deadline, b->cap*sizeof(double));
	b->heap = xrealloc(b->heap, b->cap*sizeof(int));
	b->hpos = xrealloc(b->hpos, b->cap*sizeof(int));
	b->nextb = xrealloc(b->nextb, b->cap*sizeof(int));
	b->found = xrealloc(b->found, b->cap);
}

/*--- Bonds ---*/

// Heap order: earlier expiry first, ties broken by pair (i, j)
static int bond_earlier(bonds_t *b, int p, int q){
	if (b->deadline[p] != b->deadline[q]) return b->deadline[p] < b->deadline[q];
	if (b->i[p] != b->i[q]) return b->i[p] < b->i[q];
	return b->j[p] < b->j[q];
}

static void heap_set(bonds_t *b, int at, int k){
	b->heap[at] = k;
	b->hpos[k] = at;
}

static void heap_push(bonds_t *b, int k){
	int at = b->nheap;
	b->nheap += 1;
	while (at > 0 && bond_earlier(b, k, b->heap[(at-1)/2])){
		heap_set(b, at, b->heap[(at-1)/2]);
		at = (at-1)/2;
	}
	heap_set(b, at, k);
}

static int heap_pop(bonds_t *b){
	int top = b->heap[0];
	b->nheap -= 1;
	int last = b->heap[b->nheap];
	int at = 0;
	while (2*at+1 < b->nheap){
		int child = 2*at+1;
		if (child+1 < b->nheap && bond_earlier(b, b->heap[child+1], b->heap[child])) child += 1;
		if (!bond_earlier(b, b->heap[child], last)) break;
		heap_set(b, at, b->heap[child]);
		at = child;
	}
	if (b->nheap > 0) heap_set(b, at, last);
	b->hpos[top] = -1;
	return top;
}

// Exponential lifetime from now, at least one step
static void draw_lifetime(bonds_t *b, int k){
	double life = INFINITY;
	if (b->koff > 0) life = -log1p(-rand_unif())/b->koff;
	if (life == 0) life = nextafter(b->now, INFINITY) - b->now;
	b->deadline[k] = b->now + life;
	if (b->deadline[k] <= b->now) b->deadline[k] = nextafter(b->now, INFINITY);
	heap_push(b, k);
}

// New bond on Verlet list entry `entry` with rest length r0
static void bond_add(sys_t *s, bonds_t *b, int entry, double r0){
	grow_bonds(b, b->n + 1);
	int k = b->n;
	b->n += 1;
	b->i[k] = s->pair_i[entry];
	b->j[k] = s->pair_j[entry];
	b->pair[k] = entry;
	b->r0[k] = r0;
	b->hpos[k] = -1;
	b->deadline[k] = INFINITY;
	b->formed += 1;
	b->version += 1;
	s->pair_bond[entry] = k;
	if (b->timed) draw_lifetime(b, k);
}

// Remove bond k; the last bond moves into its slot
static void bond_remove(sys_t *s, bonds_t *b, int k){
	b->n -= 1;
	int last = b->n;
	b->version += 1;
	s->pair_bond[b->pair[k]] = -1;
	if (k != last){
		b->i[k] = b->i[last];
		b->j[k] = b->j[last];
		b->pair[k] = b->pair[last];
		b->r0[k] = b->r0[last];
		b->deadline[k] = b->deadline[last];
		b->hpos[k] = b->hpos[last];
		if (b->hpos[k] >= 0) b->heap[b->hpos[k]] = k;
		s->pair_bond[b->pair[k]] = k;
	}
}

// Break every bond whose lifetime has run out
static void expire_bonds(sys_t *s, bonds_t *b){
	while (b->nheap > 0 && b->deadline[b->heap[0]] <= b->now){
		int k = heap_pop(b);
		bond_remove(s, b, k);
		b->broken += 1;
	}
}

// Rebuild the heap from every bond, in list order
static void heap_rebuild(bonds_t *b){
	int k;
	b->nheap = 0;
	for (k=0; k<b->n; k++){
		heap_push(b, k);
	}
}

static int cmp_bond_key(const void *p, const void *q){
	const long *u = p, *v = q;
	if (u[0] < v[0]) return -1;
	if (u[0] > v[0]) return 1;
	return 0;
}

// Order bonds by pair (i, j); used once after preparation, before lifetimes are drawn
static void sort_bonds(sys_t *s, bonds_t *b){
	int k, n = b->n;
	long *key = xmalloc(2*(size_t)n*sizeof(long));
	for (k=0; k<n; k++){
		key[2*k] = (long)b->i[k]*s->N + b->j[k];
		key[2*k+1] = b->pair[k];
	}
	qsort(key, n, 2*sizeof(long), cmp_bond_key);
	for (k=0; k<n; k++){
		b->i[k] = (int)(key[2*k] / s->N);
		b->j[k] = (int)(key[2*k] % s->N);
		b->pair[k] = (int)key[2*k+1];
		b->r0[k] = s->pair_sig[b->pair[k]]; // preparation forms contact bonds only
		s->pair_bond[b->pair[k]] = k;
	}
	free(key);
}

/*--- Neighbour list ---*/

// Box size and cell grid
static void set_box(sys_t *s, double Lx, double Ly){
	s->Lx = Lx;
	s->Ly = Ly;
	s->inv_Lx = 1.0/Lx;
	s->inv_Ly = 1.0/Ly;
	double max_cutoff = s->sig_max * s->r_l;
	s->ncx = (int)(Lx / max_cutoff);
	s->ncy = (int)(Ly / max_cutoff);
	s->use_cells = s->ncx >= 3 && s->ncy >= 3;
	if (s->use_cells){
		s->csx = Lx / s->ncx;
		s->csy = Ly / s->ncy;
		s->head = xrealloc(s->head, (size_t)s->ncx*s->ncy*sizeof(int));
	}
}

// Bond between i < j, or -1
static int find_bond(bonds_t *b, int i, int j){
	int k;
	for (k=b->first[i]; k>=0; k=b->nextb[k]){
		if (b->j[k] == j) return k;
	}
	return -1;
}

// Add pair i < j if it is within the list cutoff or bonded. The all-pairs build compares
// distances, the cell build squared distances.
static void add_pair(sys_t *s, bonds_t *b, int i, int j, int by_distance){
	double dx = s->x[i] - s->x[j];
	double dy = s->y[i] - s->y[j];
	double kx = nearbyint(dx * s->inv_Lx);
	double ky = nearbyint(dy * s->inv_Ly);
	dx -= s->Lx * kx;
	dy -= s->Ly * ky;
	double sigma = s->R[i] + s->R[j];
	double cutoff = sigma * s->r_l;
	double d2 = dx*dx + dy*dy;
	int near;
	if (by_distance) near = sqrt(d2) < cutoff;
	else near = d2 < cutoff*cutoff;
	int k = find_bond(b, i, j);
	if (!near && k < 0) return;

	grow_pairs(s, s->npairs + 1);
	int n = s->npairs;
	s->npairs += 1;
	s->pair_i[n] = i;
	s->pair_j[n] = j;
	s->pair_kx[n] = kx;
	s->pair_ky[n] = ky;
	s->pair_sig[n] = sigma;
	s->pair_bond[n] = k;
	if (k >= 0){
		b->pair[k] = n;
		b->found[k] = 1;
	}
}

// Rebuild the Verlet list: all pairs in index order (all_pairs = 1, or a small box), or a
// half-stencil cell list. Bonds stretched beyond the cutoff are added at the end.
static void build_list(sys_t *s, bonds_t *b, int all_pairs){
	int i, j, k, n, N = s->N;
	s->npairs = 0;

	// Bonds by their lower cell, for find_bond
	for (i=0; i<N; i++){
		b->first[i] = -1;
	}
	for (k=b->n-1; k>=0; k--){
		b->nextb[k] = b->first[b->i[k]];
		b->first[b->i[k]] = k;
		b->found[k] = 0;
	}

	if (all_pairs || !s->use_cells){
		for (i=0; i<N; i++){
			for (j=i+1; j<N; j++){
				add_pair(s, b, i, j, 1);
			}
		}
	}
	else{
		int ncx = s->ncx, ncy = s->ncy;
		for (k=0; k<ncx*ncy; k++){
			s->head[k] = -1;
		}
		// Sort cells into the grid
		for (i=0; i<N; i++){
			double xw = s->x[i] - s->Lx * nearbyint(s->x[i] * s->inv_Lx);
			double yw = s->y[i] - s->Ly * nearbyint(s->y[i] * s->inv_Ly);
			int cx = (int)((xw + 0.5*s->Lx) / s->csx);
			int cy = (int)((yw + 0.5*s->Ly) / s->csy);
			if (cx < 0) cx = 0;
			if (cy < 0) cy = 0;
			if (cx >= ncx) cx = ncx-1;
			if (cy >= ncy) cy = ncy-1;
			int c = cy*ncx + cx;
			s->next[i] = s->head[c];
			s->head[c] = i;
		}
		// Half stencil: own cell + 4 neighbours, so each pair of grid cells is visited once
		static const int off[5][2] = {{0,0},{1,0},{-1,1},{0,1},{1,1}};
		int cx, cy, pa, pb;
		for (cy=0; cy<ncy; cy++){
			for (cx=0; cx<ncx; cx++){
				int c = cy*ncx + cx;
				for (n=0; n<5; n++){
					int cx2 = (cx + off[n][0] + ncx) % ncx;
					int cy2 = (cy + off[n][1] + ncy) % ncy;
					int c2 = cy2*ncx + cx2;
					for (pa=s->head[c]; pa>=0; pa=s->next[pa]){
						// In its own grid cell, pa pairs only with the cells after it
						int first = s->head[c2];
						if (n == 0) first = s->next[pa];
						for (pb=first; pb>=0; pb=s->next[pb]){
							if (pa < pb) add_pair(s, b, pa, pb, 0);
							else add_pair(s, b, pb, pa, 0);
						}
					}
				}
			}
		}
		for (k=0; k<b->n; k++){
			if (!b->found[k]) add_pair(s, b, b->i[k], b->j[k], 0);
		}
	}

	memcpy(s->x0, s->x, N*sizeof(double));
	memcpy(s->y0, s->y, N*sizeof(double));
}

/*--- Forces ---*/

// Pair forces and binding on contact. CIL-P: rest lengths ratchet in by kr_dt sigma_ij until
// they reach sigma_ij; CIL-C bonds always have r0 = sigma_ij.
static void pair_forces(sys_t *s, bonds_t *b, double kr_dt, obs_t *o){
	int n, i, j, k, N = s->N;
	memset(s->fx, 0, N*sizeof(double));
	memset(s->fy, 0, N*sizeof(double));
	o->Wxx = 0;
	o->Wyy = 0;
	o->Wxy = 0;

	for (n=0; n<s->npairs; n++){
		i = s->pair_i[n];
		j = s->pair_j[n];
		double dx = s->x[i] - s->x[j];
		double dy = s->y[i] - s->y[j];
		dx -= s->Lx * s->pair_kx[n];
		dy -= s->Ly * s->pair_ky[n];
		double r = sqrt(dx*dx + dy*dy);
		double sigma = s->pair_sig[n];

		// If not bonded, bind on contact
		if (s->pair_bond[n] < 0){
			if (r > sigma) continue;
			bond_add(s, b, n, sigma);
		}
		k = s->pair_bond[n];
		if (!(r > 0)) fail("coincident cells: force direction undefined");

		// Unit vector j -> i
		double inv_r = 1.0 / r;
		double ux = dx * inv_r;
		double uy = dy * inv_r;

		// Linear spring: > 0 repulsive, < 0 attractive
		double r0 = sigma;
		if (b->ratchet) r0 = b->r0[k];
		double mag = r0 - r;
		double Fx = ux * mag;
		double Fy = uy * mag;
		s->fx[i] += Fx;
		s->fy[i] += Fy;
		s->fx[j] += -Fx;
		s->fy[j] += -Fy;

		// Virial
		o->Wxx += Fx*dx;
		o->Wyy += Fy*dy;
		o->Wxy += Fx*dy;

		// Ratchet the rest length
		if (b->ratchet){
			if (b->r0[k] > sigma) b->r0[k] -= kr_dt*sigma;
			else if (b->r0[k] < sigma) b->r0[k] = sigma;
		}
	}
}

// CIL-C: each bonded cell crawls with force f0 away from its bonded neighbours
static void crawl_forces(sys_t *s, bonds_t *b, double f0){
	int i, j, k, N = s->N;

	// Sum the bond unit vectors of each cell (in bond-list order)
	memset(s->mx, 0, N*sizeof(double));
	memset(s->my, 0, N*sizeof(double));
	memset(s->nb, 0, N*sizeof(int));
	for (k=0; k<b->n; k++){
		i = b->i[k];
		j = b->j[k];
		double dx = s->x[i] - s->x[j];
		double dy = s->y[i] - s->y[j];
		dx -= s->Lx * nearbyint(dx * s->inv_Lx);
		dy -= s->Ly * nearbyint(dy * s->inv_Ly);
		double r = sqrt(dx*dx + dy*dy);
		if (!(r > 0)) fail("coincident cells: force direction undefined");
		double inv_r = 1.0 / r;
		double ux = dx * inv_r;
		double uy = dy * inv_r;
		s->mx[i] += ux;
		s->my[i] += uy;
		s->mx[j] += -ux;
		s->my[j] += -uy;
		s->nb[i] += 1;
		s->nb[j] += 1;
	}

	// Normalize and crawl with force f0
	for (i=0; i<N; i++){
		if (s->nb[i] > 0){
			double m = sqrt(s->mx[i]*s->mx[i] + s->my[i]*s->my[i]);
			// Perfectly balanced bonds give no direction
			if (m > 1e-12){
				s->mx[i] /= m;
				s->my[i] /= m;
			}
			else{
				s->mx[i] = 0;
				s->my[i] = 0;
			}
		}
		s->fx[i] += f0*s->mx[i];
		s->fy[i] += f0*s->my[i];
	}
}

// Forward Euler step. Returns 1 if a cell has moved half the smallest list margin since
// the last build, so the list must be rebuilt.
static int euler_step(sys_t *s, double dt){
	int i, rebuild = 0;
	double half_skin = 0.5*s->sig_min*(s->r_l - 1);
	double lim = half_skin*half_skin;
	for (i=0; i<s->N; i++){
		s->x[i] += dt * s->fx[i];
		s->y[i] += dt * s->fy[i];
		if (!isfinite(s->x[i]) || !isfinite(s->y[i])) fail("nonfinite position");
		double dx = s->x[i] - s->x0[i];
		double dy = s->y[i] - s->y0[i];
		if (dx*dx + dy*dy >= lim) rebuild = 1;
	}
	return rebuild;
}

// Passive FIRE at fixed box to relax the random start. Contacts bind as they form
// (no lifetimes yet). Uses the all-pairs list.
static int fire_relax(sys_t *s, bonds_t *b){
	int i, N = s->N;
	obs_t o;
	double *vx = calloc(N, sizeof(double));
	double *vy = calloc(N, sizeof(double));
	if (!vx || !vy) fail("out of memory");
	double alpha0 = 0.25, finc = 1.01, fdec = 0.25, falpha = 0.99;
	double NNEGMAX = 2000, fireitmin = 50;
	double alpha = alpha0, fire_dt = 0.1, dtmax = 0.1, dtmin = 1e-2*fire_dt;
	double Ftol = 1e-14, fcheck = 10*Ftol, min_rc = 0.5;
	int npPos = 0, npNeg = 0, fireit = 0;

	// Zero reference positions force a rebuild on the first iteration
	memset(s->x0, 0, N*sizeof(double));
	memset(s->y0, 0, N*sizeof(double));
	while (fcheck > Ftol){
		if (fireit >= 1000000) fail("FIRE did not converge");
		fireit += 1;
		for (i=0; i<N; i++){
			vx[i] += 0.5*fire_dt*s->fx[i];
			vy[i] += 0.5*fire_dt*s->fy[i];
		}
		double vnorm = 0, fnorm = 0, P = 0;
		for (i=0; i<N; i++){
			vnorm += sqrt(vx[i]*vx[i] + vy[i]*vy[i]);
		}
		for (i=0; i<N; i++){
			fnorm += sqrt(s->fx[i]*s->fx[i] + s->fy[i]*s->fy[i]);
		}
		for (i=0; i<N; i++){
			P += vx[i]*s->fx[i] + vy[i]*s->fy[i];
		}

		if (P > 0){
			npPos += 1;
			npNeg = 0;
			if (npPos > fireitmin && fire_dt*finc < dtmax) fire_dt = fire_dt*finc;
			alpha = alpha*falpha;
		}
		else{
			npPos = 0;
			npNeg += 1;
			if (npNeg > NNEGMAX) fail("FIRE repeatedly failed to descend");
			if (fireit > fireitmin){
				if (fire_dt*fdec > dtmin) fire_dt = fire_dt*fdec;
				alpha = alpha0;
			}
			for (i=0; i<N; i++){
				s->x[i] -= 0.5*fire_dt*vx[i];
				s->y[i] -= 0.5*fire_dt*vy[i];
			}
			memset(vx, 0, N*sizeof(double));
			memset(vy, 0, N*sizeof(double));
		}
		if (fnorm > 0){
			for (i=0; i<N; i++){
				vx[i] = (1 - alpha)*vx[i] + alpha*(s->fx[i]/fnorm)*vnorm;
				vy[i] = (1 - alpha)*vy[i] + alpha*(s->fy[i]/fnorm)*vnorm;
			}
		}
		for (i=0; i<N; i++){
			s->x[i] += 0.5*fire_dt*vx[i];
			s->y[i] += 0.5*fire_dt*vy[i];
		}

		double disp_max = 0;
		for (i=0; i<N; i++){
			double dx = s->x0[i] - s->x[i];
			double dy = s->y0[i] - s->y[i];
			double d = sqrt(dx*dx + dy*dy);
			if (d > disp_max) disp_max = d;
		}
		if (2*disp_max > (min_rc*s->r_l) - min_rc) build_list(s, b, 1);
		pair_forces(s, b, 0, &o);

		for (i=0; i<N; i++){
			vx[i] += 0.5*fire_dt*s->fx[i];
			vy[i] += 0.5*fire_dt*s->fy[i];
		}
		fcheck = 0;
		for (i=0; i<N; i++){
			fcheck += sqrt(s->fx[i]*s->fx[i] + s->fy[i]*s->fy[i]);
		}
		fcheck = fcheck/N;
	}
	free(vx);
	free(vy);
	return fireit;
}

/*--- CIL-P binding ---*/

// Vertex of a Voronoi cell (relative to its cell) and the neighbour that owns the edge to
// the next vertex (-1: the initial square)
typedef struct {
	double x, y;
	int nb;
} vtx_t;

// A neighbour image: displacement from the cell, squared distance, index
typedef struct {
	double dx, dy, d2;
	int id;
} cand_t;

// Keep the part of convex polygon v (n vertices, counter-clockwise) with ax x + ay y <= c;
// the new edge belongs to neighbour id
static int clip_cell(vtx_t *v, int n, double ax, double ay, double c, int id, vtx_t *out){
	int k, m = 0;
	for (k=0; k<n; k++){
		vtx_t *p = &v[k];
		vtx_t *q = &v[(k+1) % n];
		double fp = ax*p->x + ay*p->y - c;
		double fq = ax*q->x + ay*q->y - c;
		if (m + 2 > MAXV) fail("Voronoi cell has too many vertices");
		if (fp <= 0){
			out[m] = *p;
			m += 1;
		}
		if ((fp <= 0) != (fq <= 0)){
			double t = fp / (fp - fq);
			out[m].x = p->x + t*(q->x - p->x);
			out[m].y = p->y + t*(q->y - p->y);
			if (fp <= 0) out[m].nb = id;
			else out[m].nb = p->nb;
			m += 1;
		}
	}
	return m;
}

// Radical (power) Voronoi neighbours in the periodic box (as voro++ container_poly).
// Cell i is cut by each neighbour image j at displacement d by the half-plane
//     2 d.x <= |d|^2 + R_i^2 - R_j^2,
// nearest first. Image j can still cut only if |d| < r + sqrt(r^2 + R_max^2 - R_i^2), with
// r the cell's farthest vertex; the search ring grows until that bound is covered.
static void voronoi_pairs(sys_t *s, bind_t *a){
	int i, k, N = s->N;
	double R_max = 0;
	for (i=0; i<N; i++){
		if (s->R[i] > R_max) R_max = s->R[i];
	}

	// Sort cells into a grid of spacing about 1.2
	a->nvx = (int)(s->Lx / 1.2);
	a->nvy = (int)(s->Ly / 1.2);
	if (a->nvx < 1) a->nvx = 1;
	if (a->nvy < 1) a->nvy = 1;
	int nvx = a->nvx, nvy = a->nvy;
	double csx = s->Lx / nvx, csy = s->Ly / nvy;
	double cmin = csx;
	if (csy < cmin) cmin = csy;
	a->vhead = xrealloc(a->vhead, (size_t)nvx*nvy*sizeof(int));
	for (k=0; k<nvx*nvy; k++){
		a->vhead[k] = -1;
	}
	for (i=0; i<N; i++){
		a->xw[i] = s->x[i] - s->Lx*floor(s->x[i]/s->Lx);
		a->yw[i] = s->y[i] - s->Ly*floor(s->y[i]/s->Ly);
		int cx = (int)(a->xw[i]/csx);
		int cy = (int)(a->yw[i]/csy);
		if (cx >= nvx) cx = nvx-1;
		if (cy >= nvy) cy = nvy-1;
		a->vnext[i] = a->vhead[cy*nvx + cx];
		a->vhead[cy*nvx + cx] = i;
	}

	int ncand = 0, cand_cap = 64;
	cand_t *cand = xmalloc(cand_cap*sizeof(cand_t));
	vtx_t va[MAXV], vb[MAXV];
	a->n = 0;
	for (i=0; i<N; i++){
		int cx = (int)(a->xw[i]/csx);
		int cy = (int)(a->yw[i]/csy);
		if (cx >= nvx) cx = nvx-1;
		if (cy >= nvy) cy = nvy-1;
		double w_i = R_max*R_max - s->R[i]*s->R[i];
		int ring, nv = 0, first_pair = a->n;
		vtx_t *v = va;
		for (ring=2; ; ring++){
			// Every image within cover of cell i
			double cover = ring * cmin;
			ncand = 0;
			int ox, oy;
			for (oy=-ring; oy<=ring; oy++){
				int ky = cy + oy;
				int wy = ((ky % nvy) + nvy) % nvy;
				double sy = s->Ly * floor((double)ky / nvy);
				for (ox=-ring; ox<=ring; ox++){
					int kx = cx + ox;
					int wx = ((kx % nvx) + nvx) % nvx;
					double sx = s->Lx * floor((double)kx / nvx);
					int j;
					for (j=a->vhead[wy*nvx + wx]; j>=0; j=a->vnext[j]){
						double dx = a->xw[j] + sx - a->xw[i];
						double dy = a->yw[j] + sy - a->yw[i];
						if (j == i && dx == 0 && dy == 0) continue;
						if (ncand == cand_cap){
							cand_cap *= 2;
							cand = xrealloc(cand, cand_cap*sizeof(cand_t));
						}
						cand[ncand].dx = dx;
						cand[ncand].dy = dy;
						cand[ncand].d2 = dx*dx + dy*dy;
						cand[ncand].id = j;
						ncand += 1;
					}
				}
			}

			// Nearest first (insertion sort; a few dozen images)
			for (k=1; k<ncand; k++){
				cand_t c = cand[k];
				int m = k - 1;
				while (m >= 0 && cand[m].d2 > c.d2){
					cand[m+1] = cand[m];
					m -= 1;
				}
				cand[m+1] = c;
			}

			// Start from a square of half side cover and clip by each image
			double h = cover;
			va[0].x = -h; va[0].y = -h; va[0].nb = -1;
			va[1].x =  h; va[1].y = -h; va[1].nb = -1;
			va[2].x =  h; va[2].y =  h; va[2].nb = -1;
			va[3].x = -h; va[3].y =  h; va[3].nb = -1;
			v = va;
			nv = 4;
			double rmax2 = 2*h*h;
			for (k=0; k<ncand && nv>0; k++){
				double rmax = sqrt(rmax2);
				double bound = rmax + sqrt(rmax2 + w_i);
				if (cand[k].d2 >= bound*bound) break;
				double c = cand[k].d2 + s->R[i]*s->R[i] - s->R[cand[k].id]*s->R[cand[k].id];
				vtx_t *out = va;
				if (v == va) out = vb;
				nv = clip_cell(v, nv, 2*cand[k].dx, 2*cand[k].dy, c, k, out);
				v = out;
				rmax2 = 0;
				int m;
				for (m=0; m<nv; m++){
					double r2 = v[m].x*v[m].x + v[m].y*v[m].y;
					if (r2 > rmax2) rmax2 = r2;
				}
			}

			// Done if no edge of the square survives and nothing beyond cover can cut
			int open = 0, m;
			for (m=0; m<nv; m++){
				if (v[m].nb < 0) open = 1;
			}
			double rmax = sqrt(rmax2);
			double bound = rmax + sqrt(rmax2 + w_i);
			if (!open && bound <= cover) break;
			if (ring > 64) fail("Voronoi cell did not close");
		}

		// Neighbours j > i of this cell, each once
		for (k=0; k<nv; k++){
			int j = cand[v[k].nb].id, m, seen = 0;
			if (j <= i) continue;
			for (m=first_pair; m<a->n; m++){
				if (a->j[m] == j) seen = 1;
			}
			if (seen) continue;
			if (a->n == a->cap){
				a->cap *= 2;
				a->i = xrealloc(a->i, a->cap*sizeof(int));
				a->j = xrealloc(a->j, a->cap*sizeof(int));
			}
			a->i[a->n] = i;
			a->j[a->n] = j;
			a->n += 1;
		}
	}
	free(cand);
}

// Bonded partners of every cell, rebuilt only when the bonds have changed
static void partner_lists(sys_t *s, bonds_t *b, bind_t *a){
	int i, k, N = s->N;
	if (a->version == b->version) return;
	for (i=0; i<=N; i++){
		a->start[i] = 0;
	}
	for (k=0; k<b->n; k++){
		a->start[b->i[k]+1] += 1;
		a->start[b->j[k]+1] += 1;
	}
	for (i=0; i<N; i++){
		a->start[i+1] += a->start[i];
	}
	if (2*b->n > a->partner_cap){
		a->partner_cap = 2*b->n;
		a->partner = xrealloc(a->partner, a->partner_cap*sizeof(int));
	}
	for (k=0; k<b->n; k++){
		a->partner[a->start[b->i[k]]] = b->j[k];
		a->start[b->i[k]] += 1;
		a->partner[a->start[b->j[k]]] = b->i[k];
		a->start[b->j[k]] += 1;
	}
	for (i=N; i>0; i--){
		a->start[i] = a->start[i-1];
	}
	a->start[0] = 0;
	a->version = b->version;
}

// 1 if direction phi from cell c lies strictly inside an angular gap wider than theta_m
// between the cell's bonds (always with fewer than two bonds)
static int open_gap(sys_t *s, bind_t *a, int c, double phi){
	int k, m, n = a->start[c+1] - a->start[c];
	double th[MAXDEG];
	if (n <= 1) return 1;
	if (n > MAXDEG) fail("too many bonds on one cell for the gap test");

	// Sorted bond angles
	for (k=0; k<n; k++){
		int j = a->partner[a->start[c] + k];
		double dx = s->x[j] - s->x[c];
		double dy = s->y[j] - s->y[c];
		dx -= s->Lx * nearbyint(dx * s->inv_Lx);
		dy -= s->Ly * nearbyint(dy * s->inv_Ly);
		double t = atan2(dy, dx);
		for (m=k; m>0 && th[m-1] > t; m--){
			th[m] = th[m-1];
		}
		th[m] = t;
	}

	// Gaps between neighbouring bonds, then the gap across +-pi
	for (k=0; k<n-1; k++){
		if (th[k+1] - th[k] > a->theta_m && th[k] < phi && phi < th[k+1]) return 1;
	}
	if (2*pi - (th[n-1] - th[0]) > a->theta_m && (th[0] > phi || phi > th[n-1])) return 1;
	return 0;
}

// Verlet list entry of i < j, added if the pair is not in the list
static int list_entry(sys_t *s, int i, int j){
	int n;
	for (n=0; n<s->npairs; n++){
		if (s->pair_i[n] == i && s->pair_j[n] == j) return n;
	}
	grow_pairs(s, s->npairs + 1);
	n = s->npairs;
	s->npairs += 1;
	double dx = s->x[i] - s->x[j];
	double dy = s->y[i] - s->y[j];
	s->pair_i[n] = i;
	s->pair_j[n] = j;
	s->pair_kx[n] = nearbyint(dx * s->inv_Lx);
	s->pair_ky[n] = nearbyint(dy * s->inv_Ly);
	s->pair_sig[n] = s->R[i] + s->R[j];
	s->pair_bond[n] = -1;
	return n;
}

// Try to bind Voronoi pair m: unbonded, with an open gap on both cells. The new bond
// starts unstressed, r0 = r.
static void try_bind(sys_t *s, bonds_t *b, bind_t *a, int m){
	int i = a->i[m], j = a->j[m], k;
	partner_lists(s, b, a);
	for (k=a->start[i]; k<a->start[i+1]; k++){
		if (a->partner[k] == j) return;
	}
	double dx = s->x[j] - s->x[i];
	double dy = s->y[j] - s->y[i];
	dx -= s->Lx * nearbyint(dx * s->inv_Lx);
	dy -= s->Ly * nearbyint(dy * s->inv_Ly);
	if (!open_gap(s, a, i, atan2(dy, dx))) return;
	if (!open_gap(s, a, j, atan2(-dy, -dx))) return;
	bond_add(s, b, list_entry(s, i, j), sqrt(dx*dx + dy*dy));
	a->binds += 1;
}

// Binding attempts during one step: a Poisson process at rate kon x (Voronoi pairs), each
// on a random pair, so every qualifying pair binds at rate kon
static void binding_step(sys_t *s, bonds_t *b, bind_t *a, double dt){
	a->clock -= a->kon * a->n * dt;
	while (a->clock <= 0){
		if (a->n > 0) try_bind(s, b, a, (int)(rand_unif() * a->n));
		a->attempts += 1;
		a->clock += -log1p(-rand_unif());
	}
}

/*--- Body-force probe ---*/

// Each cell feels, in the centre-of-mass frame (Y = mean y),
//     transverse   (--vF0):  F_x,i += F0 sin(q (y_i - Y))
//     longitudinal (--V0):   F_y,i += V0 q sin(q (y_i - Y))
// with q = 2 pi n / Ly. vf.txt columns:
//     t, a = (2/N) sum_i (v_x,i - <v_x>) sin(q (y_i - Y)), U = int a dt,
//     c1 = <cos q(y - Y)>, s1 = <sin q(y - Y)>, Y,
// and with --v_split the sin modes of the applied, pair and crawl forces, x (aF aP aC) and
// y (bF bP bC), the y velocity mode bv, then their time integrals.
typedef struct {
	double F0, V0;       // amplitudes
	int n;               // wavenumber index
	double q, U, Y;
	int on;              // a probe runs
	int stage;           // its stage: 1 (--probe 1) or 0 (during the bulk measurement)
	long every;          // steps between vf.txt rows
	int log;             // log-spaced rows per decade before the first regular row
	unsigned char *logmask; // logmask[m] = 1: write a row at step m < every
	int split;           // also write the force split
	double inst[7];      // aF aP aC bv bF bP bC
	double integ[7];     // their time integrals
	char desc[64];       // directory name after probe_
	FILE *fp;
	char path[1400], ckpath[1400];
} probe_t;

static probe_t pr;

// After the internal forces: centre of mass, the force split, then the applied force
static void probe_apply(sys_t *s, double f0){
	int i, dim, N = s->N;
	double Yc = 0;
	for (i=0; i<N; i++){
		Yc += s->y[i];
	}
	pr.Y = Yc/N;

	// Modes of the applied (F), pair (P) and crawl (C) forces, each about its own mean
	if (pr.split){
		double mF[2] = {0, 0}, mP[2] = {0, 0}, mC[2] = {0, 0};
		for (i=0; i<N; i++){
			double ph = pr.q*(s->y[i] - pr.Y);
			double cxi = f0*s->mx[i];
			double cyi = f0*s->my[i];
			mF[0] += pr.F0*sin(ph);
			mF[1] += pr.V0*pr.q*sin(ph);
			mP[0] += s->fx[i] - cxi;
			mP[1] += s->fy[i] - cyi;
			mC[0] += cxi;
			mC[1] += cyi;
		}
		for (dim=0; dim<2; dim++){
			mF[dim] /= N;
			mP[dim] /= N;
			mC[dim] /= N;
		}
		double sF[2] = {0, 0}, sP[2] = {0, 0}, sC[2] = {0, 0};
		for (i=0; i<N; i++){
			double ph = pr.q*(s->y[i] - pr.Y);
			double sn = sin(ph);
			double cxi = f0*s->mx[i];
			double cyi = f0*s->my[i];
			sF[0] += (pr.F0*sn - mF[0])*sn;
			sF[1] += (pr.V0*pr.q*sn - mF[1])*sn;
			sP[0] += (s->fx[i] - cxi - mP[0])*sn;
			sP[1] += (s->fy[i] - cyi - mP[1])*sn;
			sC[0] += (cxi - mC[0])*sn;
			sC[1] += (cyi - mC[1])*sn;
		}
		pr.inst[0] = 2*sF[0]/N;
		pr.inst[1] = 2*sP[0]/N;
		pr.inst[2] = 2*sC[0]/N;
		pr.inst[4] = 2*sF[1]/N;
		pr.inst[5] = 2*sP[1]/N;
		pr.inst[6] = 2*sC[1]/N;
	}

	// Applied force
	for (i=0; i<N; i++){
		double ph = pr.q*(s->y[i] - pr.Y);
		s->fx[i] += pr.F0*sin(ph);
		s->fy[i] += pr.V0*pr.q*sin(ph);
	}
}

// Measurement step m: mode velocity, density modes, integrals, and a vf.txt row
static void probe_record(sys_t *s, long m, int last, double dt, const char *dir){
	int i, k, N = s->N;
	double mv = 0, av = 0, c1 = 0, s1 = 0;
	for (i=0; i<N; i++){
		mv += s->fx[i];
	}
	mv /= N;
	for (i=0; i<N; i++){
		double ph = pr.q*(s->y[i] - pr.Y);
		av += (s->fx[i] - mv)*sin(ph);
		c1 += cos(ph);
		s1 += sin(ph);
	}
	av *= 2.0/N;
	c1 /= N;
	s1 /= N;
	if (m > 0) pr.U += av*dt;

	// y velocity mode and the integrals of the split
	if (pr.split){
		double my_ = 0, bv = 0;
		for (i=0; i<N; i++){
			my_ += s->fy[i];
		}
		my_ /= N;
		for (i=0; i<N; i++){
			bv += (s->fy[i] - my_)*sin(pr.q*(s->y[i] - pr.Y));
		}
		pr.inst[3] = 2*bv/N;
		if (m > 0){
			for (k=0; k<7; k++){
				pr.integ[k] += pr.inst[k]*dt;
			}
		}
	}

	// Open vf.txt at the first row
	if (!pr.fp){
		snprintf(pr.path, sizeof(pr.path), "%s/vf.txt", dir);
		pr.fp = fopen(pr.path, "w");
		if (!pr.fp) fail("cannot write vf.txt");
		fprintf(pr.fp, "# F0 %.10g V0 %.10g qn %d q %.10g Ly %.10g rho %.10g (centre-of-mass frame)\n",
		        pr.F0, pr.V0, pr.n, pr.q, s->Ly, N/(s->Lx*s->Ly));
		fprintf(pr.fp, "# t a(t) U(t) c1 s1 Y");
		if (pr.split) fprintf(pr.fp, " aF aP aC bv bF bP bC IaF IaP IaC Vy IbF IbP IbC");
		fprintf(pr.fp, "\n");
	}

	// Rows every pr.every steps, at the end, and at the log-spaced early steps
	int write_row = m % pr.every == 0 || last;
	if (pr.logmask && m < pr.every && pr.logmask[m]) write_row = 1;
	if (write_row){
		fprintf(pr.fp, "%.10g %.12g %.12g %.12g %.12g %.12g", m*dt, av, pr.U, c1, s1, pr.Y);
		if (pr.split){
			for (k=0; k<7; k++){
				fprintf(pr.fp, " %.12g", pr.inst[k]);
			}
			for (k=0; k<7; k++){
				fprintf(pr.fp, " %.12g", pr.integ[k]);
			}
		}
		fprintf(pr.fp, "\n");
	}
	if (last){
		fclose(pr.fp);
		pr.fp = NULL;
	}
}

// At every periodic checkpoint: U, the integrals and the vf.txt offset, keyed by the step
// the checkpoint resumes at
static void probe_checkpoint(long next_step){
	int k;
	long offset = -1;
	if (pr.fp){
		fflush(pr.fp);
		offset = ftell(pr.fp);
	}
	FILE *fp = fopen(pr.ckpath, "a");
	if (!fp) fail("cannot write vf_ckpt.txt");
	fprintf(fp, "%ld %.17g %ld", next_step, pr.U, offset);
	if (pr.split){
		for (k=0; k<7; k++){
			fprintf(fp, " %.17g", pr.integ[k]);
		}
	}
	fprintf(fp, "\n");
	if (fclose(fp) != 0) fail("error writing vf_ckpt.txt");
}

/*--- Output ---*/

// Grid of response times and multi-tau lags: 1..31, then 16..31 times powers of 2
static long next_lag(long n){
	if (n == 0) return 1;
	long stride = 1;
	while (n/stride >= 32) stride *= 2;
	return n + stride;
}

// Burst lags: 4 per factor of 2 (1..7, then 8..14 by 2, 16..28 by 4, ...), a subset of
// the response grid
static long next_burst_lag(long n){
	if (n == 0) return 1;
	long stride = 1;
	while (n/stride >= 8) stride *= 2;
	return n + stride;
}

// Snapshot schedule: a structural frame every frame_every steps, and bursts starting every
// burst_every steps with frames at next_burst_lag lags up to burst_span
typedef struct {
	long next_sparse, next_origin;
	long origin[MAXBURST], lag[MAXBURST]; // -1: unused
	long sparse_step, period_step, span_step;
	long next_due, frames;
} snap_t;

static snap_t snap_schedule(long frame_every, long burst_every, long burst_span){
	snap_t sn;
	int k;
	memset(&sn, 0, sizeof(sn));
	sn.sparse_step = frame_every;
	sn.period_step = burst_every;
	sn.span_step = burst_span;
	for (k=0; k<MAXBURST; k++){
		sn.origin[k] = -1;
	}
	return sn;
}

// Bonds of the last written frame (sorted pair indices); n = -1 before the first frame
typedef struct {
	int n, cap;
	uint32_t *pair, *cur, *rem, *add;
} framebonds_t;

static void grow_framebonds(framebonds_t *fb, int n){
	if (n <= fb->cap) return;
	fb->cap = 2*n;
	fb->pair = xrealloc(fb->pair, fb->cap*sizeof(uint32_t));
	fb->cur = xrealloc(fb->cur, fb->cap*sizeof(uint32_t));
	fb->rem = xrealloc(fb->rem, fb->cap*sizeof(uint32_t));
	fb->add = xrealloc(fb->add, fb->cap*sizeof(uint32_t));
}

static void write_bytes(FILE *fp, const void *p, size_t n){
	if (fwrite(p, 1, n, fp) != n) fail("output write failed");
}

// traj.bin (CILTRJ03) header: endian mark, N, model (0 CIL-C, 1 CIL-P), radii
static void write_traj_header(FILE *fp, sys_t *s, int model){
	uint32_t endian = 0x01020304;
	uint32_t n = (uint32_t)s->N;
	uint32_t m = (uint32_t)model;
	if (*(unsigned char*)&endian != 4) fail("traj.bin needs a little-endian host");
	write_bytes(fp, "CILTRJ03", 8);
	write_bytes(fp, &endian, 4);
	write_bytes(fp, &n, 4);
	write_bytes(fp, &m, 4);
	write_bytes(fp, s->R, s->N*sizeof(double));
}

// Index of pair i < j among all N(N-1)/2 pairs
static long pair_index(int i, int j, int N){
	return (long)i*N - (long)i*(i+1)/2 + (j - i - 1);
}

static int cmp_u32(const void *p, const void *q){
	uint32_t u = *(const uint32_t*)p;
	uint32_t v = *(const uint32_t*)q;
	if (u < v) return -1;
	if (u > v) return 1;
	return 0;
}

// One traj.bin frame, if one is due at this step:
//     int64 step; double t, Lx, Ly; uint32 flags (1 structural, 2 burst, 4 keyframe),
//     burst origins, bonds, removed, added, rest lengths; int64 origins[]; double centre of
//     mass (x, y); float (x, y) - centre of mass of each cell; bonds as sorted pair indices
//     (all of them on keyframes, else those removed and added since the last frame); then
//     (pair index, float r0) of each bond not at sigma_ij.
static void write_frame(FILE *fp, snap_t *sn, long step, double dt, sys_t *s, bonds_t *b, framebonds_t *fb, float *xy){
	int k;
	if (step != sn->next_due) return;
	uint32_t flags = 0, nmatch = 0;
	int64_t origins[MAXBURST];

	// Which schedules this frame belongs to
	if (step == sn->next_sparse){
		flags |= 1;
		sn->next_sparse += sn->sparse_step;
	}
	if (step == sn->next_origin){
		int slot = 0;
		while (slot < MAXBURST && sn->origin[slot] >= 0) slot += 1;
		if (slot == MAXBURST) fail("too many overlapping snapshot bursts");
		sn->origin[slot] = step;
		sn->lag[slot] = 0;
		sn->next_origin += sn->period_step;
	}
	for (k=0; k<MAXBURST; k++){
		if (sn->origin[k] >= 0 && sn->origin[k] + sn->lag[k] == step){
			flags |= 2;
			origins[nmatch] = sn->origin[k];
			nmatch += 1;
			sn->lag[k] = next_burst_lag(sn->lag[k]);
			if (sn->lag[k] > sn->span_step) sn->origin[k] = -1;
		}
	}

	// Next frame due
	if (sn->next_sparse < sn->next_origin) sn->next_due = sn->next_sparse;
	else sn->next_due = sn->next_origin;
	for (k=0; k<MAXBURST; k++){
		if (sn->origin[k] >= 0 && sn->origin[k] + sn->lag[k] < sn->next_due) sn->next_due = sn->origin[k] + sn->lag[k];
	}

	// Bonds now, sorted, and what changed since the last frame
	grow_framebonds(fb, b->n + 1);
	for (k=0; k<b->n; k++){
		fb->cur[k] = (uint32_t)pair_index(b->i[k], b->j[k], s->N);
	}
	qsort(fb->cur, b->n, sizeof(uint32_t), cmp_u32);
	int key = (flags & 1) || fb->n < 0;
	uint32_t nrem = 0, nadd = 0;
	if (key) flags |= 4;
	else{
		// Merge the two sorted lists
		int ip = 0, ic = 0;
		while (ip < fb->n || ic < b->n){
			if (ic >= b->n || (ip < fb->n && fb->pair[ip] < fb->cur[ic])){
				fb->rem[nrem] = fb->pair[ip];
				nrem += 1;
				ip += 1;
			}
			else if (ip >= fb->n || fb->cur[ic] < fb->pair[ip]){
				fb->add[nadd] = fb->cur[ic];
				nadd += 1;
				ic += 1;
			}
			else{
				ip += 1;
				ic += 1;
			}
		}
	}
	uint32_t nr0 = 0;
	for (k=0; k<b->n; k++){
		if (b->r0[k] != s->pair_sig[b->pair[k]]) nr0 += 1;
	}

	// Positions relative to the centre of mass
	double cx = 0, cy = 0;
	for (k=0; k<s->N; k++){
		cx += s->x[k];
		cy += s->y[k];
	}
	cx /= s->N;
	cy /= s->N;
	for (k=0; k<s->N; k++){
		xy[2*k] = (float)(s->x[k] - cx);
		xy[2*k+1] = (float)(s->y[k] - cy);
	}

	int64_t index = step;
	uint32_t count = (uint32_t)b->n;
	double t = step*dt;
	write_bytes(fp, &index, 8);
	write_bytes(fp, &t, 8);
	write_bytes(fp, &s->Lx, 8);
	write_bytes(fp, &s->Ly, 8);
	write_bytes(fp, &flags, 4);
	write_bytes(fp, &nmatch, 4);
	write_bytes(fp, &count, 4);
	write_bytes(fp, &nrem, 4);
	write_bytes(fp, &nadd, 4);
	write_bytes(fp, &nr0, 4);
	write_bytes(fp, origins, nmatch*8);
	write_bytes(fp, &cx, 8);
	write_bytes(fp, &cy, 8);
	write_bytes(fp, xy, 2*(size_t)s->N*sizeof(float));
	if (key){
		write_bytes(fp, fb->cur, b->n*sizeof(uint32_t));
	}
	else{
		write_bytes(fp, fb->rem, nrem*sizeof(uint32_t));
		write_bytes(fp, fb->add, nadd*sizeof(uint32_t));
	}
	for (k=0; k<b->n; k++){
		if (b->r0[k] != s->pair_sig[b->pair[k]]){
			uint32_t pk = (uint32_t)pair_index(b->i[k], b->j[k], s->N);
			float r0 = (float)b->r0[k];
			write_bytes(fp, &pk, 4);
			write_bytes(fp, &r0, 4);
		}
	}
	memcpy(fb->pair, fb->cur, b->n*sizeof(uint32_t));
	fb->n = b->n;
	sn->frames += 1;
}

/*--- End states (equilibrated.txt, endpoint.txt), read back by --start ---*/

// Text formats.
//   CIL-C, CILC_BD_STATE_2: header; RNG, bond count; radius and position of each cell;
//          bonds as pair indices; EXPIRY and each bond's remaining lifetime.
//   CIL-P, CILP_BD_STATE_1: header; RNG, bond count, binding clock; radius and position
//          of each cell; each bond as pair index, rest length and remaining lifetime.
static const char *endpoint_magic(int model){
	if (model == 0) return "CILC_BD_STATE_2";
	return "CILP_BD_STATE_1";
}

static void write_endpoint(const char *path, int model, sys_t *s, bonds_t *b, bind_t *a, double f0, double kr, double dt, double age){
	int k;
	FILE *fp = fopen(path, "w");
	if (!fp) fail("could not write endpoint");
	if (model == 0){
		fprintf(fp, "CILC_BD_STATE_2\n%d %.17g %.17g %.17g %.17g %.17g %.17g\n", s->N, s->Lx, s->Ly, f0, b->koff, dt, age);
		fprintf(fp, "%" PRIu64 " %d\n", rng_state, b->n);
	}
	else{
		fprintf(fp, "CILP_BD_STATE_1\n%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
		        s->N, s->Lx, s->Ly, a->kon, kr, a->theta_m, b->koff, dt, age);
		fprintf(fp, "%" PRIu64 " %d %.17g\n", rng_state, b->n, a->clock);
	}
	for (k=0; k<s->N; k++){
		fprintf(fp, "%.17g %.17g %.17g\n", s->R[k], s->x[k], s->y[k]);
	}
	if (model == 0){
		for (k=0; k<b->n; k++){
			fprintf(fp, "%ld\n", pair_index(b->i[k], b->j[k], s->N));
		}
		fprintf(fp, "EXPIRY %d\n", b->n);
		for (k=0; k<b->n; k++){
			fprintf(fp, "%ld %.17g\n", pair_index(b->i[k], b->j[k], s->N), b->deadline[k] - b->now);
		}
	}
	else{
		for (k=0; k<b->n; k++){
			fprintf(fp, "%ld %.17g %.17g\n", pair_index(b->i[k], b->j[k], s->N), b->r0[k], b->deadline[k] - b->now);
		}
	}
	if (ferror(fp) || fclose(fp) != 0) fail("error writing endpoint");
}

// Number of cells in an end state of this model
static int endpoint_N(const char *path, int model){
	char magic[64];
	int N;
	FILE *fp = fopen(path, "r");
	if (!fp) fail("could not open --start file");
	if (fscanf(fp, "%63s %d", magic, &N) != 2) fail("--start file is not an end state");
	if (strcmp(magic, endpoint_magic(model)) != 0 || N < 1) fail("--start file is not an end state of this model");
	fclose(fp);
	return N;
}

// Pair index -> (i, j)
static void pair_cells(long p, int N, int *i, int *j){
	int lo = 0, hi = N-2;
	while (lo < hi){
		int mid = (lo + hi + 1)/2;
		if (pair_index(mid, mid+1, N) <= p) lo = mid;
		else hi = mid - 1;
	}
	*i = lo;
	*j = (int)(p - pair_index(lo, lo+1, N)) + lo + 1;
}

// Read an end state: box, radii, positions, bonds (lifetimes from now = 0), RNG and
// (CIL-P) binding clock. CIL-C bonds get r0 = sigma_ij.
static void read_endpoint(const char *path, int model, sys_t *s, bonds_t *b, bind_t *a, double f0, double kr, double *Lx, double *Ly){
	char magic[64], method[16];
	int N, k, nb, count, ok;
	double f0_in, kon_in, kr_in, theta_in, koff_in, dt_in, age;
	FILE *fp = fopen(path, "r");
	if (!fp) fail("could not open --start file");

	// Header: must match this run's parameters
	if (model == 0){
		ok = fscanf(fp, "%63s %d %lf %lf %lf %lf %lf %lf", magic, &N, Lx, Ly, &f0_in, &koff_in, &dt_in, &age) == 8;
		ok = ok && N == s->N && f0_in == f0 && koff_in == b->koff;
		ok = ok && fscanf(fp, "%" SCNu64 " %d", &rng_state, &nb) == 2 && nb >= 0;
		if (!ok) fail("--start file does not match f0 and koff");
	}
	else{
		ok = fscanf(fp, "%63s %d %lf %lf %lf %lf %lf %lf %lf %lf", magic, &N, Lx, Ly, &kon_in, &kr_in, &theta_in, &koff_in, &dt_in, &age) == 10;
		ok = ok && N == s->N && kon_in == a->kon && kr_in == kr && theta_in == a->theta_m && koff_in == b->koff;
		ok = ok && fscanf(fp, "%" SCNu64 " %d %lf", &rng_state, &nb, &a->clock) == 3 && nb >= 0;
		if (!ok) fail("--start file does not match kon, kr, theta_m and koff");
	}

	// Cells
	for (k=0; k<N; k++){
		if (fscanf(fp, "%lf %lf %lf", &s->R[k], &s->x[k], &s->y[k]) != 3) fail("--start file truncated");
		if (!(s->R[k] >= 0.5)) fail("radii must be at least 0.5");
	}

	// Bonds
	b->n = 0;
	for (k=0; k<nb; k++){
		long p;
		grow_bonds(b, b->n + 1);
		if (model == 0){
			if (fscanf(fp, "%ld", &p) != 1) fail("--start file has a bad bond");
		}
		else{
			if (fscanf(fp, "%ld %lf %lf", &p, &b->r0[k], &b->deadline[k]) != 3) fail("--start file has a bad bond");
			if (!(b->deadline[k] > 0)) fail("--start file has a bad bond");
		}
		if (p < 0 || p >= (long)N*(N-1)/2) fail("--start file has a bad bond");
		pair_cells(p, N, &b->i[b->n], &b->j[b->n]);
		if (model == 0) b->r0[k] = s->R[b->i[k]] + s->R[b->j[k]];
		b->hpos[b->n] = -1;
		b->n += 1;
	}

	// CIL-C lifetimes
	if (model == 0){
		if (fscanf(fp, "%15s %d", method, &count) != 2 || strcmp(method, "EXPIRY") != 0 || count != nb) fail("--start file has no bond lifetimes");
		for (k=0; k<nb; k++){
			long p;
			if (fscanf(fp, "%ld %lf", &p, &b->deadline[k]) != 2) fail("--start file has a bad bond lifetime");
			if (p != pair_index(b->i[k], b->j[k], N) || !(b->deadline[k] > 0)) fail("--start file has a bad bond lifetime");
		}
	}
	fclose(fp);
}

/*--- Checkpoint ---*/

// SIGTERM (e.g. HTCondor eviction) asks for a checkpoint at the end of the current step
static volatile sig_atomic_t got_sigterm = 0;
static void on_sigterm(int sig){
	(void)sig;
	got_sigterm = 1;
}

// Run parameters; a checkpoint must match them
typedef struct {
	int model, N, init_fire, strain_type, traj, from_start, strain_only;
	unsigned long seed;
	double phi, f0, kon, kr, theta_m, koff, dt, skin, t_equil, t_run, strain, t_strain;
	long log_every, origin_every, frame_every, burst_every, burst_span, voronoi_every;
} run_t;

// Everything else needed to continue exactly; outputs are cut back to the offsets
typedef struct {
	int stage, done, rebuild, npairs, nbonds, nheap, frame_bonds, nvor;
	long next_step, next_response;
	long off_log, off_resp, off_traj;
	long formed, broken, attempts, binds;
	double Lx, Ly, clock;
	uint64_t rng;
	snap_t snap;
} ckpt_t;

#define CKPT_MAGIC "CILBDCK1"

static void sync_file(FILE *fp){
	if (!fp) return;
	if (fflush(fp) != 0 || fsync(fileno(fp)) != 0) fail("output flush failed");
}

// Read n items; 1 if all were read
static int read_items(FILE *fp, void *ptr, size_t size, size_t n){
	return fread(ptr, size, n, fp) == n;
}

static void write_checkpoint(const char *dir, const char *name, run_t *p, ckpt_t *c, sys_t *s, bonds_t *b, bind_t *a, framebonds_t *fb, MTCorrelator *mt, MTScalarCorrelator *sc){
	char path[1200], tmp[1300];
	int N = s->N;
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	snprintf(tmp, sizeof(tmp), "%s.tmp", path);
	c->npairs = s->npairs;
	c->nbonds = b->n;
	c->nheap = b->nheap;
	c->frame_bonds = fb->n;
	c->nvor = a->n;
	c->formed = b->formed;
	c->broken = b->broken;
	c->attempts = a->attempts;
	c->binds = a->binds;
	c->clock = a->clock;
	c->Lx = s->Lx;
	c->Ly = s->Ly;
	c->rng = rng_state;

	// Write to a temporary file, then rename, so a crash never leaves half a checkpoint
	FILE *fp = fopen(tmp, "wb");
	if (!fp) fail("could not write checkpoint");
	fwrite(CKPT_MAGIC, 1, 8, fp);
	fwrite(p, sizeof(run_t), 1, fp);
	fwrite(c, sizeof(ckpt_t), 1, fp);
	fwrite(s->x, sizeof(double), N, fp);
	fwrite(s->y, sizeof(double), N, fp);
	fwrite(s->x0, sizeof(double), N, fp);
	fwrite(s->y0, sizeof(double), N, fp);
	fwrite(s->pair_i, sizeof(int), s->npairs, fp);
	fwrite(s->pair_j, sizeof(int), s->npairs, fp);
	fwrite(s->pair_bond, sizeof(int), s->npairs, fp);
	fwrite(s->pair_sig, sizeof(double), s->npairs, fp);
	fwrite(s->pair_kx, sizeof(double), s->npairs, fp);
	fwrite(s->pair_ky, sizeof(double), s->npairs, fp);
	fwrite(b->i, sizeof(int), b->n, fp);
	fwrite(b->j, sizeof(int), b->n, fp);
	fwrite(b->pair, sizeof(int), b->n, fp);
	fwrite(b->deadline, sizeof(double), b->n, fp);
	fwrite(b->hpos, sizeof(int), b->n, fp);
	fwrite(b->heap, sizeof(int), b->nheap, fp);
	fwrite(b->r0, sizeof(double), b->n, fp);
	fwrite(a->i, sizeof(int), a->n, fp);
	fwrite(a->j, sizeof(int), a->n, fp);
	MTCorrelator copy = *mt;
	copy.buf = NULL;
	fwrite(&copy, sizeof(copy), 1, fp);
	fwrite(mt->buf, sizeof(double), (size_t)MT_NLEV*MT_MBUF*N*2, fp);
	fwrite(sc, sizeof(*sc), 1, fp);
	if (fb->n > 0) fwrite(fb->pair, sizeof(uint32_t), fb->n, fp);
	if (ferror(fp) || fflush(fp) != 0 || fsync(fileno(fp)) != 0 || fclose(fp) != 0) fail("error writing checkpoint");
	if (rename(tmp, path) != 0) fail("could not rename checkpoint");
}

// Run length stored in a checkpoint (-1 if there is none)
static double checkpoint_t_run(const char *dir, const char *name){
	char path[1200], magic[8];
	run_t p_ck;
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	FILE *fp = fopen(path, "rb");
	if (!fp) return -1;
	int ok = read_items(fp, magic, 1, 8) && memcmp(magic, CKPT_MAGIC, 8) == 0;
	ok = ok && read_items(fp, &p_ck, sizeof(run_t), 1);
	fclose(fp);
	if (!ok) return -1;
	return p_ck.t_run;
}

// Returns 1 and fills the state if a checkpoint exists; exits if it belongs to another run.
// extend: the bulk-end state of a shorter run, so t_run (and the log interval and strain
// length that follow from it) may differ.
static int read_checkpoint(const char *dir, const char *name, int extend, run_t *p, ckpt_t *c, sys_t *s, bonds_t *b, bind_t *a, framebonds_t *fb, MTCorrelator *mt, MTScalarCorrelator *sc){
	char path[1200], magic[8];
	int N = s->N;
	run_t p_ck;
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	FILE *fp = fopen(path, "rb");
	if (!fp) return 0;

	// Parameters
	int ok = read_items(fp, magic, 1, 8) && memcmp(magic, CKPT_MAGIC, 8) == 0;
	ok = ok && read_items(fp, &p_ck, sizeof(run_t), 1);
	if (ok && extend){
		run_t q = *p;
		q.t_run = p_ck.t_run;
		q.t_strain = p_ck.t_strain;
		q.log_every = p_ck.log_every;
		ok = memcmp(&p_ck, &q, sizeof(run_t)) == 0;
	}
	else{
		ok = ok && memcmp(&p_ck, p, sizeof(run_t)) == 0;
	}
	if (!ok){
		fprintf(stderr, "%s does not match this run's parameters; delete it to start over\n", path);
		exit(1);
	}

	// Sizes
	ok = read_items(fp, c, sizeof(ckpt_t), 1);
	ok = ok && c->npairs >= 0 && c->nbonds >= 0 && c->nheap >= 0 && c->nheap <= c->nbonds && c->nvor >= 0;
	if (ok){
		grow_pairs(s, c->npairs);
		grow_bonds(b, c->nbonds);
		s->npairs = c->npairs;
		b->n = c->nbonds;
		b->nheap = c->nheap;
		if (c->nvor > a->cap){
			a->cap = c->nvor;
			a->i = xrealloc(a->i, a->cap*sizeof(int));
			a->j = xrealloc(a->j, a->cap*sizeof(int));
		}
		a->n = c->nvor;
	}

	// Arrays, in the order written
	ok = ok && read_items(fp, s->x, sizeof(double), N);
	ok = ok && read_items(fp, s->y, sizeof(double), N);
	ok = ok && read_items(fp, s->x0, sizeof(double), N);
	ok = ok && read_items(fp, s->y0, sizeof(double), N);
	ok = ok && read_items(fp, s->pair_i, sizeof(int), s->npairs);
	ok = ok && read_items(fp, s->pair_j, sizeof(int), s->npairs);
	ok = ok && read_items(fp, s->pair_bond, sizeof(int), s->npairs);
	ok = ok && read_items(fp, s->pair_sig, sizeof(double), s->npairs);
	ok = ok && read_items(fp, s->pair_kx, sizeof(double), s->npairs);
	ok = ok && read_items(fp, s->pair_ky, sizeof(double), s->npairs);
	ok = ok && read_items(fp, b->i, sizeof(int), b->n);
	ok = ok && read_items(fp, b->j, sizeof(int), b->n);
	ok = ok && read_items(fp, b->pair, sizeof(int), b->n);
	ok = ok && read_items(fp, b->deadline, sizeof(double), b->n);
	ok = ok && read_items(fp, b->hpos, sizeof(int), b->n);
	ok = ok && read_items(fp, b->heap, sizeof(int), b->nheap);
	ok = ok && read_items(fp, b->r0, sizeof(double), b->n);
	ok = ok && read_items(fp, a->i, sizeof(int), a->n);
	ok = ok && read_items(fp, a->j, sizeof(int), a->n);
	ok = ok && read_items(fp, mt, sizeof(*mt), 1);
	if (ok){
		mt->buf = xmalloc((size_t)MT_NLEV*MT_MBUF*N*2*sizeof(double));
		ok = read_items(fp, mt->buf, sizeof(double), (size_t)MT_NLEV*MT_MBUF*N*2);
	}
	ok = ok && read_items(fp, sc, sizeof(*sc), 1);
	if (ok){
		fb->n = c->frame_bonds;
		if (fb->n > 0){
			grow_framebonds(fb, fb->n + 1);
			ok = read_items(fp, fb->pair, sizeof(uint32_t), fb->n);
		}
	}
	fclose(fp);
	if (!ok){
		fprintf(stderr, "%s is truncated or corrupt\n", path);
		exit(1);
	}

	rng_state = c->rng;
	b->formed = c->formed;
	b->broken = c->broken;
	a->attempts = c->attempts;
	a->binds = c->binds;
	a->clock = c->clock;
	a->version = -1; // partner lists are rebuilt from the bonds
	set_box(s, c->Lx, c->Ly);
	return 1;
}

// Wall-clock seconds (a node made slow by a busy host still reports normal CPU time)
static double wall_seconds(void){
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1e-9*ts.tv_nsec;
}

// Number of times this job has left a slow node, kept beside the checkpoint
static int read_slow_moves(const char *dir){
	char path[1200];
	int n = 0;
	snprintf(path, sizeof(path), "%s/slow_moves.txt", dir);
	FILE *fp = fopen(path, "r");
	if (fp){
		if (fscanf(fp, "%d", &n) != 1) n = 0;
		fclose(fp);
	}
	return n;
}

static void write_slow_moves(const char *dir, int n){
	char path[1200];
	snprintf(path, sizeof(path), "%s/slow_moves.txt", dir);
	FILE *fp = fopen(path, "w");
	if (fp){
		fprintf(fp, "%d\n", n);
		fclose(fp);
	}
}

// Open an output file for appending, cut back to the checkpointed offset
static FILE *reopen_at(const char *path, long off){
	struct stat st;
	FILE *fp = fopen(path, "r+");
	int ok = fp != NULL;
	ok = ok && fstat(fileno(fp), &st) == 0 && st.st_size >= off;
	ok = ok && ftruncate(fileno(fp), off) == 0 && fseek(fp, 0, SEEK_END) == 0;
	if (!ok){
		fprintf(stderr, "could not resume %s\n", path);
		exit(1);
	}
	return fp;
}

static FILE *open_write(const char *dir, const char *name){
	char path[1300];
	snprintf(path, sizeof(path), "%s/%s", dir, name);
	FILE *fp = fopen(path, "w");
	if (!fp){
		fprintf(stderr, "could not write %s\n", path);
		exit(1);
	}
	return fp;
}

// mkdir -p
static void make_dir(const char *path){
	char tmp[1200];
	char *c;
	snprintf(tmp, sizeof(tmp), "%s", path);
	for (c=tmp+1; *c; c++){
		if (*c == '/'){
			*c = 0;
			mkdir(tmp, 0755);
			*c = '/';
		}
	}
	if (mkdir(tmp, 0755) != 0 && errno != EEXIST){
		fprintf(stderr, "could not create %s\n", path);
		exit(1);
	}
}

static void usage(void){
	fprintf(stderr,
		"usage: cil_bd --model C|P [options]\n"
		"  --model C|P     CIL-C (crawling) or CIL-P (pulling) (C)\n"
		"  --seed S        RNG seed for radii, positions, lifetimes and binding (1)\n"
		"  --N N           number of cells (256)\n"
		"  --phi X         packing fraction (0.4)\n"
		"  --koff X        unbinding rate, tau_b = 1/koff (0.001)\n"
		"  --f0 X          CIL-C: crawling force (1)\n"
		"  --kon X         CIL-P: binding rate per open-gap Voronoi pair (0.021)\n"
		"  --kr X          CIL-P: ratchet rate, dr0/dt = -kr sigma_ij (0.02)\n"
		"  --theta_m X     CIL-P: smallest open angular gap for binding, radians (2)\n"
		"  --voronoi_every n  CIL-P: steps between Voronoi neighbour updates (500)\n"
		"  --dt X          time step (0.02 CIL-C, 0.01 CIL-P)\n"
		"  --init M        fire (passive FIRE at fixed phi) or random (fire)\n"
		"  --skin X        Verlet list margin as a fraction of sigma_ij (0.4)\n"
		"  --t_equil X     dynamics before the bulk measurement (0)\n"
		"  --t_run X       measurement time of the bulk stage (100)\n"
		"  --t_strain X    measurement time of the strain or probe stage (t_run)\n"
		"  --shear G       after the bulk stage, step pure shear G and measure again\n"
		"  --dilation G    after the bulk stage, step area log-strain G and measure again\n"
		"  --strain_only 1 with --shear/--dilation (or --probe): equilibrate as the bulk job does,\n"
		"                  then go straight to the strain stage, in <seed>/shear_G or dilation_G\n"
		"  --start FILE    skip the bulk stage: strain (or probe) this saved bulk end state\n"
		"                  (writes to <dir of FILE>/shear_G, dilation_G or probe_...)\n"
		"  --extend 1      continue a finished job's bulk stage to the (longer) --t_run, from\n"
		"                  bulk/extend.bin; other arguments as the job's. Writes the bulk\n"
		"                  outputs as one run of that length would, and runs no strain stage\n"
		"  --probe 1       body-force probe in place of the strain step, with --vF0 or --V0 and\n"
		"                  --start FILE or --strain_only 1: rebuild the lists as a zero-amplitude\n"
		"                  step does, then run for t_strain in <seed>/probe_<T|L>_n<n>_<F|V><amp>\n"
		"  --vF0 X         transverse probe amplitude F0\n"
		"  --V0 X          longitudinal probe amplitude V0 (peak force V0 q)\n"
		"  --vqn n         probe wavevector q = 2 pi n / Ly (1)\n"
		"  --v_every n     steps between vf.txt rows (50)\n"
		"  --v_log n       also n log-spaced vf.txt rows per decade before the first (0)\n"
		"  --v_split 0|1   also write the applied, pair and crawl force modes (0)\n"
		"  --log_every n   steps between log lines (t_run/1e4 or t_equil/5e3, at least 1/dt)\n"
		"  --origin_every n  steps between short-lag correlator origins (16)\n"
		"  --frame_every n   steps between structural frames (0.1 tau_b)\n"
		"  --burst_every n   steps between snapshot bursts (0.5 tau_b)\n"
		"  --burst_span n    longest burst lag in steps (3 tau_b)\n"
		"  --traj 0|1      write traj.bin (1)\n"
		"  --out DIR       output root directory (data)\n"
		"  --ckpt_minutes X   wall-clock minutes between checkpoints (30; 0 = only at the end)\n"
		"  --ckpt_exit_code C exit with code C after each checkpoint (0 = keep running)\n"
		"  --slow_ns X        leave a slow node: after --slow_check_minutes of wall time in a phase,\n"
		"                     if a step takes more than X ns per cell, checkpoint and exit with\n"
		"                     --slow_exit_code so the scheduler can move the job (0 = off)\n"
		"  --slow_check_minutes M   when to measure (5)\n"
		"  --slow_max_moves K       stop moving after K slow exits and finish where it is (3)\n"
		"  --slow_exit_code C       exit code for a slow-node exit (99)\n"
		"A run resumes from <dir>/checkpoint.bin if present, and checkpoints on SIGTERM.\n");
	exit(1);
}

/*--- Stages ---*/

// Output files and correlators of the current stage
typedef struct {
	char dir[1200];
	FILE *resp, *traj;
	MTCorrelator mt;
	MTScalarCorrelator sc;
} stage_t;

// Directory name of a step: shear_<G>, dilation_<G>, or probe_<T|L>_n<n>_<F|V><amp>
static void kind_name(char *out, size_t n, int strain_type, double strain){
	if (strain_type == 3) snprintf(out, n, "probe_%s", pr.desc);
	else if (strain_type == 1) snprintf(out, n, "shear_%g", strain);
	else snprintf(out, n, "dilation_%g", strain);
}

// Directory of a stage
static void stage_dir(char *out, size_t n, const char *job_dir, int stage, int strain_type, double strain, int from_start){
	char kn[96];
	if (stage == 0){
		snprintf(out, n, "%s/bulk", job_dir);
	}
	else if (from_start){
		snprintf(out, n, "%s", job_dir);
	}
	else{
		kind_name(kn, sizeof(kn), strain_type, strain);
		snprintf(out, n, "%s/%s", job_dir, kn);
	}
}

// Fresh correlators and output files for a stage
static void open_stage(stage_t *st, sys_t *s, run_t *p){
	make_dir(st->dir);
	mt_init(&st->mt, s->N, 2*pi*sqrt(s->N/(s->Lx*s->Ly)), 1, p->dt);
	st->mt.origin_stride = p->origin_every;
	mts_init(&st->sc, 1, p->dt);
	st->resp = open_write(st->dir, "response.txt");
	fprintf(st->resp, "# t Pxx Pyy Pxy z\n");
	st->traj = NULL;
	if (p->traj){
		st->traj = open_write(st->dir, "traj.bin");
		write_traj_header(st->traj, s, p->model);
	}
}

// Flush the log and stage files and note their sizes in the checkpoint state
static void mark_outputs(ckpt_t *ck, FILE *fp_log, stage_t *st){
	sync_file(fp_log);
	sync_file(st->resp);
	sync_file(st->traj);
	ck->off_log = ftell(fp_log);
	ck->off_resp = 0;
	ck->off_traj = 0;
	if (st->resp) ck->off_resp = ftell(st->resp);
	if (st->traj) ck->off_traj = ftell(st->traj);
}

// Affine step on the bulk end state; bond lifetimes are re-based to now = 0
static void apply_strain(sys_t *s, bonds_t *b, int strain_type, double strain){
	int i, k;
	double sx = exp(strain/2);
	double sy;
	if (strain_type == 2) sy = exp(strain/2);  // dilation
	else sy = exp(-strain/2);                  // pure shear (or a zero-amplitude probe step)
	for (k=0; k<b->n; k++){
		b->deadline[k] -= b->now;
	}
	b->now = 0;
	heap_rebuild(b);
	double Lx = s->Lx*sx, Ly = s->Ly*sy;
	for (i=0; i<s->N; i++){
		s->x[i] *= sx;
		s->y[i] *= sy;
	}
	set_box(s, Lx, Ly);
	build_list(s, b, 0);
}

// One log.txt line
static void write_log_line(FILE *fp, run_t *p, sys_t *s, bonds_t *b, bind_t *a, obs_t *o, long step, long m, int phase){
	int i, N = s->N;
	double vol_box = s->Lx*s->Ly;
	double mean_v = 0, max_v = 0;
	for (i=0; i<N; i++){
		double v = hypot(s->fx[i], s->fy[i]);
		if (!isfinite(v)) fail("nonfinite force");
		mean_v += v;
		if (v > max_v) max_v = v;
	}
	fprintf(fp, "%ld\t%.17g\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g", step, m*p->dt, phase,
	        o->Wxx/vol_box, o->Wyy/vol_box, o->Wxy/vol_box, 2.0*b->n/N, mean_v/N, max_v);
	if (p->model == 1) fprintf(fp, "\t%ld\t%ld\t%ld\t%ld\t%d", b->formed, b->broken, a->binds, a->attempts, a->n);
	fprintf(fp, "\n");
}

int main(int argc, char *argv[]){
	int i, k;

	// Default parameters //
	run_t p;
	memset(&p, 0, sizeof(p));
	p.model = 0;
	p.seed = 1;
	p.N = 256;
	p.phi = 0.4;
	p.f0 = 1;
	p.kon = 0.021;
	p.kr = 0.02;
	p.theta_m = 2.0;
	p.koff = 0.001;
	p.dt = -1; // model default below
	p.init_fire = 1;
	p.skin = 0.4;
	p.t_equil = 0;
	p.t_run = 100;
	p.t_strain = -1;
	p.origin_every = 16;
	p.voronoi_every = 500;
	p.traj = 1;
	pr.n = 1;
	pr.every = 50;
	const char *out_root = "data";
	const char *start = NULL;
	double ckpt_minutes = 30;
	int ckpt_exit_code = 0;
	int extend = 0, probe_mode = 0;
	int c_opt = 0, p_opt = 0; // options of one model were given
	double slow_ns = 0, slow_check_minutes = 5;
	int slow_max_moves = 3, slow_exit_code = 99;

	// Read arguments //
	for (k=1; k<argc; k++){
		if (k+1 >= argc) usage();
		char *opt = argv[k];
		char *val = argv[k+1];
		k += 1;
		if (!strcmp(opt, "--model")){
			if (!strcmp(val, "C") || !strcmp(val, "c")) p.model = 0;
			else if (!strcmp(val, "P") || !strcmp(val, "p")) p.model = 1;
			else usage();
		}
		else if (!strcmp(opt, "--seed"))          p.seed = strtoul(val, NULL, 10);
		else if (!strcmp(opt, "--N"))             p.N = atoi(val);
		else if (!strcmp(opt, "--phi"))           p.phi = atof(val);
		else if (!strcmp(opt, "--koff"))          p.koff = atof(val);
		else if (!strcmp(opt, "--f0")){           p.f0 = atof(val); c_opt = 1; }
		else if (!strcmp(opt, "--kon")){          p.kon = atof(val); p_opt = 1; }
		else if (!strcmp(opt, "--kr")){           p.kr = atof(val); p_opt = 1; }
		else if (!strcmp(opt, "--theta_m")){      p.theta_m = atof(val); p_opt = 1; }
		else if (!strcmp(opt, "--voronoi_every")){ p.voronoi_every = atol(val); p_opt = 1; }
		else if (!strcmp(opt, "--dt"))            p.dt = atof(val);
		else if (!strcmp(opt, "--init"))          p.init_fire = strcmp(val, "random") != 0;
		else if (!strcmp(opt, "--skin"))          p.skin = atof(val);
		else if (!strcmp(opt, "--t_equil"))       p.t_equil = atof(val);
		else if (!strcmp(opt, "--t_run"))         p.t_run = atof(val);
		else if (!strcmp(opt, "--t_strain"))      p.t_strain = atof(val);
		else if (!strcmp(opt, "--shear")){        p.strain_type = 1; p.strain = atof(val); }
		else if (!strcmp(opt, "--dilation")){     p.strain_type = 2; p.strain = atof(val); }
		else if (!strcmp(opt, "--strain_only"))   p.strain_only = atoi(val);
		else if (!strcmp(opt, "--start"))         start = val;
		else if (!strcmp(opt, "--extend"))        extend = atoi(val);
		else if (!strcmp(opt, "--probe"))         probe_mode = atoi(val);
		else if (!strcmp(opt, "--vF0"))           pr.F0 = atof(val);
		else if (!strcmp(opt, "--V0"))            pr.V0 = atof(val);
		else if (!strcmp(opt, "--vqn"))           pr.n = atoi(val);
		else if (!strcmp(opt, "--v_every"))       pr.every = atol(val);
		else if (!strcmp(opt, "--v_log"))         pr.log = atoi(val);
		else if (!strcmp(opt, "--v_split"))       pr.split = atoi(val);
		else if (!strcmp(opt, "--log_every"))     p.log_every = atol(val);
		else if (!strcmp(opt, "--origin_every"))  p.origin_every = atol(val);
		else if (!strcmp(opt, "--frame_every"))   p.frame_every = atol(val);
		else if (!strcmp(opt, "--burst_every"))   p.burst_every = atol(val);
		else if (!strcmp(opt, "--burst_span"))    p.burst_span = atol(val);
		else if (!strcmp(opt, "--traj"))          p.traj = atoi(val);
		else if (!strcmp(opt, "--out"))           out_root = val;
		else if (!strcmp(opt, "--ckpt_minutes"))  ckpt_minutes = atof(val);
		else if (!strcmp(opt, "--ckpt_exit_code")) ckpt_exit_code = atoi(val);
		else if (!strcmp(opt, "--slow_ns"))       slow_ns = atof(val);
		else if (!strcmp(opt, "--slow_check_minutes")) slow_check_minutes = atof(val);
		else if (!strcmp(opt, "--slow_max_moves")) slow_max_moves = atoi(val);
		else if (!strcmp(opt, "--slow_exit_code")) slow_exit_code = atoi(val);
		else usage();
	}

	// The other model's parameters are refused, and zeroed so they enter neither the
	// dynamics nor the checkpoint
	if (p.model == 0){
		if (p_opt) fail("--kon, --kr, --theta_m and --voronoi_every are CIL-P options (--model P)");
		p.kon = 0;
		p.kr = 0;
		p.theta_m = 0;
		p.voronoi_every = 0;
	}
	else{
		if (c_opt) fail("--f0 is a CIL-C option (--model C)");
		p.f0 = 0;
	}
	if (p.dt < 0){
		if (p.model == 0) p.dt = 0.02;
		else p.dt = 0.01;
	}
	rng_state = p.seed; // --start replaces it with the saved stream

	// Probe: its own stage, entered like a strain step of zero amplitude. Without
	// --probe 1, a probe runs during the bulk measurement.
	pr.on = pr.F0 != 0 || pr.V0 != 0;
	if (probe_mode){
		if (!pr.on) fail("--probe needs --vF0 or --V0");
		if (p.strain_type) fail("--probe replaces --shear/--dilation");
		if (!start && !p.strain_only) fail("--probe needs --start FILE or --strain_only 1");
		p.strain_type = 3;
		p.strain = 0;
		pr.stage = 1;
		if (pr.F0 != 0 && pr.V0 != 0) snprintf(pr.desc, sizeof(pr.desc), "TL_n%d_F%g_V%g", pr.n, pr.F0, pr.V0);
		else if (pr.F0 != 0) snprintf(pr.desc, sizeof(pr.desc), "T_n%d_F%g", pr.n, pr.F0);
		else snprintf(pr.desc, sizeof(pr.desc), "L_n%d_V%g", pr.n, pr.V0);
	}
	char kn[96];
	kind_name(kn, sizeof(kn), p.strain_type, p.strain);
	if (start){
		p.from_start = 1;
		p.N = endpoint_N(start, p.model);
		if (!p.strain_type) fail("--start needs --shear, --dilation or --probe");
	}

	// Check parameters //
	int N = p.N;
	if (N < 1 || p.dt <= 0 || p.koff < 0 || p.phi <= 0 || p.t_equil < 0 || p.t_run < 0) fail("bad parameters");
	if (p.f0 < 0 || p.kon < 0 || p.kr < 0 || p.theta_m < 0) fail("bad parameters");
	if (!(p.skin > 0 && p.skin <= 1) || p.origin_every < 1 || fabs(p.strain) > 2) fail("bad parameters");
	if (p.model == 1 && p.voronoi_every < 1) fail("bad parameters");
	if (p.t_strain < 0) p.t_strain = p.t_run;
	if (extend && start) fail("--extend continues a bulk stage; it cannot be combined with --start");
	if (p.strain_only && (!p.strain_type || start || extend)) fail("--strain_only needs --shear, --dilation or --probe, and no --start or --extend");
	long n_equil = lround(p.t_equil/p.dt);
	long n_run = lround(p.t_run/p.dt);
	long n_strain = lround(p.t_strain/p.dt);
	if (fabs(n_equil*p.dt - p.t_equil) > 1e-9*fmax(1, p.t_equil)) fail("t_equil must be a multiple of dt");
	if (fabs(n_run*p.dt - p.t_run) > 1e-9*fmax(1, p.t_run)) fail("t_run must be a multiple of dt");
	if (fabs(n_strain*p.dt - p.t_strain) > 1e-9*fmax(1, p.t_strain)) fail("t_strain must be a multiple of dt");
	if (fmax(p.t_run, p.t_strain) > p.dt*(MT_MBUF-1)*(double)(1L << (MT_NLEV-1))) fail("t_run exceeds the multi-tau lag range; increase MT_NLEV");

	// Output schedules: frames every 0.1 tau_b, bursts every 0.5 tau_b spanning 3 tau_b
	double tau_ref = fmax(p.t_run, 1);
	if (p.koff > 0) tau_ref = 1/p.koff;
	long unit = lround(0.1*tau_ref/p.dt);
	if (unit < 1) unit = 1;
	if (p.frame_every == 0) p.frame_every = unit;
	if (p.burst_every == 0) p.burst_every = 5*unit;
	if (p.burst_span == 0) p.burst_span = 30*unit;
	if (p.log_every == 0) p.log_every = (long)ceil(fmax(1., fmax(p.t_run/10000, p.t_equil/5000))/p.dt);
	if (p.frame_every < 1 || p.burst_every < 1 || p.burst_span < 1 || p.log_every < 1) fail("bad output schedule");
	if ((double)p.burst_span/p.burst_every >= MAXBURST - 1) fail("burst_span / burst_every must be below 15");

	// Allocate //
	sys_t s;
	bonds_t b;
	bind_t a;
	memset(&s, 0, sizeof(s));
	memset(&b, 0, sizeof(b));
	memset(&a, 0, sizeof(a));
	s.N = N;
	s.r_l = 1.0 + p.skin;
	s.x = xmalloc(N*sizeof(double));
	s.y = xmalloc(N*sizeof(double));
	s.R = xmalloc(N*sizeof(double));
	s.fx = xmalloc(N*sizeof(double));
	s.fy = xmalloc(N*sizeof(double));
	s.mx = calloc(N, sizeof(double));
	s.my = calloc(N, sizeof(double));
	if (!s.mx || !s.my) fail("out of memory");
	s.x0 = xmalloc(N*sizeof(double));
	s.y0 = xmalloc(N*sizeof(double));
	s.nb = xmalloc(N*sizeof(int));
	s.next = xmalloc(N*sizeof(int));
	s.cap = 16*N;
	s.pair_i = xmalloc(s.cap*sizeof(int));
	s.pair_j = xmalloc(s.cap*sizeof(int));
	s.pair_bond = xmalloc(s.cap*sizeof(int));
	s.pair_kx = xmalloc(s.cap*sizeof(double));
	s.pair_ky = xmalloc(s.cap*sizeof(double));
	s.pair_sig = xmalloc(s.cap*sizeof(double));
	b.koff = p.koff;
	b.ratchet = p.model == 1;
	b.cap = 4*N;
	b.i = xmalloc(b.cap*sizeof(int));
	b.j = xmalloc(b.cap*sizeof(int));
	b.pair = xmalloc(b.cap*sizeof(int));
	b.r0 = xmalloc(b.cap*sizeof(double));
	b.deadline = xmalloc(b.cap*sizeof(double));
	b.heap = xmalloc(b.cap*sizeof(int));
	b.hpos = xmalloc(b.cap*sizeof(int));
	b.nextb = xmalloc(b.cap*sizeof(int));
	b.found = xmalloc(b.cap);
	b.first = xmalloc(N*sizeof(int));
	a.kon = p.kon;
	a.theta_m = p.theta_m;
	a.version = -1;
	a.cap = 4*N;
	a.i = xmalloc(a.cap*sizeof(int));
	a.j = xmalloc(a.cap*sizeof(int));
	a.start = xmalloc((N + 1)*sizeof(int));
	a.partner_cap = 8*N;
	a.partner = xmalloc(a.partner_cap*sizeof(int));
	a.vnext = xmalloc(N*sizeof(int));
	a.xw = xmalloc(N*sizeof(double));
	a.yw = xmalloc(N*sizeof(double));
	double kr_dt = p.kr*p.dt;
	float *xy = xmalloc(2*(size_t)N*sizeof(float));
	framebonds_t fb;
	memset(&fb, 0, sizeof(fb));
	fb.n = -1;
	double (*sample)[2] = xmalloc(2*(size_t)N*sizeof(double));

	// Radii and job directory //
	double Lx = 0, Ly = 0;
	char job_dir[1100];
	if (start){
		// From a saved end state, written beside it
		read_endpoint(start, p.model, &s, &b, &a, p.f0, p.kr, &Lx, &Ly);
		char base[1024];
		snprintf(base, sizeof(base), "%s", start);
		char *slash = strrchr(base, '/');
		if (slash) *slash = 0;
		else snprintf(base, sizeof(base), ".");
		snprintf(job_dir, sizeof(job_dir), "%s/%s", base, kn);
	}
	else{
		// Radii uniform in [0.5, 0.7), from their own stream so positions do not depend on them
		uint64_t radius_rng = p.seed*0x2545F4914F6CDD1DULL + 12345;
		for (i=0; i<N; i++){
			s.R[i] = 0.5*(1.0 + 0.4*splitmix_unif(&radius_rng));
		}
		if (p.model == 0) snprintf(job_dir, sizeof(job_dir), "%s/phi_%g_f0_%g_koff_%g_N_%d_dt_%g/seed_%lu", out_root, p.phi, p.f0, p.koff, N, p.dt, p.seed);
		else snprintf(job_dir, sizeof(job_dir), "%s/phi_%g_kon_%g_kr_%g_koff_%g_N_%d_dt_%g/seed_%lu", out_root, p.phi, p.kon, p.kr, p.koff, N, p.dt, p.seed);
		// A strain-only job gets its own directory, beside the bulk job's stages
		if (p.strain_only){
			size_t l = strlen(job_dir);
			snprintf(job_dir + l, sizeof(job_dir) - l, "/%s", kn);
		}
	}

	// Box size, and the smallest and largest sigma_ij //
	double vol = 0, R_lo = INFINITY, R_lo2 = INFINITY, R_hi = 0, R_hi2 = 0;
	for (i=0; i<N; i++){
		vol += pi*s.R[i]*s.R[i];
		if (s.R[i] < R_lo){
			R_lo2 = R_lo;
			R_lo = s.R[i];
		}
		else if (s.R[i] < R_lo2){
			R_lo2 = s.R[i];
		}
		if (s.R[i] > R_hi){
			R_hi2 = R_hi;
			R_hi = s.R[i];
		}
		else if (s.R[i] > R_hi2){
			R_hi2 = s.R[i];
		}
	}
	if (N > 1){
		s.sig_min = R_lo + R_lo2;
		s.sig_max = R_hi + R_hi2;
	}
	else{
		s.sig_min = 2*s.R[0];
		s.sig_max = 2*s.R[0];
	}
	if (!start){
		Lx = sqrt(vol / p.phi);
		Ly = Lx;
	}
	set_box(&s, Lx, Ly);

	// Probe wavevector and early rows (m = round(10^(j/log)) below the first regular row)
	pr.q = 2*pi*pr.n/s.Ly;
	if (pr.every < 1 || pr.log < 0) fail("bad --v_every or --v_log");
	if (pr.log > 0 && pr.every > 1){
		pr.logmask = calloc((size_t)pr.every, 1);
		if (!pr.logmask) fail("out of memory");
		int jl;
		for (jl=0; ; jl++){
			long ml = lround(pow(10.0, (double)jl/pr.log));
			if (ml >= pr.every) break;
			pr.logmask[ml] = 1;
		}
	}
	if (pr.on && !probe_mode && (p.strain_type || start || extend)) fail("body-force probe: bulk runs or --probe 1 only");
	if (probe_mode && extend) fail("--probe cannot be combined with --extend");
	make_dir(job_dir);

	// Checkpoint directory: the job directory, or for an extension the bulk stage
	// directory, so the finished job's own files stay as they are
	char ckdir[1200];
	const char *ckname = "checkpoint.bin";
	if (extend){
		ckname = "checkpoint_extend.bin";
		snprintf(ckdir, sizeof(ckdir), "%s/bulk", job_dir);
	}
	else{
		snprintf(ckdir, sizeof(ckdir), "%s", job_dir);
	}

	// Parameter and radii files //
	const char *strain_name = "none";
	if (p.strain_type == 1) strain_name = "shear";
	if (p.strain_type == 2) strain_name = "dilation";
	if (p.strain_type == 3) strain_name = "probe";
	FILE *fp;
	if (extend) fp = open_write(ckdir, "params_extend.txt");
	else fp = open_write(ckdir, "params.txt");
	fprintf(fp, "model %s\nseed %lu\nN %d\nphi %.17g\n", p.model == 0 ? "C" : "P", p.seed, N, p.phi);
	if (p.model == 0) fprintf(fp, "f0 %.17g\n", p.f0);
	else fprintf(fp, "kon %.17g\nkr %.17g\ntheta_m %.17g\n", p.kon, p.kr, p.theta_m);
	fprintf(fp, "koff %.17g\ndt %.17g\ninit %s\nskin %.17g\n", p.koff, p.dt, p.init_fire ? "fire" : "random", p.skin);
	fprintf(fp, "t_equil %.17g\nt_run %.17g\nt_strain %.17g\n", p.t_equil, p.t_run, p.t_strain);
	fprintf(fp, "strain_type %s\nstrain %.17g\nstrain_only %d\nstart %s\n", strain_name, p.strain, p.strain_only, start ? start : "none");
	fprintf(fp, "log_every %ld\norigin_every %ld\nframe_every %ld\nburst_every %ld\nburst_span %ld\n",
	        p.log_every, p.origin_every, p.frame_every, p.burst_every, p.burst_span);
	if (p.model == 1) fprintf(fp, "voronoi_every %ld\n", p.voronoi_every);
	fprintf(fp, "traj %d\n", p.traj);
	fclose(fp);
	fp = open_write(job_dir, "radii.txt");
	for (i=0; i<N; i++){
		fprintf(fp, "%.17g\n", s.R[i]);
	}
	fclose(fp);
	if (p.model == 0){
		printf("CIL-C  N = %d  Lx = %.4f  Ly = %.4f  koff = %g  tau_b = %g  f0 = %g\n", N, s.Lx, s.Ly, p.koff, 1/p.koff, p.f0);
	}
	else{
		printf("CIL-P  N = %d  Lx = %.4f  Ly = %.4f  koff = %g  tau_b = %g  kon = %g  kr = %g  theta_m = %g\n",
		       N, s.Lx, s.Ly, p.koff, 1/p.koff, p.kon, p.kr, p.theta_m);
	}
	printf("output: %s\n", job_dir);

	// Resume from a checkpoint, if there is one //
	ckpt_t ck;
	stage_t st;
	memset(&ck, 0, sizeof(ck));
	memset(&st, 0, sizeof(st));
	int resumed = 0, extend_fresh = 0;
	if (extend){
		// An extension in progress (or done) at this length resumes from its own checkpoint;
		// otherwise it starts from the bulk-end state of the last completed bulk stage
		if (checkpoint_t_run(ckdir, ckname) == p.t_run){
			resumed = read_checkpoint(ckdir, ckname, 0, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
		}
		else{
			double t_done = checkpoint_t_run(ckdir, "extend.bin");
			if (t_done < 0) fail("--extend needs bulk/extend.bin from a finished bulk stage");
			if (t_done >= p.t_run){
				printf("bulk stage already measured for %g >= t_run (%s/extend.bin)\n", t_done, ckdir);
				return 0;
			}
			resumed = read_checkpoint(ckdir, "extend.bin", 1, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
			extend_fresh = resumed;
			if (ck.stage != 0) fail("bulk/extend.bin is not a bulk state");
			ck.done = 0;
			printf("extending the bulk stage from t = %g to %g\n", t_done, p.t_run);
		}
	}
	else{
		resumed = read_checkpoint(ckdir, ckname, 0, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
	}
	if (resumed && ck.done){
		printf("run already complete (%s/%s)\n", ckdir, ckname);
		return 0;
	}
	signal(SIGTERM, on_sigterm);
	signal(SIGINT, on_sigterm);
	b.timed = 1;
	int rebuild = 0;
	obs_t o;
	FILE *fp_log;
	int from_dir = p.from_start || p.strain_only; // second stage written to job_dir itself

	if (resumed){
		// Reopen the outputs at the checkpointed sizes
		rebuild = ck.rebuild;
		stage_dir(st.dir, sizeof(st.dir), job_dir, ck.stage, p.strain_type, p.strain, from_dir);
		char path[1300];
		if (extend_fresh){
			fp_log = open_write(ckdir, "log_extend.txt");
			fprintf(fp_log, "# columns as in log.txt, which holds the run up to the bulk-end state this extension starts from\n");
		}
		else{
			if (extend) snprintf(path, sizeof(path), "%s/log_extend.txt", ckdir);
			else snprintf(path, sizeof(path), "%s/log.txt", ckdir);
			fp_log = reopen_at(path, ck.off_log);
		}
		st.resp = NULL;
		st.traj = NULL;
		// A strain-only job's equilibration has no outputs
		if (!(p.strain_only && ck.stage == 0)){
			snprintf(path, sizeof(path), "%s/response.txt", st.dir);
			st.resp = reopen_at(path, ck.off_resp);
			if (p.traj){
				snprintf(path, sizeof(path), "%s/traj.bin", st.dir);
				st.traj = reopen_at(path, ck.off_traj);
			}
		}
		printf("resuming from checkpoint: stage %d, step %ld\n", ck.stage, ck.next_step);
	}
	else{
		fp_log = open_write(job_dir, "log.txt");
		fprintf(fp_log, "# step\tt\tphase\tPxx\tPyy\tPxy\tz\tmean_v\tmax_v");
		if (p.model == 1) fprintf(fp_log, "\tformed\tbroken\tbinds\tattempts\tvoronoi_pairs");
		fprintf(fp_log, "   (phase 0 equil, 1 bulk, 2 strain; t from the start of the measurement");
		if (p.model == 1) fprintf(fp_log, "; counts since the start of the job");
		fprintf(fp_log, ")\n");

		if (start){
			// The bulk stage was run by another job; lifetimes are already from now = 0
			b.now = 0;
			apply_strain(&s, &b, p.strain_type, p.strain);
			ck.stage = 1;
		}
		else{
			// Random (overlapping) positions at phi, relaxed by passive FIRE at the same box
			b.timed = 0;
			for (i=0; i<N; i++){
				s.x[i] = Lx*rand_unif() - Lx/2.;
				s.y[i] = Lx*rand_unif() - Lx/2.;
			}
			build_list(&s, &b, 1);
			pair_forces(&s, &b, 0, &o);
			int fireit = 0;
			if (p.init_fire) fireit = fire_relax(&s, &b);
			if (!p.init_fire && b.n == 0 && p.f0 > 0) fail("random start has no contacts; choose another seed or a higher phi");

			// Then bond lifetimes and the first binding attempt
			sort_bonds(&s, &b);
			build_list(&s, &b, 0);
			b.timed = 1;
			b.now = 0;
			for (k=0; k<b.n; k++){
				draw_lifetime(&b, k);
			}
			if (a.kon > 0) a.clock = -log1p(-rand_unif());
			printf("prepared: %s, %d FIRE iterations, %d bonds\n", p.init_fire ? "fire" : "random", fireit, b.n);
			ck.stage = 0;
		}
		stage_dir(st.dir, sizeof(st.dir), job_dir, ck.stage, p.strain_type, p.strain, from_dir);
		if (p.strain_only && ck.stage == 0){
			// Equilibration only: correlators for the checkpoint, no files
			mt_init(&st.mt, s.N, 2*pi*sqrt(s.N/(s.Lx*s.Ly)), 1, p.dt);
			mts_init(&st.sc, 1, p.dt);
			st.resp = NULL;
			st.traj = NULL;
		}
		else{
			open_stage(&st, &s, &p);
		}
		ck.snap = snap_schedule(p.frame_every, p.burst_every, p.burst_span);
		ck.next_step = 0;
		ck.next_response = 0;
	}
	fflush(stdout);

	// A resumed probe restores U, the integrals and the vf.txt offset of the matching step
	snprintf(pr.ckpath, sizeof(pr.ckpath), "%s/vf_ckpt.txt", ckdir);
	if (pr.on && resumed){
		FILE *cf = fopen(pr.ckpath, "r");
		long use_off = -2;
		char line[1024];
		while (cf && fgets(line, sizeof(line), cf)){
			long stp, off;
			double U, I7[7];
			int nr = sscanf(line, "%ld %lf %ld %lf %lf %lf %lf %lf %lf %lf", &stp, &U, &off, &I7[0], &I7[1], &I7[2], &I7[3], &I7[4], &I7[5], &I7[6]);
			if (nr >= 3 && stp == ck.next_step){
				pr.U = U;
				use_off = off;
				if (pr.split){
					if (nr != 10) fail("body-force probe: vf_ckpt.txt has no --v_split integrals");
					for (k=0; k<7; k++){
						pr.integ[k] = I7[k];
					}
				}
			}
		}
		if (cf) fclose(cf);
		long probe_start = 0;
		if (pr.stage == 0) probe_start = n_equil;
		if (use_off == -2 && ck.stage == pr.stage && ck.next_step > probe_start) fail("body-force probe: no vf_ckpt.txt entry for this checkpoint");
		if (use_off >= 0){
			snprintf(pr.path, sizeof(pr.path), "%s/vf.txt", st.dir);
			pr.fp = reopen_at(pr.path, use_off);
		}
	}
	double t_last_ckpt = wall_seconds();
	int slow_moves = read_slow_moves(ckdir);
	int slow_check = slow_ns > 0 && slow_moves < slow_max_moves;
	double t_wall0 = wall_seconds();
	long step0 = ck.next_step;
	clock_t t_start = clock();

	// Run the stages //
	for (;;){
		// Steps of this stage: equilibration (bulk only), then measurement
		long n_eq = 0;
		long n_end;
		if (ck.stage == 0){
			n_eq = n_equil;
			n_end = n_equil;
			if (!p.strain_only) n_end += n_run;
		}
		else{
			n_end = n_strain;
		}

		long step;
		for (step=ck.next_step; step<=n_end; step++){
			b.now = step*p.dt;

			// Unbinding, then CIL-P binding (Voronoi pairs every voronoi_every steps)
			if (step > 0) expire_bonds(&s, &b);
			if (a.kon > 0){
				if (step % p.voronoi_every == 0) voronoi_pairs(&s, &a);
				if (step > 0) binding_step(&s, &b, &a, p.dt);
			}

			// Forces
			if (rebuild) build_list(&s, &b, 0);
			pair_forces(&s, &b, kr_dt, &o);
			if (p.f0 > 0) crawl_forces(&s, &b, p.f0);
			int probing = pr.on && ck.stage == pr.stage && step >= n_eq;
			if (probing) probe_apply(&s, p.f0);

			// The equilibrated state, the common start of the bulk and strain jobs
			if (ck.stage == 0 && step == n_equil && !p.from_start && !extend){
				char path_eq[1300];
				snprintf(path_eq, sizeof(path_eq), "%s/equilibrated.txt", job_dir);
				write_endpoint(path_eq, p.model, &s, &b, &a, p.f0, p.kr, p.dt, p.t_equil);
			}

			// Measurements
			long m = step - n_eq; // measurement step
			if (probing && st.resp) probe_record(&s, m, step == n_end, p.dt, st.dir);
			if (step % p.log_every == 0 || step == n_eq || step == n_end){
				int phase = 0;
				if (m >= 0) phase = 1 + ck.stage;
				write_log_line(fp_log, &p, &s, &b, &a, &o, step, m, phase);
			}
			if (m >= 0 && st.resp){
				double vol_box = s.Lx*s.Ly;
				// Pair stress on the multi-tau grid
				if (m == ck.next_response || step == n_end){
					fprintf(st.resp, "%.17g %.17g %.17g %.17g %.17g\n", m*p.dt, o.Wxx/vol_box, o.Wyy/vol_box, o.Wxy/vol_box, 2.0*b.n/N);
					ck.next_response = next_lag(m);
				}
				if (st.traj) write_frame(st.traj, &ck.snap, m, p.dt, &s, &b, &fb, xy);

				// Correlators, centre of mass removed
				double cx = 0, cy = 0;
				for (i=0; i<N; i++){
					cx += s.x[i];
					cy += s.y[i];
				}
				cx /= N;
				cy /= N;
				for (i=0; i<N; i++){
					sample[i][0] = s.x[i] - cx;
					sample[i][1] = s.y[i] - cy;
				}
				mt_push(&st.mt, sample);
				mts_push(&st.sc, o.Wxy/vol_box);
			}
			if (step == n_end) break;

			// Move
			rebuild = euler_step(&s, p.dt);

			// Slow node? Measured once per phase, slow_check_minutes into it. The window
			// restarts with the measurement, which costs more per step than equilibration.
			if (step + 1 == n_eq){
				t_wall0 = wall_seconds();
				step0 = step + 1;
				slow_check = slow_ns > 0 && slow_moves < slow_max_moves;
			}
			int slow_exit = 0;
			if (slow_check && step % 1000 == 0){
				double el = wall_seconds() - t_wall0;
				if (el >= 60*slow_check_minutes){
					double ns = 1e9*el / ((double)(step + 1 - step0) * N);
					printf("speed check: %.0f ns per cell-step (threshold %g), %s\n", ns, slow_ns, m < 0 ? "equilibration" : "measurement");
					if (ns > slow_ns){
						slow_exit = 1;
						slow_moves += 1;
						write_slow_moves(ckdir, slow_moves);
						printf("slow node: moving (move %d of %d)\n", slow_moves, slow_max_moves);
					}
					fflush(stdout);
					slow_check = 0;
				}
			}

			// Bulk-end state for a later --extend: everything at the start of the last bulk step
			if (ck.stage == 0 && step + 1 == n_end && st.resp){
				mark_outputs(&ck, fp_log, &st);
				ck.next_step = step + 1;
				ck.done = 0;
				ck.rebuild = rebuild;
				write_checkpoint(st.dir, "extend.bin", &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
			}

			// Checkpoint on a wall-clock schedule, when asked to stop (SIGTERM), or to leave
			// a slow node
			int due = got_sigterm || slow_exit;
			if (ckpt_minutes > 0 && wall_seconds() - t_last_ckpt >= 60*ckpt_minutes) due = 1;
			if (due){
				mark_outputs(&ck, fp_log, &st);
				ck.next_step = step + 1;
				ck.done = 0;
				ck.rebuild = rebuild;
				write_checkpoint(ckdir, ckname, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
				if (pr.on && ck.stage == pr.stage) probe_checkpoint(ck.next_step);
				t_last_ckpt = wall_seconds();
				if (got_sigterm || slow_exit || ckpt_exit_code > 0){
					printf("checkpointed at stage %d, step %ld, exiting\n", ck.stage, step + 1);
					fflush(stdout);
					if (got_sigterm) return 143;
					if (slow_exit) return slow_exit_code;
					return ckpt_exit_code;
				}
			}
		}

		// End of stage: correlators //
		char path[1300], path2[1300], path3[1300];
		snprintf(path, sizeof(path), "%s/isf.txt", st.dir);
		snprintf(path2, sizeof(path2), "%s/chi4.txt", st.dir);
		snprintf(path3, sizeof(path3), "%s/msd.txt", st.dir);
		if (st.resp){
			mt_dump(&st.mt, path, path2, path3);
			snprintf(path, sizeof(path), "%s/stress_acf.txt", st.dir);
			mts_dump(&st.sc, path);
			if (fclose(st.resp) != 0) fail("error closing response.txt");
			if (st.traj && fclose(st.traj) != 0) fail("error closing traj.bin");
		}

		// Bulk end state, for branching other amplitudes with --start
		if (ck.stage == 0 && !extend && !p.strain_only){
			snprintf(path, sizeof(path), "%s/endpoint.txt", job_dir);
			write_endpoint(path, p.model, &s, &b, &a, p.f0, p.kr, p.dt, p.t_equil + p.t_run);
		}
		if (ck.stage == 1 || !p.strain_type || extend) break;
		mt_free(&st.mt);

		// Strain stage, from the same state //
		apply_strain(&s, &b, p.strain_type, p.strain);
		rebuild = 0;
		ck.stage = 1;
		stage_dir(st.dir, sizeof(st.dir), job_dir, 1, p.strain_type, p.strain, from_dir);
		open_stage(&st, &s, &p);
		ck.snap = snap_schedule(p.frame_every, p.burst_every, p.burst_span);
		fb.n = -1;
		ck.next_step = 0;
		ck.next_response = 0;
		t_wall0 = wall_seconds();
		step0 = 0;
		slow_check = slow_ns > 0 && slow_moves < slow_max_moves;
		mark_outputs(&ck, fp_log, &st);
		ck.rebuild = 0;
		write_checkpoint(ckdir, ckname, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
	}

	// Final checkpoint marks the run complete, so a requeued job exits at once
	sync_file(fp_log);
	ck.done = 1;
	ck.off_log = ftell(fp_log);
	write_checkpoint(ckdir, ckname, &p, &ck, &s, &b, &a, &fb, &st.mt, &st.sc);
	fclose(fp_log);
	mt_free(&st.mt);
	printf("cpu time %.1f s\n", (double)(clock()-t_start)/CLOCKS_PER_SEC);
	return 0;
}

/* Run attractive particles with CIL-C activity */
/* Alex Grigas */
/*
 * Overdamped dynamics at T = 0 with drag 1: dr_i/dt = F_i, forward Euler.
 * Each cell has a repulsive core of radius R_i (diameters uniform in [1, 1.4]) and an
 * attractive shell of radius lambda R_i. For a pair at distance r, with sigma_ij = R_i + R_j,
 *     F_ij = (sigma_ij - r) r_hat_ij   if r < lambda sigma_ij,   0 otherwise,
 * so the pair is repulsive while the cores overlap and attractive within the shells. Such a
 * pair is a contact. CIL-C: each cell crawls with force v0 directly away from the sum of
 * its contact directions.
 *
 * The shell packing fraction phi_a sets lambda = sqrt(phi_a/phi). The largest attractive
 * pair force is (lambda - 1) sigma_ij.
 *
 * A run starts from random positions, relaxes the cores without attraction or crawling for
 * t_relax, then runs for t_run. Outputs: log.txt (contacts, speeds, bond turnover, pair
 * stress), window.txt (fraction of cells moving relative to their contacts, per window tau),
 * traj.txt (frames of positions and crawling directions). Runs checkpoint and resume exactly.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>

#define pi 3.14159265358979323846
#define MAXC 32 // contacts recorded per cell for bond turnover

// RNG (splitmix64), used only for the initial condition since T = 0
static uint64_t rng_state;
static double rand_unif(void){
	uint64_t z = (rng_state += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	z ^= z >> 31;
	return (z >> 11) * 0x1.0p-53;
}

// Cells, forces and the Verlet list
typedef struct {
	int N;
	double L;            // box side (square, periodic)
	double lambda;       // shell radius / core radius
	double skin;         // Verlet skin
	double *x, *y;       // unwrapped positions
	double *R;           // core radii
	double *fx, *fy;     // total force (= velocity, drag 1)
	double *mx, *my;     // summed contact unit vectors, pointing i -> j
	double *ex, *ey;     // crawling direction (0 if no contacts)
	int *nc;             // contacts per cell
	double *x0, *y0;     // positions at the last list build

	// Verlet list (half list), built with a cell list. Each pair stores its periodic image
	// shift and sigma_ij.
	int npairs, cap;
	int *pair_i, *pair_j;
	double *pair_sx, *pair_sy, *pair_sig;

	// Cell list
	int ncell;
	int *head, *next;

	int nbuilds;
} sys_t;

// Contact count and pair virial sums, sum over contacts of r_ij,a f_ij,b
typedef struct {
	int ncontacts;
	double Wxx, Wyy, Wxy;
} obs_t;

// Contacts, for bond turnover: pair (a < b) stored in row a
typedef struct {
	int *n;
	int *list;
} bonds_t;

static void *xmalloc(size_t n){
	void *p = malloc(n);
	if (!p){
		fprintf(stderr, "out of memory\n");
		exit(1);
	}
	return p;
}

// Add pair i, j to the Verlet list if within the shell plus skin
static void add_pair(sys_t *s, int i, int j){
	double dx = s->x[i] - s->x[j];
	double dy = s->y[i] - s->y[j];
	double sx = s->L * nearbyint(dx / s->L);
	double sy = s->L * nearbyint(dy / s->L);
	dx -= sx;
	dy -= sy;
	double sigma = s->R[i] + s->R[j];
	double rl = s->lambda*sigma + s->skin;
	if (dx*dx + dy*dy < rl*rl){
		if (s->npairs == s->cap){
			s->cap *= 2;
			s->pair_i = realloc(s->pair_i, s->cap * sizeof(int));
			s->pair_j = realloc(s->pair_j, s->cap * sizeof(int));
			s->pair_sx = realloc(s->pair_sx, s->cap * sizeof(double));
			s->pair_sy = realloc(s->pair_sy, s->cap * sizeof(double));
			s->pair_sig = realloc(s->pair_sig, s->cap * sizeof(double));
		}
		s->pair_i[s->npairs] = i;
		s->pair_j[s->npairs] = j;
		s->pair_sx[s->npairs] = sx;
		s->pair_sy[s->npairs] = sy;
		s->pair_sig[s->npairs] = sigma;
		s->npairs += 1;
	}
}

// Rebuild the Verlet list: a cell list when the box holds >= 3 cells per side, otherwise
// all pairs
static void build_list(sys_t *s){
	int i, j, c, cx, cy, k, N = s->N, nc = s->ncell;
	double L = s->L;
	s->npairs = 0;

	if (nc < 3){
		for (i=0; i<N; i++){
			for (j=i+1; j<N; j++){
				add_pair(s, i, j);
			}
		}
	}
	else{
		// Sort cells into the grid
		for (c=0; c<nc*nc; c++){
			s->head[c] = -1;
		}
		for (i=0; i<N; i++){
			double xw = s->x[i] - L*floor(s->x[i]/L);
			double yw = s->y[i] - L*floor(s->y[i]/L);
			cx = (int)(xw / L * nc);
			cy = (int)(yw / L * nc);
			if (cx >= nc) cx = nc-1;
			if (cy >= nc) cy = nc-1;
			c = cx + nc*cy;
			s->next[i] = s->head[c];
			s->head[c] = i;
		}

		// Half stencil: own grid cell + 4 neighbours, so each pair of grid cells is visited once
		int off[4][2] = {{1,0},{1,1},{0,1},{-1,1}};
		for (cy=0; cy<nc; cy++){
			for (cx=0; cx<nc; cx++){
				c = cx + nc*cy;
				// Own grid cell
				for (i=s->head[c]; i>=0; i=s->next[i]){
					for (j=s->next[i]; j>=0; j=s->next[j]){
						add_pair(s, i, j);
					}
				}
				// Neighbouring grid cells
				for (k=0; k<4; k++){
					int cx2 = (cx + off[k][0] + nc) % nc;
					int cy2 = (cy + off[k][1] + nc) % nc;
					int c2 = cx2 + nc*cy2;
					for (i=s->head[c]; i>=0; i=s->next[i]){
						for (j=s->head[c2]; j>=0; j=s->next[j]){
							add_pair(s, i, j);
						}
					}
				}
			}
		}
	}

	memcpy(s->x0, s->x, N*sizeof(double));
	memcpy(s->y0, s->y, N*sizeof(double));
	s->nbuilds += 1;
}

// Pair forces and CIL-C crawling. attract = 0 gives repulsive cores only (the initial
// relaxation). If bonds != NULL, also record the contacts.
static void compute_forces(sys_t *s, double v0, int attract, obs_t *o, bonds_t *bonds){
	int n, i, j, N = s->N;
	double lam = 1.0;
	if (attract) lam = s->lambda;

	memset(s->fx, 0, N*sizeof(double));
	memset(s->fy, 0, N*sizeof(double));
	memset(s->mx, 0, N*sizeof(double));
	memset(s->my, 0, N*sizeof(double));
	memset(s->nc, 0, N*sizeof(int));
	if (bonds) memset(bonds->n, 0, N*sizeof(int));
	o->ncontacts = 0;
	o->Wxx = 0;
	o->Wyy = 0;
	o->Wxy = 0;

	for (n=0; n<s->npairs; n++){
		i = s->pair_i[n];
		j = s->pair_j[n];
		double dx = s->x[i] - s->x[j];
		double dy = s->y[i] - s->y[j];
		dx -= s->pair_sx[n];
		dy -= s->pair_sy[n];
		double sigma = s->pair_sig[n];

		// Contact if within the shell
		double ra = lam*sigma;
		double r2 = dx*dx + dy*dy;
		if (r2 >= ra*ra) continue;

		// Unit vector j -> i
		double r = sqrt(r2);
		double inv_r = 1.0 / r;
		double ux = dx * inv_r;
		double uy = dy * inv_r;

		// Linear spring: > 0 repulsive, < 0 attractive
		double mag = sigma - r;
		s->fx[i] += mag*ux;
		s->fy[i] += mag*uy;
		s->fx[j] -= mag*ux;
		s->fy[j] -= mag*uy;

		// Contact directions for CIL (i -> j is -u)
		s->mx[i] -= ux;
		s->my[i] -= uy;
		s->mx[j] += ux;
		s->my[j] += uy;
		s->nc[i] += 1;
		s->nc[j] += 1;

		// Virial
		o->ncontacts += 1;
		o->Wxx += dx * mag*ux;
		o->Wyy += dy * mag*uy;
		o->Wxy += dx * mag*uy;

		// Record the contact in the row of its lower cell
		if (bonds){
			int a = j;
			int b = i;
			if (i < j){
				a = i;
				b = j;
			}
			if (bonds->n[a] < MAXC){
				bonds->list[a*MAXC + bonds->n[a]] = b;
				bonds->n[a] += 1;
			}
		}
	}

	// CIL-C: crawl opposite the net contact direction with force v0
	for (i=0; i<N; i++){
		s->ex[i] = 0;
		s->ey[i] = 0;
		if (s->nc[i] > 0){
			double m = sqrt(s->mx[i]*s->mx[i] + s->my[i]*s->my[i]);
			// Perfectly balanced contacts give no direction
			if (m > 1e-12){
				s->ex[i] = -s->mx[i] / m;
				s->ey[i] = -s->my[i] / m;
			}
		}
		s->fx[i] += v0 * s->ex[i];
		s->fy[i] += v0 * s->ey[i];
	}
}

// Rates of change of contact lengths, d|r_ij|/dt = (v_i - v_j) . r_hat_ij. They vanish for
// any rigid-body motion, so they measure internal motion only. Call after compute_forces.
static void bond_rates(sys_t *s, double *mean_rate, double *max_rate){
	int n, i, j, count = 0;
	*mean_rate = 0;
	*max_rate = 0;
	for (n=0; n<s->npairs; n++){
		i = s->pair_i[n];
		j = s->pair_j[n];
		double dx = s->x[i] - s->x[j] - s->pair_sx[n];
		double dy = s->y[i] - s->y[j] - s->pair_sy[n];
		double ra = s->lambda*s->pair_sig[n];
		double r2 = dx*dx + dy*dy;
		if (r2 >= ra*ra) continue;
		double rate = fabs(((s->fx[i]-s->fx[j])*dx + (s->fy[i]-s->fy[j])*dy) / sqrt(r2));
		*mean_rate += rate;
		if (rate > *max_rate) *max_rate = rate;
		count += 1;
	}
	if (count > 0) *mean_rate /= count;
}

static int cmp_int(const void *p, const void *q){
	int u = *(const int*)p;
	int v = *(const int*)q;
	if (u < v) return -1;
	if (u > v) return 1;
	return 0;
}

// Sort each cell's contact row
static void sort_bonds(int N, bonds_t *b){
	int a;
	for (a=0; a<N; a++){
		qsort(b->list + a*MAXC, b->n[a], sizeof(int), cmp_int);
	}
}

// Count contacts broken (in prev, not cur) and formed (in cur, not prev); both sorted
static void bond_changes(int N, bonds_t *prev, bonds_t *cur, long *broken, long *formed){
	int a;
	*broken = 0;
	*formed = 0;
	for (a=0; a<N; a++){
		int *p = prev->list + a*MAXC;
		int *c = cur->list + a*MAXC;
		int np = prev->n[a], ncur = cur->n[a];
		int ip = 0, ic = 0;
		while (ip < np || ic < ncur){
			if (ic >= ncur || (ip < np && p[ip] < c[ic])){
				*broken += 1;
				ip += 1;
			}
			else if (ip >= np || c[ic] < p[ip]){
				*formed += 1;
				ic += 1;
			}
			else{
				ip += 1;
				ic += 1;
			}
		}
	}
}

// Window measures over [t_w, t], given positions xw, yw and the sorted contacts bw at t_w:
//   f_act  = fraction of cells whose displacement, relative to the mean displacement of
//            their contacts at t_w, exceeds delta (rigid motion of a cluster does not count)
//   q_bond = fraction of contacts at t_w that are gone at t
static void window_measures(sys_t *s, double *xw, double *yw, bonds_t *bw, bonds_t *bcur,
                            double delta, double *f_act, double *q_bond, double *max_rel){
	int i, k, N = s->N;
	double *sx = calloc(N, sizeof(double));
	double *sy = calloc(N, sizeof(double));
	int *cnt = calloc(N, sizeof(int));
	long nb = 0, broken, formed;

	// Summed displacement of each cell's contacts
	for (i=0; i<N; i++){
		for (k=0; k<bw->n[i]; k++){
			int j = bw->list[i*MAXC + k];
			sx[i] += s->x[j] - xw[j];
			sy[i] += s->y[j] - yw[j];
			sx[j] += s->x[i] - xw[i];
			sy[j] += s->y[i] - yw[i];
			cnt[i] += 1;
			cnt[j] += 1;
			nb += 1;
		}
	}

	// Displacement relative to the contacts
	int n_act = 0;
	*max_rel = 0;
	for (i=0; i<N; i++){
		double rx = s->x[i] - xw[i];
		double ry = s->y[i] - yw[i];
		if (cnt[i] > 0){
			rx -= sx[i] / cnt[i];
			ry -= sy[i] / cnt[i];
		}
		double rel = sqrt(rx*rx + ry*ry);
		if (rel > delta) n_act += 1;
		if (rel > *max_rel) *max_rel = rel;
	}
	*f_act = (double)n_act / N;

	// Contacts lost
	bond_changes(N, bw, bcur, &broken, &formed);
	*q_bond = 0;
	if (nb > 0) *q_bond = (double)broken / nb;
	free(sx);
	free(sy);
	free(cnt);
}

// SIGTERM (e.g. HTCondor eviction) asks for a checkpoint at the end of the current step
static volatile sig_atomic_t got_sigterm = 0;
static void on_sigterm(int sig){
	(void)sig;
	got_sigterm = 1;
}

// Everything else needed to continue a run exactly. The positions at the last list build
// rebuild the same list; outputs are cut back to the offsets on resume.
typedef struct {
	long next_step;
	int done, rebuild, have_prev, have_win, quiet_windows, nbuilds;
	long off_log, off_traj, off_win;
} ckpt_t;

#define CKPT_MAGIC "CILBDCK1"

static void write_checkpoint(const char *dir, sys_t *s, ckpt_t *c, double *xw, double *yw,
                             bonds_t *bw, bonds_t *bp, unsigned long seed, double dt, double v0){
	char path[1200], tmp[1300];
	int N = s->N;
	snprintf(path, sizeof(path), "%s/checkpoint.bin", dir);
	snprintf(tmp, sizeof(tmp), "%s.tmp", path);

	// Write to a temporary file, then rename, so a crash never leaves half a checkpoint
	FILE *fp = fopen(tmp, "wb");
	if (!fp){
		fprintf(stderr, "could not write %s\n", tmp);
		exit(1);
	}
	fwrite(CKPT_MAGIC, 1, 8, fp);
	fwrite(&N, sizeof(int), 1, fp);
	fwrite(&seed, sizeof(seed), 1, fp);
	fwrite(&dt, sizeof(double), 1, fp);
	fwrite(&v0, sizeof(double), 1, fp);
	fwrite(&s->lambda, sizeof(double), 1, fp);
	fwrite(c, sizeof(ckpt_t), 1, fp);
	fwrite(s->x, sizeof(double), N, fp);
	fwrite(s->y, sizeof(double), N, fp);
	fwrite(s->x0, sizeof(double), N, fp);
	fwrite(s->y0, sizeof(double), N, fp);
	fwrite(xw, sizeof(double), N, fp);
	fwrite(yw, sizeof(double), N, fp);
	fwrite(bw->n, sizeof(int), N, fp);
	fwrite(bw->list, sizeof(int), (size_t)N*MAXC, fp);
	fwrite(bp->n, sizeof(int), N, fp);
	fwrite(bp->list, sizeof(int), (size_t)N*MAXC, fp);
	if (fflush(fp) != 0 || fsync(fileno(fp)) != 0 || fclose(fp) != 0){
		fprintf(stderr, "error writing %s\n", tmp);
		exit(1);
	}
	if (rename(tmp, path) != 0){
		fprintf(stderr, "could not rename %s\n", tmp);
		exit(1);
	}
}

// Read n items; 1 if all were read
static int read_items(FILE *fp, void *ptr, size_t size, size_t n){
	return fread(ptr, size, n, fp) == n;
}

// Returns 1 and fills the state if a checkpoint exists; exits if it belongs to another run
static int read_checkpoint(const char *dir, sys_t *s, ckpt_t *c, double *xw, double *yw,
                           bonds_t *bw, bonds_t *bp, unsigned long seed, double dt, double v0){
	char path[1200], magic[8];
	int N = s->N, n_ck;
	unsigned long seed_ck;
	double dt_ck, v0_ck, lam_ck;
	snprintf(path, sizeof(path), "%s/checkpoint.bin", dir);
	FILE *fp = fopen(path, "rb");
	if (!fp) return 0;

	// Parameters
	int ok = read_items(fp, magic, 1, 8) && memcmp(magic, CKPT_MAGIC, 8) == 0;
	ok = ok && read_items(fp, &n_ck, sizeof(int), 1);
	ok = ok && read_items(fp, &seed_ck, sizeof(seed_ck), 1);
	ok = ok && read_items(fp, &dt_ck, sizeof(double), 1);
	ok = ok && read_items(fp, &v0_ck, sizeof(double), 1);
	ok = ok && read_items(fp, &lam_ck, sizeof(double), 1);
	ok = ok && n_ck == N && seed_ck == seed && dt_ck == dt && v0_ck == v0 && lam_ck == s->lambda;
	if (!ok){
		fprintf(stderr, "%s does not match this run's parameters; delete it to start over\n", path);
		exit(1);
	}

	// State, in the order written
	ok = read_items(fp, c, sizeof(ckpt_t), 1);
	ok = ok && read_items(fp, s->x, sizeof(double), N);
	ok = ok && read_items(fp, s->y, sizeof(double), N);
	ok = ok && read_items(fp, s->x0, sizeof(double), N);
	ok = ok && read_items(fp, s->y0, sizeof(double), N);
	ok = ok && read_items(fp, xw, sizeof(double), N);
	ok = ok && read_items(fp, yw, sizeof(double), N);
	ok = ok && read_items(fp, bw->n, sizeof(int), N);
	ok = ok && read_items(fp, bw->list, sizeof(int), (size_t)N*MAXC);
	ok = ok && read_items(fp, bp->n, sizeof(int), N);
	ok = ok && read_items(fp, bp->list, sizeof(int), (size_t)N*MAXC);
	fclose(fp);
	if (!ok){
		fprintf(stderr, "%s is truncated or corrupt\n", path);
		exit(1);
	}
	return 1;
}

// Wall-clock seconds (a node made slow by a busy host still reports normal CPU time)
static double wall_seconds(void){
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1e-9*ts.tv_nsec;
}

// Number of times this run has left a slow node, kept beside the checkpoint
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
	FILE *fp = fopen(path, "r+");
	int ok = fp != NULL;
	ok = ok && ftruncate(fileno(fp), off) == 0 && fseek(fp, 0, SEEK_END) == 0;
	if (!ok){
		fprintf(stderr, "could not resume %s\n", path);
		exit(1);
	}
	return fp;
}

// Forward Euler step. Returns 1 if a cell has moved more than half the skin since the last
// list build, so the list must be rebuilt.
static int euler_step(sys_t *s, double dt){
	int i, rebuild = 0;
	double lim = 0.25*s->skin*s->skin; // (skin/2)^2
	for (i=0; i<s->N; i++){
		s->x[i] += dt * s->fx[i];
		s->y[i] += dt * s->fy[i];
		double dx = s->x[i] - s->x0[i];
		double dy = s->y[i] - s->y0[i];
		if (dx*dx + dy*dy > lim) rebuild = 1;
	}
	return rebuild;
}

// One traj.txt frame: a line with the box size, time and step, then x, y, ex, ey per cell
static void write_frame(FILE *fp, sys_t *s, double t, long step){
	int i;
	fprintf(fp, "%.10g\t%.10g\t%.10g\t%ld\n", s->L, s->L, t, step);
	for (i=0; i<s->N; i++){
		fprintf(fp, "%.10g\t%.10g\t%.6g\t%.6g\n", s->x[i], s->y[i], s->ex[i], s->ey[i]);
	}
	fflush(fp);
}

// mkdir -p
static void make_dir(const char *path){
	char tmp[1024];
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
		"usage: attractive_cil [options]\n"
		"  --seed S        RNG seed for initial condition (1)\n"
		"  --N N           number of cells (256)\n"
		"  --phi X         core packing fraction (0.4)\n"
		"  --phi_a X       shell packing fraction, sets lambda = sqrt(phi_a/phi) (0.9)\n"
		"  --lambda X      shell radius / core radius, overrides --phi_a\n"
		"  --v0 X          CIL crawling force = free speed (0.1)\n"
		"  --dt X          time step (0.005)\n"
		"  --t_relax X     repulsive-only relaxation time before CIL (50)\n"
		"  --t_run X       run time (1e4)\n"
		"  --log_every n   steps between log lines (200)\n"
		"  --traj_every n  steps between trajectory frames (20000)\n"
		"  --skin X        Verlet skin (0.1)\n"
		"  --tau X         window for f_act, q_bond (100)\n"
		"  --delta X       displacement threshold for f_act, relative to contacts (0.1)\n"
		"  --stop_windows n  stop once f_act = 0 for n consecutive windows (0 = never)\n"
		"  --out DIR       output root directory (data)\n"
		"  --ckpt_minutes X   wall-clock minutes between checkpoints (30; 0 = only at the end)\n"
		"  --ckpt_exit_code C exit with code C after each checkpoint, for HTCondor\n"
		"                     exit-driven self-checkpointing (0 = keep running)\n"
		"  --slow_ns X        leave a slow node: after --slow_check_minutes of wall time, if a step\n"
		"                     takes more than X ns per cell, checkpoint and exit with\n"
		"                     --slow_exit_code so the scheduler can move the job (0 = off)\n"
		"  --slow_check_minutes M   when to measure (5)\n"
		"  --slow_max_moves K       stop moving after K slow exits and finish where it is (3)\n"
		"  --slow_exit_code C       exit code for a slow-node exit (99)\n"
		"A run resumes from <dir>/checkpoint.bin if present, and checkpoints on SIGTERM.\n");
	exit(1);
}

// Flush the outputs and note their sizes in the checkpoint state
static void mark_outputs(ckpt_t *ck, FILE *fp_log, FILE *fp_traj, FILE *fp_win){
	fflush(fp_log);
	fflush(fp_traj);
	fflush(fp_win);
	ck->off_log = ftell(fp_log);
	ck->off_traj = ftell(fp_traj);
	ck->off_win = ftell(fp_win);
}

int main(int argc, char *argv[]){
	int i, k;

	// Default parameters //
	unsigned long seed = 1;
	int N = 256;
	double phi = 0.4, phi_a = 0.9, lambda_in = -1, v0 = 0.1;
	double dt = 0.005, t_relax = 50, t_run = 1e4, skin = 0.1;
	long log_every = 200, traj_every = 20000;
	double tau = 100, delta = 0.1;
	int stop_windows = 0;
	const char *out_root = "data";
	double ckpt_minutes = 30;
	int ckpt_exit_code = 0;
	double slow_ns = 0, slow_check_minutes = 5;
	int slow_max_moves = 3, slow_exit_code = 99;

	// Read arguments //
	for (k=1; k<argc; k++){
		if (k+1 >= argc) usage();
		char *opt = argv[k];
		char *val = argv[k+1];
		k += 1;
		if      (!strcmp(opt, "--seed"))          seed = strtoul(val, NULL, 10);
		else if (!strcmp(opt, "--N"))             N = atoi(val);
		else if (!strcmp(opt, "--phi"))           phi = atof(val);
		else if (!strcmp(opt, "--phi_a"))         phi_a = atof(val);
		else if (!strcmp(opt, "--lambda"))        lambda_in = atof(val);
		else if (!strcmp(opt, "--v0"))            v0 = atof(val);
		else if (!strcmp(opt, "--dt"))            dt = atof(val);
		else if (!strcmp(opt, "--t_relax"))       t_relax = atof(val);
		else if (!strcmp(opt, "--t_run"))         t_run = atof(val);
		else if (!strcmp(opt, "--log_every"))     log_every = atol(val);
		else if (!strcmp(opt, "--traj_every"))    traj_every = atol(val);
		else if (!strcmp(opt, "--skin"))          skin = atof(val);
		else if (!strcmp(opt, "--tau"))           tau = atof(val);
		else if (!strcmp(opt, "--delta"))         delta = atof(val);
		else if (!strcmp(opt, "--stop_windows"))  stop_windows = atoi(val);
		else if (!strcmp(opt, "--out"))           out_root = val;
		else if (!strcmp(opt, "--ckpt_minutes"))  ckpt_minutes = atof(val);
		else if (!strcmp(opt, "--ckpt_exit_code")) ckpt_exit_code = atoi(val);
		else if (!strcmp(opt, "--slow_ns"))       slow_ns = atof(val);
		else if (!strcmp(opt, "--slow_check_minutes")) slow_check_minutes = atof(val);
		else if (!strcmp(opt, "--slow_max_moves")) slow_max_moves = atoi(val);
		else if (!strcmp(opt, "--slow_exit_code")) slow_exit_code = atoi(val);
		else usage();
	}
	rng_state = seed * 0x2545F4914F6CDD1DULL + 12345;

	// Allocate //
	sys_t s;
	memset(&s, 0, sizeof(s));
	s.N = N;
	s.skin = skin;
	s.x = xmalloc(N*sizeof(double));
	s.y = xmalloc(N*sizeof(double));
	s.R = xmalloc(N*sizeof(double));
	s.fx = xmalloc(N*sizeof(double));
	s.fy = xmalloc(N*sizeof(double));
	s.mx = xmalloc(N*sizeof(double));
	s.my = xmalloc(N*sizeof(double));
	s.ex = xmalloc(N*sizeof(double));
	s.ey = xmalloc(N*sizeof(double));
	s.x0 = xmalloc(N*sizeof(double));
	s.y0 = xmalloc(N*sizeof(double));
	s.nc = xmalloc(N*sizeof(int));
	s.next = xmalloc(N*sizeof(int));
	s.cap = 16*N;
	s.pair_i = xmalloc(s.cap*sizeof(int));
	s.pair_j = xmalloc(s.cap*sizeof(int));
	s.pair_sx = xmalloc(s.cap*sizeof(double));
	s.pair_sy = xmalloc(s.cap*sizeof(double));
	s.pair_sig = xmalloc(s.cap*sizeof(double));

	// Radii (diameters uniform in [1, 1.4]) and box //
	double sum_R = 0, sum_R2 = 0, R_max = 0;
	for (i=0; i<N; i++){
		s.R[i] = 0.5*(1.0 + 0.4*rand_unif());
		sum_R += s.R[i];
		sum_R2 += s.R[i]*s.R[i];
		if (s.R[i] > R_max) R_max = s.R[i];
	}
	s.L = sqrt(pi*sum_R2 / phi);

	// Shell size from phi_a: shell area = lambda^2 core area
	if (lambda_in > 0){
		s.lambda = lambda_in;
	}
	else{
		if (phi_a < phi){
			fprintf(stderr, "phi_a must be >= phi\n");
			exit(1);
		}
		s.lambda = sqrt(phi_a/phi);
	}
	double phi_a_true = phi*s.lambda*s.lambda;
	double Fa_mean = (s.lambda - 1.0)*2.0*sum_R/N; // mean largest attractive pair force, (lambda-1)<sigma_ij>

	// Cell list sized for the largest list cutoff
	double rlist_max = 2.0*s.lambda*R_max + skin;
	s.ncell = (int)floor(s.L / rlist_max);
	if (s.ncell >= 3){
		s.head = xmalloc(s.ncell*s.ncell*sizeof(int));
	}

	// Random (overlapping) initial positions
	for (i=0; i<N; i++){
		s.x[i] = s.L*rand_unif();
		s.y[i] = s.L*rand_unif();
	}

	// Output directory, parameters and radii //
	char dir[1024], buf[1200];
	if (lambda_in > 0){
		snprintf(dir, sizeof(dir), "%s/phi_%g_lambda_%g_v0_%g_N_%d_dt_%g/seed_%lu", out_root, phi, s.lambda, v0, N, dt, seed);
	}
	else{
		snprintf(dir, sizeof(dir), "%s/phi_%g_phia_%g_v0_%g_N_%d_dt_%g/seed_%lu", out_root, phi, phi_a, v0, N, dt, seed);
	}
	make_dir(dir);

	snprintf(buf, sizeof(buf), "%s/params.txt", dir);
	FILE *fp = fopen(buf, "w");
	fprintf(fp, "seed %lu\nN %d\nphi %.10g\nphi_a %.10g\nlambda %.10g\nFa_mean %.10g\nv0 %.10g\ndt %.10g\n",
	        seed, N, phi, phi_a_true, s.lambda, Fa_mean, v0, dt);
	fprintf(fp, "t_relax %.10g\nt_run %.10g\nlog_every %ld\ntraj_every %ld\nskin %.10g\nL %.10g\nncell %d\ntau %.10g\ndelta %.10g\n",
	        t_relax, t_run, log_every, traj_every, skin, s.L, s.ncell, tau, delta);
	fclose(fp);

	snprintf(buf, sizeof(buf), "%s/radii.txt", dir);
	fp = fopen(buf, "w");
	for (i=0; i<N; i++){
		fprintf(fp, "%.10g\n", s.R[i]);
	}
	fclose(fp);

	printf("N = %d  L = %.4f  lambda = %.4f  phi_a = %.4f  <F_a,max> = %.4f  v0/<F_a,max> = %.3f  ncell = %d\n",
	       N, s.L, s.lambda, phi_a_true, Fa_mean, v0/Fa_mean, s.ncell);
	printf("output: %s\n", dir);

	// Resume from a checkpoint, if there is one //
	bonds_t bonds_prev, bonds_cur, bonds_win;
	bonds_prev.n = xmalloc(N*sizeof(int));
	bonds_prev.list = xmalloc(N*MAXC*sizeof(int));
	bonds_cur.n = xmalloc(N*sizeof(int));
	bonds_cur.list = xmalloc(N*MAXC*sizeof(int));
	bonds_win.n = xmalloc(N*sizeof(int));
	bonds_win.list = xmalloc(N*MAXC*sizeof(int));
	double *xw = xmalloc(N*sizeof(double));
	double *yw = xmalloc(N*sizeof(double));
	ckpt_t ck;
	memset(&ck, 0, sizeof(ck));
	int resumed = read_checkpoint(dir, &s, &ck, xw, yw, &bonds_win, &bonds_prev, seed, dt, v0);
	if (resumed && ck.done){
		printf("run already complete (%s/checkpoint.bin)\n", dir);
		return 0;
	}
	signal(SIGTERM, on_sigterm);

	obs_t o;
	clock_t t_start = clock();
	long step;
	int rebuild = 0;
	char path_log[1200], path_traj[1200], path_win[1200];
	snprintf(path_log, sizeof(path_log), "%s/log.txt", dir);
	snprintf(path_traj, sizeof(path_traj), "%s/traj.txt", dir);
	snprintf(path_win, sizeof(path_win), "%s/window.txt", dir);
	FILE *fp_log, *fp_traj, *fp_win;

	if (resumed){
		// Rebuild the identical Verlet list from the positions at the last build
		double *xs = xmalloc(N*sizeof(double));
		double *ys = xmalloc(N*sizeof(double));
		memcpy(xs, s.x, N*sizeof(double));
		memcpy(ys, s.y, N*sizeof(double));
		memcpy(s.x, s.x0, N*sizeof(double));
		memcpy(s.y, s.y0, N*sizeof(double));
		build_list(&s);
		memcpy(s.x, xs, N*sizeof(double));
		memcpy(s.y, ys, N*sizeof(double));
		free(xs);
		free(ys);
		s.nbuilds = ck.nbuilds;
		rebuild = ck.rebuild;
		fp_log = reopen_at(path_log, ck.off_log);
		fp_traj = reopen_at(path_traj, ck.off_traj);
		fp_win = reopen_at(path_win, ck.off_win);
		printf("resuming from checkpoint at t = %g\n", ck.next_step*dt);
	}
	else{
		// Repulsive-only relaxation from the random start (no CIL, no attraction)
		long n_relax = (long)(t_relax / dt);
		build_list(&s);
		for (step=0; step<n_relax; step++){
			if (rebuild) build_list(&s);
			compute_forces(&s, 0.0, 0, &o, NULL);
			rebuild = euler_step(&s, dt);
		}
		fp_log = fopen(path_log, "w");
		fprintf(fp_log, "# step\tt\tncontacts\tmean_v\tmax_v\tmean_dv\tmax_dv\tvcm_x\tvcm_y\t"
		                "mean_bond_rate\tmax_bond_rate\tn_isolated\tbroken\tformed\tWxx\tWyy\tWxy\n");
		fp_traj = fopen(path_traj, "w");
		fp_win = fopen(path_win, "w");
		fprintf(fp_win, "# t_start\tt_end\tf_act\tq_bond\tmax_rel_disp   (tau = %g, delta = %g)\n", tau, delta);
	}

	// Windows start on log steps
	long n_run = (long)(t_run / dt);
	long win_every = (long)(tau / dt);
	if (win_every % log_every != 0){
		win_every = (win_every / log_every + 1) * log_every;
	}
	long broken = 0, formed = 0;
	int quiet_windows = ck.quiet_windows;
	int have_prev = ck.have_prev;
	int have_win = ck.have_win;
	time_t t_last_ckpt = time(NULL);

	// Slow-node check: measured once per process, over the first slow_check_minutes
	int slow_moves = read_slow_moves(dir);
	int slow_check = slow_ns > 0 && slow_moves < slow_max_moves;
	double t_wall0 = wall_seconds();
	long step0 = ck.next_step;

	// Run //
	for (step=ck.next_step; step<=n_run; step++){
		int do_log = (step % log_every == 0);

		// Forces, recording the contacts on log steps
		if (rebuild) build_list(&s);
		if (do_log) compute_forces(&s, v0, 1, &o, &bonds_cur);
		else compute_forces(&s, v0, 1, &o, NULL);

		if (do_log){
			// Centre-of-mass velocity: pair forces cancel, so this is the mean CIL force
			double vcm_x = 0, vcm_y = 0;
			for (i=0; i<N; i++){
				vcm_x += s.fx[i];
				vcm_y += s.fy[i];
			}
			vcm_x /= N;
			vcm_y /= N;

			// Speeds, with and without the centre-of-mass motion, and isolated cells
			double mean_v = 0, max_v = 0, mean_dv = 0, max_dv = 0;
			int n_iso = 0;
			for (i=0; i<N; i++){
				double v = sqrt(s.fx[i]*s.fx[i] + s.fy[i]*s.fy[i]);
				double dvx = s.fx[i] - vcm_x;
				double dvy = s.fy[i] - vcm_y;
				double dv = sqrt(dvx*dvx + dvy*dvy);
				mean_v += v;
				mean_dv += dv;
				if (v > max_v) max_v = v;
				if (dv > max_dv) max_dv = dv;
				if (s.nc[i] == 0) n_iso += 1;
			}
			mean_v /= N;
			mean_dv /= N;

			double mean_rate, max_rate;
			bond_rates(&s, &mean_rate, &max_rate);

			// Contacts broken and formed since the last log step
			sort_bonds(N, &bonds_cur);
			if (have_prev){
				bond_changes(N, &bonds_prev, &bonds_cur, &broken, &formed);
			}
			have_prev = 1;

			fprintf(fp_log, "%ld\t%.10g\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t%ld\t%ld\t%.10g\t%.10g\t%.10g\n",
			        step, step*dt, o.ncontacts, mean_v, max_v, mean_dv, max_dv, vcm_x, vcm_y,
			        mean_rate, max_rate, n_iso, broken, formed, o.Wxx, o.Wyy, o.Wxy);

			// Window measures, every tau
			int stop = 0;
			if (step % win_every == 0){
				if (have_win){
					double f_act, q_bond, max_rel;
					window_measures(&s, xw, yw, &bonds_win, &bonds_cur, delta, &f_act, &q_bond, &max_rel);
					fprintf(fp_win, "%.10g\t%.10g\t%.6e\t%.6e\t%.6e\n", (step-win_every)*dt, step*dt, f_act, q_bond, max_rel);
					if (f_act == 0) quiet_windows += 1;
					else quiet_windows = 0;
					if (stop_windows > 0 && quiet_windows >= stop_windows) stop = 1;
				}
				memcpy(xw, s.x, N*sizeof(double));
				memcpy(yw, s.y, N*sizeof(double));
				memcpy(bonds_win.n, bonds_cur.n, N*sizeof(int));
				memcpy(bonds_win.list, bonds_cur.list, N*MAXC*sizeof(int));
				have_win = 1;
			}

			// These contacts become the previous ones
			bonds_t tmp = bonds_prev;
			bonds_prev = bonds_cur;
			bonds_cur = tmp;

			if (step % traj_every == 0 || stop){
				write_frame(fp_traj, &s, step*dt, step);
				printf("t = %-10g  z = %.3f  <|v-vcm|> = %.3e  max d|r_ij|/dt = %.3e  bonds -%ld +%ld\n",
				       step*dt, 2.0*o.ncontacts/N, mean_dv, max_rate, broken, formed);
				fflush(stdout);
				fflush(fp_log);
				fflush(fp_win);
			}

			// Trapped: no cell has moved relative to its neighbours by more than delta in
			// stop_windows consecutive windows of length tau
			if (stop){
				printf("trapped: f_act = 0 for %d windows of tau = %g, stopping at t = %g\n",
				       stop_windows, tau, step*dt);
				break;
			}
		}

		// Move
		rebuild = euler_step(&s, dt);

		// Slow node? Checked every 1000 steps until the measuring window has passed
		int slow_exit = 0;
		if (slow_check && step % 1000 == 0){
			double el = wall_seconds() - t_wall0;
			if (el >= 60*slow_check_minutes){
				double ns = 1e9*el / ((double)(step + 1 - step0) * N);
				printf("speed check: %.0f ns per cell-step (threshold %g)\n", ns, slow_ns);
				if (ns > slow_ns){
					slow_exit = 1;
					write_slow_moves(dir, slow_moves + 1);
					printf("slow node: moving (move %d of %d)\n", slow_moves + 1, slow_max_moves);
				}
				slow_check = 0;
			}
		}

		// Checkpoint on a wall-clock schedule, when asked to stop (SIGTERM), or to leave a
		// slow node
		int due = got_sigterm || slow_exit;
		if (ckpt_minutes > 0 && difftime(time(NULL), t_last_ckpt) >= 60*ckpt_minutes) due = 1;
		if (due && step < n_run){
			ck.next_step = step + 1;
			ck.done = 0;
			ck.rebuild = rebuild;
			ck.have_prev = have_prev;
			ck.have_win = have_win;
			ck.quiet_windows = quiet_windows;
			ck.nbuilds = s.nbuilds;
			mark_outputs(&ck, fp_log, fp_traj, fp_win);
			write_checkpoint(dir, &s, &ck, xw, yw, &bonds_win, &bonds_prev, seed, dt, v0);
			t_last_ckpt = time(NULL);
			if (got_sigterm || slow_exit || ckpt_exit_code > 0){
				fclose(fp_log);
				fclose(fp_traj);
				fclose(fp_win);
				printf("checkpointed at t = %g, exiting\n", (step+1)*dt);
				fflush(stdout);
				if (got_sigterm) return 143;
				if (slow_exit) return slow_exit_code;
				return ckpt_exit_code;
			}
		}
	}

	// Final checkpoint marks the run complete, so a requeued job exits at once
	ck.done = 1;
	ck.next_step = step;
	mark_outputs(&ck, fp_log, fp_traj, fp_win);
	write_checkpoint(dir, &s, &ck, xw, yw, &bonds_win, &bonds_prev, seed, dt, v0);
	fclose(fp_log);
	fclose(fp_traj);
	fclose(fp_win);
	printf("cpu time %.1f s, %d list builds\n", (double)(clock()-t_start)/CLOCKS_PER_SEC, s.nbuilds);
	return 0;
}

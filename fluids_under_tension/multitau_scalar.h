/* multitau_scalar.h
 * On-the-fly multi-tau correlator for a scalar signal
 * (e.g. shear stress autocorrelation for Green-Kubo viscosity)
 *
 * Usage:
 *   #include "multitau_scalar.h"
 *   MTScalarCorrelator sc;
 *   mts_init(&sc, sample_interval, dt);
 *   // in simulation loop:
 *   if (step % sample_interval == 0) mts_push(&sc, sigma_xy);
 *   // after simulation:
 *   mts_dump(&sc, "stress_acf.txt");
 */

#ifndef MULTITAU_SCALAR_H
#define MULTITAU_SCALAR_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef MTS_NLEV
#define MTS_NLEV 40
#endif
#ifndef MTS_MBUF
#define MTS_MBUF 32
#endif

typedef struct {
	int    sample_interval;
	double dt;

	double buf[MTS_NLEV][MTS_MBUF];

	double acf_sum[MTS_NLEV][MTS_MBUF];
	long   acf_cnt[MTS_NLEV][MTS_MBUF];

	int    head[MTS_NLEV];
	int    n_stored[MTS_NLEV];
	long   insert_count[MTS_NLEV];
} MTScalarCorrelator;


void mts_init(MTScalarCorrelator *c, int sample_interval, double md_dt){
	c->sample_interval = sample_interval;
	c->dt = md_dt;

	memset(c->buf,          0, sizeof(c->buf));
	memset(c->acf_sum,      0, sizeof(c->acf_sum));
	memset(c->acf_cnt,      0, sizeof(c->acf_cnt));
	memset(c->head,         0, sizeof(c->head));
	memset(c->n_stored,     0, sizeof(c->n_stored));
	memset(c->insert_count, 0, sizeof(c->insert_count));
}


void mts_push(MTScalarCorrelator *c, double val){
	int lev, j, idx;

	for (lev = 0; lev < MTS_NLEV; lev++){

		int max_lag = c->n_stored[lev];
		if (max_lag > MTS_MBUF - 1)
			max_lag = MTS_MBUF - 1;

		for (j = (lev==0 ? 1 : MTS_MBUF/2); j <= max_lag; j++){
			idx = (c->head[lev] - j + MTS_MBUF) % MTS_MBUF;
			c->acf_sum[lev][j] += val * c->buf[lev][idx];
			c->acf_cnt[lev][j] += 1;
		}

		/* Store new value */
		c->buf[lev][c->head[lev]] = val;
		c->head[lev] = (c->head[lev] + 1) % MTS_MBUF;
		if (c->n_stored[lev] < MTS_MBUF)
			c->n_stored[lev]++;
		c->insert_count[lev]++;

		/* Propagate to next level every 2 insertions */
		if (c->insert_count[lev] % 2 != 0)
			break;
	}
}


void mts_dump(MTScalarCorrelator *c, const char *outfile){
	FILE *fp = fopen(outfile, "w");
	if (!fp){
		fprintf(stderr, "mts_dump: cannot open %s\n", outfile);
		exit(1);
	}

	int lev, j;
	double sample_dt = (double)c->sample_interval * c->dt;

	for (lev = 0; lev < MTS_NLEV; lev++){
		int j_start = (lev == 0) ? 1 : MTS_MBUF / 2;
		int j_max   = c->n_stored[lev];
		if (j_max > MTS_MBUF - 1)
			j_max = MTS_MBUF - 1;

		for (j = j_start; j <= j_max; j++){
			if (c->acf_cnt[lev][j] == 0) continue;
			double lag_time = (double)j * sample_dt * (double)(1L << lev);
			double acf = c->acf_sum[lev][j] / c->acf_cnt[lev][j];
			fprintf(fp, "%.10g\t%.10g\t%ld\n", lag_time, acf,c->acf_cnt[lev][j]);
		}
	}

	if(ferror(fp) || fclose(fp)!=0){ fprintf(stderr,"Correlator write failed\n"); exit(1); }
}

#endif /* MULTITAU_SCALAR_H */
/* multitau_isf.h
 * On-the-fly multi-tau correlator for self-ISF, chi4, and MSD.
 * 2D particles.
 * Self-ISF uses J0(q*dr) for isotropic averaging.
 */

#ifndef MULTITAU_ISF_H
#define MULTITAU_ISF_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef MT_NLEV
#define MT_NLEV 40
#endif
#ifndef MT_MBUF
#define MT_MBUF 32
#endif

typedef struct {
	int    N;
	long   origin_stride;     /* short-lag origins every origin_stride samples */
	double q;
	int    sample_interval;
	double dt;

	double *buf;                 /* coords [NLEV][MBUF][N*2]               */

	/* self-ISF accumulators */
	double isf_sum[MT_NLEV][MT_MBUF];
	double isf_sqsum[MT_NLEV][MT_MBUF];
    double isf_mean[MT_NLEV][MT_MBUF], isf_m2[MT_NLEV][MT_MBUF];
	long   isf_cnt[MT_NLEV][MT_MBUF];

	/* MSD */
	double msd_sum[MT_NLEV][MT_MBUF];
	long   msd_cnt[MT_NLEV][MT_MBUF];

	int    head[MT_NLEV];
	int    n_stored[MT_NLEV];
	long   insert_count[MT_NLEV];
} MTCorrelator;


static inline double *mt_slot(MTCorrelator *mt, int lev, int slot){
	return mt->buf + ((long)lev * MT_MBUF + slot) * mt->N * 2;
}


void mt_init(MTCorrelator *mt, int num_atoms, double q,
             int sample_interval, double md_dt){
	mt->N = num_atoms;
	mt->origin_stride = 1;
	mt->q = q;
	mt->sample_interval = sample_interval;
	mt->dt = md_dt;

	long buf_size = (long)MT_NLEV * MT_MBUF * num_atoms * 2;
	mt->buf = (double *)calloc(buf_size, sizeof(double));

	if (!mt->buf){
		fprintf(stderr, "mt_init: allocation failed\n");
		exit(1);
	}

	memset(mt->isf_sum,   0, sizeof(mt->isf_sum));
	memset(mt->isf_sqsum, 0, sizeof(mt->isf_sqsum));
	memset(mt->isf_mean,0,sizeof(mt->isf_mean));
    memset(mt->isf_m2,0,sizeof(mt->isf_m2));
	memset(mt->isf_cnt,   0, sizeof(mt->isf_cnt));
	memset(mt->msd_sum,   0, sizeof(mt->msd_sum));
	memset(mt->msd_cnt,   0, sizeof(mt->msd_cnt));
	memset(mt->head,         0, sizeof(mt->head));
	memset(mt->n_stored,     0, sizeof(mt->n_stored));
	memset(mt->insert_count, 0, sizeof(mt->insert_count));
}


void mt_free(MTCorrelator *mt){
	free(mt->buf);
	mt->buf = NULL;
}


/* J0(x)=sum_k (-x^2/4)^k/(k!)^2. For x^2<=1, truncating after
 * k=8 has remainder <=(1/4)^9/(9!)^2 < 3e-17 (before roundoff).
 * This avoids sqrt and a libm call for small particle displacements. */
static inline double mt_j0_small(double x2){
    double s=0.25*x2;
    return 1+s*(-1+s*(0.25+s*(-1.0/36+s*(1.0/576+s*(-1.0/14400+
        s*(1.0/518400+s*(-1.0/25401600+s/1625702400)))))));
}

void mt_push(MTCorrelator *mt, double coords[][2]){
	int lev, j, p, idx;
	int N = mt->N;
	double q = mt->q;
    double q2=q*q;
	double Fs, dr2_sum, dx, dy;
	double *stored, *dest;

	for (lev = 0; lev < MT_NLEV; lev++){

		/* 1. Correlate new snapshot against stored snapshots */
		int max_lag = mt->n_stored[lev];
		if (max_lag > MT_MBUF - 1)
			max_lag = MT_MBUF - 1;

		long stride=mt->origin_stride >> lev;
        if(stride<1) stride=1;
        if(mt->insert_count[lev]%stride==0)
        for (j = (lev==0 ? 1 : MT_MBUF/2); j <= max_lag; j++){
			idx = (mt->head[lev] - j + MT_MBUF) % MT_MBUF;
			stored = mt_slot(mt, lev, idx);

			Fs = 0.0;
			dr2_sum = 0.0;
			for (p = 0; p < N; p++){
				dx = coords[p][0] - stored[p * 2];
				dy = coords[p][1] - stored[p * 2 + 1];
                double dr2=dx*dx+dy*dy, x2=q2*dr2;
                Fs += x2<=1 ? mt_j0_small(x2) : j0(q*sqrt(dr2));
				dr2_sum += dr2;
			}
			Fs /= N;
			dr2_sum /= N;

			mt->isf_sum[lev][j]   += Fs;
			mt->isf_sqsum[lev][j] += Fs * Fs;
			mt->isf_cnt[lev][j]   += 1;
            double delta=Fs-mt->isf_mean[lev][j];
            mt->isf_mean[lev][j]+=delta/mt->isf_cnt[lev][j];
            mt->isf_m2[lev][j]+=delta*(Fs-mt->isf_mean[lev][j]);

			mt->msd_sum[lev][j]   += dr2_sum;
			mt->msd_cnt[lev][j]   += 1;
		}

		/* 2. Store new snapshot at head */
		dest = mt_slot(mt, lev, mt->head[lev]);
		for (p = 0; p < N; p++){
			dest[p * 2]     = coords[p][0];
			dest[p * 2 + 1] = coords[p][1];
		}
		mt->head[lev] = (mt->head[lev] + 1) % MT_MBUF;
		if (mt->n_stored[lev] < MT_MBUF)
			mt->n_stored[lev]++;
		mt->insert_count[lev]++;

		/* 3. Propagate to next level every 2 insertions */
		if (mt->insert_count[lev] % 2 != 0)
			break;
	}
}


void mt_dump(MTCorrelator *mt, const char *isf_file, const char *chi4_file,
             const char *msd_file){
	FILE *fp_isf  = fopen(isf_file,  "w");
	FILE *fp_chi4 = fopen(chi4_file, "w");
	FILE *fp_msd  = fopen(msd_file,  "w");
	if (!fp_isf || !fp_chi4 || !fp_msd){
		fprintf(stderr, "mt_dump: cannot open output files\n");
		exit(1);
	}

	int lev, j;
	double sample_dt = (double)mt->sample_interval * mt->dt;

	for (lev = 0; lev < MT_NLEV; lev++){
		int j_start = (lev == 0) ? 1 : MT_MBUF / 2;
		int j_max   = mt->n_stored[lev];
		if (j_max > MT_MBUF - 1)
			j_max = MT_MBUF - 1;

		for (j = j_start; j <= j_max; j++){
			if (mt->isf_cnt[lev][j] == 0) continue;

			double lag_time = (double)j * sample_dt * (double)(1L << lev);
			double avg_isf  = mt->isf_sum[lev][j]  / mt->isf_cnt[lev][j];
			
			double chi4 = mt->isf_cnt[lev][j]<2 ? NAN :
                (double)mt->N*fmax(0.,mt->isf_m2[lev][j]/mt->isf_cnt[lev][j]);
			double avg_msd  = mt->msd_sum[lev][j]  / mt->msd_cnt[lev][j];

			fprintf(fp_isf,  "%.10g\t%.10g\t%ld\n", lag_time, avg_isf,mt->isf_cnt[lev][j]);
			fprintf(fp_chi4, "%.10g\t%.10g\t%ld\n", lag_time, chi4,mt->isf_cnt[lev][j]);
			fprintf(fp_msd,  "%.10g\t%.10g\t%ld\n", lag_time, avg_msd,mt->msd_cnt[lev][j]);
		}
	}

	if(ferror(fp_isf) || fclose(fp_isf)!=0){ fprintf(stderr,"Correlator write failed\n"); exit(1); }
	if(ferror(fp_chi4) || fclose(fp_chi4)!=0){ fprintf(stderr,"Correlator write failed\n"); exit(1); }
	if(ferror(fp_msd) || fclose(fp_msd)!=0){ fprintf(stderr,"Correlator write failed\n"); exit(1); }
}

#endif /* MULTITAU_ISF_H */

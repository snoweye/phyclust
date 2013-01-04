/*
 * PAM := Partitioning Around Medoids
 *
 * original Id: pam.f,v 1.16 2003/06/03 13:40:56 maechler translated by
 * f2c (version 20031025) and run through f2c-clean,v 1.10 2002/03/28
 */

/* copied and modified from pam.c in R's cluster package. I also renamed it
   to partmedoids (both function and filename to avoid conflicts 

   copyright of this software should rightly belong to the R software package 
   folks. I have just modified it to obtain a standalone C program. 

   Ranjan Maitra, Ames, IA 50014, 2008/07/29

   Simplified version: only return class id (ncluv) and center id (med),
   and all other stats are took off. Wei-Chen Chen.

*/

#include <stdio.h>

#include "phyclust_constant.h"
#include "phyclust_pam_ind_2.h"


/* bswap(): the clustering algorithm in 2 parts:  I. build,	II. swap */
void bswap(int kk, int n, int *nrepr, int med_given, int do_swap, int trace_lev,
	   /* nrepr[]: here is boolean (0/1): 1 = "is representative object"  */
	   double *dysma, double *dysmb, double *beter,
	   double *dys, double *sky, double s, double *obj)
{
	int i, j, ij, k,h;
	
	/* Parameter adjustments */
	--nrepr;
	--beter;
	
	--dysma; --dysmb;
	
	if (trace_lev) printf("pam()'s bswap(*, s=%g): ", s);
	
	s = s * 1.1 + 1.;/* larger than all dys[];
			    replacing by DBL_MAX  changes result - why ? */
	
/* IDEA: when n is large compared to k (= kk),
 * ----  rather use a "sparse" representation:
 * instead of boolean vector nrepr[] , use  ind_repr <- which(nrepr) !!
 */
	for (i = 1; i <= n; ++i)
		dysma[i] = s;
	
	if (med_given) {
		if  (trace_lev) printf("medoids given\n");
		
		/* compute dysma[] : dysma[j] = D(j, nearest_representative) */
		for (i = 1; i <= n; ++i) {
			if (nrepr[i] == 1)
				for (j = 1; j <= n; ++j) {
					ij = ind_2(i, j);
					if (dysma[j] > dys[ij])
						dysma[j] = dys[ij];
				}
		}
	}
	else {
		
/*  ====== first algorithm: BUILD. ====== */
		
		if(trace_lev) printf("build %d medoids:\n", kk);
		
		/* find  kk  representatives  aka medoids :  */
		
		for (k = 1; k <= kk; ++k) {
			
			/* compute beter[i] for all non-representatives:
			 * also find ammax := max_{..} and nmax := argmax_i{beter[i]} ... */
			int nmax = -1; /* -Wall */
			double ammax, cmd;
			ammax = 0.;
			for (i = 1; i <= n; ++i) {
				if (nrepr[i] == 0) {
					beter[i] = 0.;
					for (j = 1; j <= n; ++j) {
						cmd = dysma[j] - dys[ind_2(i, j)];
						if (cmd > 0.)
							beter[i] += cmd;
					}
					if (ammax <= beter[i]) {
						/*  does < (instead of <= ) work too? -- NO! */
						ammax = beter[i];
						nmax = i;
					}
				}
			}
			
			nrepr[nmax] = 1;/* = .true. : found new representative */
			if (trace_lev >= 2)
				printf("    new repr. %d\n", nmax);
	    
			/* update dysma[] : dysma[j] = D(j, nearest_representative) */
			for (j = 1; j <= n; ++j) {
				ij = ind_2(nmax, j);
				if (dysma[j] > dys[ij])
					dysma[j] = dys[ij];
			}
		}
		/* output of the above loop:  nrepr[], dysma[], ... */
	}
	
	if(trace_lev) /* >= 2 (?) */ {
		printf("  after build: medoids are");
		for (i = 1; i <= n; ++i)
			if(nrepr[i] == 1) printf(" %2d", i);
		if(trace_lev >= 3) {
			printf("\n  and min.dist dysma[1:n] are\n");
			for (i = 1; i <= n; ++i) {
				printf(" %6.3g", dysma[i]);
				if(i % 10 == 0) printf("\n");
			}
			if(n % 10 != 0) printf("\n");
		} else printf("\n");
	}
	
	*sky = 0.;
	for (j = 1; j <= n; ++j)
		*sky += dysma[j];
	obj[0] = *sky / n;
	
	if (do_swap && (kk > 1 || med_given)) {
		
		double dzsky;
		int hbest = -1, nbest = -1;/* init: -Wall*/
		
/* ====== second algorithm: SWAP. ====== */
		
		/* Hmm: In the following, we RE-compute dysma[];
		 *      don't need it first time; then only need *update* after swap */
		
/*--   Loop : */
	L60:
		for (j = 1; j <= n; ++j) {
			/*  dysma[j] := D_j  d(j, <closest medi>)  [KR p.102, 104]
			 *  dysmb[j] := E_j  d(j, <2-nd cl.medi>)  [p.103] */
			dysma[j] = s;
			dysmb[j] = s;
			for (i = 1; i <= n; ++i) {
				if (nrepr[i]) {
					ij = ind_2(i, j);
					if (dysma[j] > dys[ij]) {
						dysmb[j] = dysma[j];
						dysma[j] = dys[ij];
					} else if (dysmb[j] > dys[ij]) {
						dysmb[j] = dys[ij];
					}
				}
			}
		}
		
		dzsky = 1.; /* 1 is arbitrary > 0; only dzsky < 0 matters in the end */
		for (h = 1; h <= n; ++h) if (!nrepr[h]) {
				for (i = 1; i <= n; ++i) if (nrepr[i]) {
						double dz = 0.;
						/* dz := T_{ih} := sum_j C_{jih}  [p.104] : */
						for (j = 1; j <= n; ++j) { /* if (!nrepr[j]) { */
							int hj = ind_2(h, j);
							ij = ind_2(i, j);
							if (dys[ij] == dysma[j]) {
								double small = dysmb[j] > dys[hj]? dys[hj] : dysmb[j];
								dz += (- dysma[j] + small);
							} else if (dys[hj] < dysma[j]) /* 1c. */
								dz += (- dysma[j] + dys[hj]);
						}
						if (dzsky > dz) {
							dzsky = dz; /* dzsky := min_{i,h} T_{i,h} */
							hbest = h;
							nbest = i;
						}
					}
			}
		if (dzsky < 0.) { /* found an improving swap */
			if(trace_lev >= 2)
				printf( "   swp new %d <-> %d old; decreasing diss. by %g\n",
					hbest, nbest, dzsky);
			nrepr[hbest] = 1;
			nrepr[nbest] = 0;
			*sky += dzsky;
			goto L60;
		}
	}
	obj[1] = *sky / n;
	return;
} /* bswap */


/* -----------------------------------------------------------
 Modified from the pam package.
 cstat(): Compute STATistics (numerical output) concerning each partition
*/
void phyclust_cstat(int kk, int nn, int *nsend, int *nrepr, double *radus,
	   double *s, double *dys, int *ncluv, int *nelem, int *med)
{
    int j, k, ja, jk, nplac, ksmal = -1/* -Wall */;
    double ss = *s * 1.1 + 1.;

    /* Parameter adjustments */
    --med;
    --nelem;
    --ncluv;
    --radus;
    --nrepr;
    --nsend;

    /* nsend[j] := i,  where x[i,] is the medoid to which x[j,] belongs */
    for (j = 1; j <= nn; ++j) {
	if (nrepr[j] == 0) {
	    double dsmal = ss;
	    for (k = 1; k <= nn; ++k) {
		if (nrepr[k] == 1) {
		    int kj_ = ind_2(k, j);
		    if (dsmal > dys[kj_]) {
			dsmal = dys[kj_];
			ksmal = k;
		    }
		}
	    }
	    nsend[j] = ksmal;
	} else {
	    nsend[j] = j;
	}
    }
    /* ncluv[j] := k , the cluster number  (k = 1..kk) */
    jk = 1;
    nplac = nsend[1];
    for (j = 1; j <= nn; ++j) {
	ncluv[j] = 0;
	if (nsend[j] == nplac)
	    ncluv[j] = 1;
    }
    for (ja = 2; ja <= nn; ++ja) {
	nplac = nsend[ja];
	if (ncluv[nplac] == 0) {
	    ++jk;
	    for (j = 2; j <= nn; ++j) {
		if (nsend[j] == nplac)
		    ncluv[j] = jk;
	    }
	    if (jk == kk)
		break;
	}
    }

	for (k = 1; k <= kk; ++k) {
	    int ntt = 0, m = -1/* -Wall */;
	    double ttt = 0.;
	    radus[k] = -1.;
	    for (j = 1; j <= nn; ++j) {
		if (ncluv[j] == k) {
		    double djm;
		    ++ntt;
		    m = nsend[j];
		    nelem[ntt] = j;
		    djm = dys[ind_2(j, m)];
		    ttt += djm;
		    if (radus[k] < djm)
			radus[k] = djm;
		}
	    }
	    if(ntt == 0) printf("bug in C cstat(): ntt=0 !!!\n");
	    med[k] = m;
	}
    return;
} /* End of phyclust_cstat(). */




/*  ### Presume cluster.only = FALSE, do.swap = TRUE, trace.lev = 0.
    res <- .C("pam",
	      as.integer(n),
	      as.integer(jp),					# taken off
	      k,
	      x = x2,						# taken off
	      dys = dv,
	      jdyss = as.integer(diss),				# taken off
	      if(mdata)valmd else double(1),			# taken off
	      if(mdata) jtmd else integer(jp),			# taken off
	      as.integer(ndyst),				# taken off
	      integer(n),		# nsend[]
	      logical(n),		# nrepr[]
	      integer(if(cluster.only) 1 else n), # nelem[]
	      double(n),		# radus[]
	      double(n),		# damer[]
	      avsil = double(n),	# `ttd'
	      double(n),		# separ[]
	      ttsil = as.double(0),
	      obj = as.double(c(cluster.only, trace.lev)),# in & out!
	      med = medID,# in & out(if !cluster.only)
	      clu = integer(n),						# taken off
	      clusinf = if(cluster.only) 0. else matrix(0., k, 5),	# taken off
	      silinf  = if(cluster.only) 0. else matrix(0., n, 4),	# taken off
	      isol = nisol,
	      DUP = FALSE, # care!!
	      PACKAGE = "cluster")
 */
void phyclust_pam(int nn, int kk, double *dys, int *nsend, 
		 int/*logical*/ *nrepr, int *nelem, 
		 double *radus, double *damer, double *avsyl,
		 double *ttsyl, double *obj, int *med, int *ncluv, int *nisol)
{
/*  nn = sample size 
	    ***DEL*** p = dimensionality (not needed, taken off)
    kk = number of clusters
	    ***DEL*** x = dataset (not needed, taken off)
    dys = distance matrix in lower triangular form, double array(1 + nn*(nn-1)/2), the first elment is 0.
	    ***DEL*** valmd, jtmd (not needed, taken off)  
	    ***DEL*** ndyst = specifies metric (not needed, taken off)
    nsend = integer array(nn)
    nrepr = integer array(nn)
    nelem = integer (array(1) if cluster.only, array(nn) o/w)
    radus = double array(nn)
    damer = double array(nn)
    avsyl = double array(nn)
            ***DEL*** separ = double array(nn)
    ttsyl = double array(1), set to 0
    obj = double array(2), set to (0, 0), as.double(c(cluster.only, trace.lev)), # in & out!
    med = center id, integer array(k), set med[0] = 0 for non-initial values.   # in & out(if !cluster.only).
    ncluv = class id, integer array(nn),
	    ***DEL*** clusinf = double (array(1) set to 0 if(cluster.only), matrix(0., k, 5) o/w)
	    ***DEL*** sylinf  = double (array(1) set to 0 if(cluster.only), matrix(0., n, 4) o/w)
    isol = nisol, integer (array(1) if(cluster.only), array(k) o/w), set isol[0] = 1 for swap.
 */
	
	/* Local variables */
	int med_given = (med[0] != 0), /* if true, med[] contain initial medoids */
		do_swap = (nisol[0] != 0);
	
	int k, i, nhalf, trace_lev = (int) obj[1];
	double s, sky;
	
	nhalf = nn * (nn - 1) / 2 + 1; /* nhalf := #{distances}+1 = length(dys) */
	
	/* s := max( dys[.] ), the largest distance */
	for (i = 1, s = 0.; i < nhalf; ++i) /* dys[0] == 0. not used here */
		if (s < dys[i])
			s = dys[i];

	/* FIXME: work with med[] = (i_1, i_2, ..., i_k)
	 * ----- instead nrepr[] = (b_1, ... b_n)   b_i in {0,1} */
	for (i = 0; i < nn; ++i)
		nrepr[i] = 0;
	if(med_given) { /* if true, med[] contain initial medoids */
		
		/* for the moment, translate these to nrepr[] 0/1 :
		 * not assuming that the med[] indices are sorted */
		for (k = 0; k < kk; k++)
			nrepr[med[k] - 1] = 1;
	}
	
/*     Build + Swap [but no build if(med_given); swap only if(do_swap) : */
	bswap(kk, nn, nrepr, med_given, do_swap, trace_lev,
	      radus, damer, avsyl, dys, &sky, s, obj);

	if(trace_lev) printf("end{bswap()}, ");
	phyclust_cstat(kk, nn, nsend, nrepr, radus, &s, dys, ncluv, nelem, med);
	if(trace_lev) printf("end{cstat()}\n");
	return;
} /* End of phyclust_pam(). */


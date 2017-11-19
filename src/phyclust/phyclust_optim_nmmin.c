/* This file contains functions for Nelder-Mead. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "phyclust_optim_nmmin.h"

nm_struct* initialize_nm_struct(int n){
	int i;
	nm_struct *nms;

	nms = (nm_struct*) malloc(sizeof(nm_struct));
	nms->n_param = n;
	nms->Bvec = NULL; 
	nms->X = (double*) malloc(n * sizeof(double));
	nms->Fmin = (double*) malloc(sizeof(double));
	nms->fminfn = NULL;
	nms->fail = (int*) malloc(sizeof(int));
	nms->abstol = 1e-16;
	nms->reltol = 1e-8;
	nms->ex = NULL;
	nms->alpha = 1.0;
	nms->beta = 0.5;
	nms->gamma = 2.0;
	nms->trace = 0;
	nms->fncount = (int*) malloc(sizeof(int));
	nms->maxit = 500;

	for(i = 0; i < n; i++){
		nms->X[i] = 0.0;
	}
	*nms->Fmin = 0.0;
	*nms->fail = 0;
	*nms->fncount = 0;
	return(nms);
} /* End of initialize_nm_struct(). */

void free_nm_struct(nm_struct *nms){
	free(nms->X);
	free(nms->Fmin);
	free(nms->fail);
	free(nms->fncount);
	free(nms);
} /* End of free_nm_struct(). */

int phyclust_optim_nmmin(nm_struct *nms){
	int ret_stop = 0;
	ret_stop = optim_nmmin(nms->n_param, nms->Bvec, nms->X, nms->Fmin, nms->fminfn,
			nms->fail, nms->abstol, nms->reltol, nms->ex,
			nms->alpha, nms->beta, nms->gamma, nms->trace,
			nms->fncount, nms->maxit);
	return(ret_stop);
} /* End of phyclust_optim_nmmin(). */




/*
 * Copy from R 2.7.1:src/main/optim.c
 * With modification to avoid conflict with R.
 * Add a struct and an initial function to control the nmmin().
 */

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2007  the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#define big             1.0e+35   /*a very large number*/

#ifdef Rprintf
#undef printf
#define printf Rprintf
#endif

#ifndef Rboolean
#define Rboolean int
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef R_FINITE
#define R_PosInf DBL_MAX
#define R_NegInf DBL_MIN
static int R_FINITE(double x){
    return ((x != R_PosInf) & (x != R_NegInf));
} /* End of R_FINITE(). */
#endif

#ifndef matrix
static double** matrix(int nrh, int nch){
    int i;
    double **m;

    m = (double**) malloc((nrh + 1) * sizeof(double*));
    if(m == NULL){
        printf("Memory allocation fails!\n");
        exit(1);
    }
    for (i = 0; i <= nrh; i++){
	m[i] = (double*) malloc((nch + 1) * sizeof(double));
        if(m[i] == NULL){
            printf("Memory allocation fails!\n");
            exit(1);
        }
    }
    return(m);
} /* End of matrix(). */
#endif


/* Nelder-Mead */
int optim_nmmin(int n, double *Bvec, double *X, double *Fmin,
	   /* optimfn fminfn, */
	   double (*fminfn)(int, double*, void*),
	   /* Rboolean *fail, */
	   int *fail,
	   double abstol, double intol, void *ex,
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit)
/* int n		= number of parameters
   double *Bvec		= initial values (pointer to 1-D array of n elements)
   double *X		= final value (pointer to 1-D array of n elements)
   double *Fmin		= value at located minimum
   optimfn fminfn	= objective function
                          the same as double (*fminfn)(int, double*, void*)
   int *fail		= 1 if no convergence in maximum number of iterations
   double abstol	= 1e-16
   double intol		= 1e-8
   void *ex		= pointer to external data (perhaps as a struct)
   double alpha		= 1.0 default  reflection factor
   double bet		= 0.5 default contraction factor 
   double gamm		= 2.0 default expansion factor 
   int trace		= 0; default tracing on if 1.
   int *fncount		= number of function evaluations/iterations
   int maxit		= maximum number of iterations
*/
{
    char action[50];
    int C;
    Rboolean calcvert;
    double convtol, f;
    int funcount=0, H, i, j, L=0;
    int n1=0;
    double oldsize;
    double **P;
    double size, step, temp, trystep;
    char tstr[24];
    double VH, VL, VR;

    if (maxit <= 0) {
	*Fmin = fminfn(n, Bvec, ex);
	*fncount = 0;
	/* *fail = FALSE; */
	*fail = 0;
	return(1);
    }
    if (trace) printf("  Nelder-Mead direct search function minimizer\n");
    P = matrix(n, n+1);
    /* *fail = FALSE; */
    *fail = 0;
    f = fminfn(n, Bvec, ex);
    if (!R_FINITE(f)) {
	/* error(_("function cannot be evaluated at initial parameters")); */
	*fail = 1;
	if (trace) printf("function cannot be evaluated at initial parameters\n");
	return(1);
    } else {
	if (trace) printf("function value for initial parameters = %f\n", f);
	funcount = 1;
	convtol = intol * (fabs(f) + intol);
	if (trace) printf("  Scaled convergence tolerance is %g\n", convtol);
	n1 = n + 1;
	C = n + 2;
	P[n1 - 1][0] = f;
	for (i = 0; i < n; i++)
	    P[i][0] = Bvec[i];

	L = 1;
	size = 0.0;

	step = 0.0;
	for (i = 0; i < n; i++) {
	    if (0.1 * fabs(Bvec[i]) > step)
		step = 0.1 * fabs(Bvec[i]);
	}
	if (step == 0.0) step = 0.1;
	if (trace) printf("Stepsize computed as %f\n", step);
	for (j = 2; j <= n1; j++) {
	    strcpy(action, "BUILD          ");
	    for (i = 0; i < n; i++)
		P[i][j - 1] = Bvec[i];

	    trystep = step;
	    while (P[j - 2][j - 1] == Bvec[j - 2]) {
		P[j - 2][j - 1] = Bvec[j - 2] + trystep;
		trystep *= 10;
	    }
	    size += trystep;
	}
	oldsize = size;
	calcvert = TRUE;
	do {
	    if (calcvert) {
		for (j = 0; j < n1; j++) {
		    if (j + 1 != L) {
			for (i = 0; i < n; i++)
			    Bvec[i] = P[i][j];
			f = fminfn(n, Bvec, ex);
			if (!R_FINITE(f)) f = big;
			funcount++;
			P[n1 - 1][j] = f;
		    }
		}
		calcvert = FALSE;
	    }

	    VL = P[n1 - 1][L - 1];
	    VH = VL;
	    H = L;

	    for (j = 1; j <= n1; j++) {
		if (j != L) {
		    f = P[n1 - 1][j - 1];
		    if (f < VL) {
			L = j;
			VL = f;
		    }
		    if (f > VH) {
			H = j;
			VH = f;
		    }
		}
	    }

	    if (VH <= VL + convtol || VL <= abstol) break;

	    sprintf(tstr, "%5d", funcount);
	    if (trace) printf("%s%s %f %f\n", action, tstr, VH, VL);

	    for (i = 0; i < n; i++) {
		temp = -P[i][H - 1];
		for (j = 0; j < n1; j++)
		    temp += P[i][j];
		P[i][C - 1] = temp / n;
	    }
	    for (i = 0; i < n; i++)
		Bvec[i] = (1.0 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
	    f = fminfn(n, Bvec, ex);
	    if (!R_FINITE(f)) f = big;
	    funcount++;
	    strcpy(action, "REFLECTION     ");
	    VR = f;
	    if (VR < VL) {
		P[n1 - 1][C - 1] = f;
		for (i = 0; i < n; i++) {
		    f = gamm * Bvec[i] + (1 - gamm) * P[i][C - 1];
		    P[i][C - 1] = Bvec[i];
		    Bvec[i] = f;
		}
		f = fminfn(n, Bvec, ex);
		if (!R_FINITE(f)) f = big;
		funcount++;
		if (f < VR) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		    strcpy(action, "EXTENSION      ");
		} else {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = P[i][C - 1];
		    P[n1 - 1][H - 1] = VR;
		}
	    } else {
		strcpy(action, "HI-REDUCTION   ");
		if (VR < VH) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = VR;
		    strcpy(action, "LO-REDUCTION   ");
		}

		for (i = 0; i < n; i++)
		    Bvec[i] = (1 - bet) * P[i][H - 1] + bet * P[i][C - 1];
		f = fminfn(n, Bvec, ex);
		if (!R_FINITE(f)) f = big;
		funcount++;

		if (f < P[n1 - 1][H - 1]) {
		    for (i = 0; i < n; i++)
			P[i][H - 1] = Bvec[i];
		    P[n1 - 1][H - 1] = f;
		} else {
		    if (VR >= VH) {
			strcpy(action, "SHRINK         ");
			calcvert = TRUE;
			size = 0.0;
			for (j = 0; j < n1; j++) {
			    if (j + 1 != L) {
				for (i = 0; i < n; i++) {
				    P[i][j] = bet * (P[i][j] - P[i][L - 1])
					+ P[i][L - 1];
				    size += fabs(P[i][j] - P[i][L - 1]);
				}
			    }
			}
			if (size < oldsize) {
			    oldsize = size;
			} else {
			    if (trace) printf("Polytope size measure not decreased in shrink\n");
			    *fail = 10;
			    break;
			}
		    }
		}
	    }

	} while (funcount <= maxit);

    }

    if (trace) {
	printf("Exiting from Nelder Mead minimizer\n");
	printf("    %d function evaluations used\n", funcount);
    }
    *Fmin = P[n1 - 1][L - 1];
    for (i = 0; i < n; i++){
	X[i] = P[i][L - 1];
	free(P[i]);
    }
    free(P[i]);
    free(P);
    if (funcount > maxit) *fail = 1;
    *fncount = funcount;
    return(0);
} /* End of optim_nmmin(). */


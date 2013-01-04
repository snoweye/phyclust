/* This file contains declarations for Nelder-Mead. */

#ifdef __HAVE_R_
	#include "phyclust.h"
#endif	/* End of __HAVE_R_. */


#ifndef __PHYCLUST_OPTIM_NMMIN_
#define __PHYCLUST_OPTIM_NMMIN_

typedef struct _nm_struct	nm_struct;


/* For NM subroutine. */
struct _nm_struct{
	int		n_param;			/* number of parameters. */
	double		*Bvec;				/* initial values (pointer to 1-D array of n elements). */
	double		*X;				/* final value (pointer to 1-D array of n elements). */
	double		*Fmin;				/* value at located minimum. */
	/* optimfn	*fminfn; */			/* objective function. */
	double		(*fminfn)(int, double*, void*);
	int		*fail;				/* 1 if no convergence in maximum number of iterations. */
	double		abstol;				/* absolute tolerance, default 1e-16. */
	double		reltol;				/* relative tolerance, default 1e-8. */
	void		*ex;				/* pointer external data for fminfn(). */
	double		alpha;				/* reflection factor, default 1.0. */
	double		beta;				/* contraction factor, default 0.5. */
	double		gamma;				/* expansion factor, default 2.0. */
	int		trace;				/* tracing, default 0. */
	int		*fncount;			/* number of iterations of function evaluations. */
	int		maxit;				/* maximum number of iterations, default 500. */
};

/* File: phyclust_optim_nmmin.c */
nm_struct* initialize_nm_struct(int n);
void free_nm_struct(nm_struct *nms);
int phyclust_optim_nmmin(nm_struct *nms);




/*
 * Copy from R 2.7.1:src/include/R_ext/Applic.h
 * With modification to avoid conflict with R.
 * Add a struct and an initial function to control the nmmin().
 */

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2007   Robert Gentleman, Ross Ihaka
 *                             and the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 *
 *
 * Application Routines, typically implemented in  ../appl/
 * ----------------------------------------------  ========
 */

/* main/optim.c */
/*
typedef double optimfn(int, double *, void *);
*/


/* Nelder-Mead */
int optim_nmmin(int n, double *Bvec, double *X, double *Fmin,
	/* optimfn fminfn, */
	double (*fminfn)(int, double*, void*),
	/* Rboolean *fail,*/
	int *fail,
	double abstol, double reltol, void *ex,
	double alpha, double bet, double gamm, int trace,
	int *fncount, int maxit);
/* int n		= number of parameters
   double *Bvec		= initial values (pointer to 1-D array of n elements)
   double *X		= final value (pointer to 1-D array of n elements)
   double *Fmin		= value at located minimum
   optimfn fminfn	= objective function
                          the same as double (*fminfn)(int, double*, void*)
   int *fail		= 1 if no convergence in maximum number of iterations
   double abstol	= 1e-16
   double reltol	= 1e-8
   void *ex		= pointer to external data (perhaps as a struct)
   double alpha		= 1.0 default  reflection factor
   double bet		= 0.5 default contraction factor 
   double gamm		= 2.0 default expansion factor 
   int trace		= 0; default tracing on if 1.
   int *fncount		= number of function evaluations/iterations
   int maxit		= maximum number of iterations
*/

#endif	/* End of __PHYCLUST_OPTIM_NMMIN_. */


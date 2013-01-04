/* This file contains declarations for Nelder-Mead. */

#include "phyclust_constant.h"

#ifndef __PHYCLUST_APE_NJ_
#define __PHYCLUST_APE_NJ_

typedef struct _nj_struct	nj_struct;

/* For APE NJ subroutine. */
struct _nj_struct{
	/* CAUTION: D requires to point on an 1D array for ape_nj(). */
	double		*D;				/* evolution distance vector, length=N(N-1)/2. */

/* Dynamical variables. */
	int		N;				/* number of sequences. */
	int		n_edge;				/* number of edges, 2N-3. */
	int		n_internal_edge;		/* number of edges, N-3. */
	int		*edge1;				/* id of from. */
	int		*edge2;				/* id of to. */
	double		*edge_length;			/* length between from and to. */
};

/* File: phyclust_ape_nj.c */
nj_struct* initialize_nj_struct(int n);
void free_nj_struct(nj_struct *njs);
void print_njs(int n, nj_struct *njs);
void phyclust_ape_nj(nj_struct *njs);
int check_njs(nj_struct *njs);

/* With modification to avoid conflict with R and ape.
 * Matrix D: a double lower triangular without diagonal format (LTWOD) and has
 *           dimension N (length  N * (N - 1) / 2).
 * edge1, edge2: int length 2 * N - 3 indicating from and to.
 * edge_length: double length 2 * N - 3 indicating edge lengths.
 */
void ape_nj(double *D, int *N, int *edge1, int *edge2, double *edge_length);

#endif	/* End of __PHYCLUST_APE_NJ_. */


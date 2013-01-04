/* This file contains functions for the neighbor-joining method. */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "phyclust_ape_nj.h"

#define Inf DBL_MAX

nj_struct* initialize_nj_struct(int n){
	int i;
	nj_struct *njs;

	njs = (nj_struct*) malloc(sizeof(nj_struct));
	njs->D = NULL;		/* D requires to point on an 1D array for ape_nj(). */
	njs->N = n;
	njs->n_edge = 2 * n - 3;
	njs->n_internal_edge = n - 3;
	njs->edge1 = (int*) malloc(njs->n_edge * sizeof(int));
	njs->edge2 = (int*) malloc(njs->n_edge * sizeof(int));
	njs->edge_length = (double*) malloc(njs->n_edge * sizeof(double));

	for(i = 0; i < njs->n_edge; i++){
		njs->edge1[i] = 0;
		njs->edge2[i] = 0;
		njs->edge_length[i] = 0.0;
	}
	return(njs);
} /* End of initialize_nj_struct(). */

void free_nj_struct(nj_struct *njs){
	free(njs->edge1);
	free(njs->edge2);
	free(njs->edge_length);
	free(njs);
} /* End of free_nj_struct(). */

void print_njs(int n, nj_struct *njs){
	int i;

	if(n > njs->n_edge){
		n = njs->n_edge;
	}
	printf("id  edge1\tedge2\t  length\n");
	for(i = 0; i < n; i++){
		if(njs->edge_length[i] < 1e+8 && njs->edge_length[i] > -1e+8){
			printf("%2d  %5d\t%5d\t%8.4f\n", i, njs->edge1[i], njs->edge2[i], njs->edge_length[i]);
		} else{
			printf("%2d  %5d\t%5d\t%8.4e\n", i, njs->edge1[i], njs->edge2[i], njs->edge_length[i]);
		}
	}
} /* End of print_njs(). */

void phyclust_ape_nj(nj_struct *njs){
	ape_nj(njs->D, &njs->N, njs->edge1, njs->edge2, njs->edge_length);
	/* Force to 0 if length is less than 0.
	for(i = 0; i < njs->n_edge; i++){
		if(njs->edge_length[i] < 0){
			njs->edge_length[i] = 0.0;
		}
	}
	*/
} /* End of phyclust_ape_nj(). */

int check_njs(nj_struct *njs){
	int i, ret = 1;

	for(i = 0; i < njs->n_edge; i++){
		/*
		if(njs->edge_length[i] < 0 || njs->edge_length[i] > Inf){
		*/
		if(njs->edge_length[i] < -Inf || njs->edge_length[i] > Inf){
			ret = 0;
			break;
		}
	}
	return(ret);
} /* End of check_njs(). */




/* With modification to avoid conflict with R and ape.
 *
 * void ape_nj(double *D, int *N, int *edge1, int *edge2, double *edge_length)
 * Matrix D: a double lower triangular without diagonal format (LTWOD) and has
 *           dimension N (length  N * (N - 1) / 2).
 * edge1, edge2: int length 2 * N - 3 indicating from and to.
 * edge_length: double length 2 * N - 3 indicating edge lengths.
 *
 * Original R function call C.
 *
## nj.R (2009-07-10)
nj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C("nj", as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3),
              DUP = FALSE, NAOK = TRUE, PACKAGE = "ape")
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
 */


/* nj.c       2009-07-17 */

/* Copyright 2006-2009 Emmanuel Paradis. */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

/*
#include <R.h>
*/

#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

int give_index(int i, int j, int n)
{
	if (i > j) return(DINDEX(j, i));
	else return(DINDEX(i, j));
} /* End of give_index(). */

double sum_dist_to_i(int n, double *D, int i)
/* returns the sum of all distances D_ij between i and j
   with j = 1...n and j != i */
{
/* we use the fact that the distances are arranged sequentially
   in the lower triangle, e.g. with n = 6 the 15 distances are
   stored as (the C indices are indicated):

           i
     1  2  3  4  5

  2  0
  3  1  5
j 4  2  6  9
  5  3  7 10 12
  6  4  8 11 13 14

  so that we sum the values of the ith column-1st loop-and those of
  (i - 1)th row (labelled 'i')-2nd loop */

	double sum = 0;
	int j, start, end;

	if (i < n) {
		/* the expression below CANNOT be factorized
		   because of the integer operations (it took
		   me a while to find out...) */
		start = n*(i - 1) - i*(i - 1)/2;
		end = start + n - i;
		for (j = start; j < end; j++) sum += D[j];
	}

	if (i > 1) {
		start = i - 2;
		for (j = 1; j <= i - 1; j++) {
			sum += D[start];
			start += n - j - 1;
		}
	}

	return(sum);
} /* End of sum_dist_to_i(). */

void ape_nj(double *D, int *N, int *edge1, int *edge2, double *edge_length)
{
	double *S, Sdist, Ndist, *new_dist, A, B, smallest_S, *DI, d_i, x, y;
	int n, i, j, k, ij, smallest, OTU1, OTU2, cur_nod, o_l, *otu_label;

	OTU1 = 0;
	OTU2 = 0;
	smallest = 0;

	S = &Sdist;
	new_dist = &Ndist;
	otu_label = &o_l;
	DI = &d_i;

	n = *N;
	cur_nod = 2*n - 2;

	/*
	S = (double*)R_alloc(n, sizeof(double));
	new_dist = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
	otu_label = (int*)R_alloc(n, sizeof(int));
	DI = (double*)R_alloc(n - 2, sizeof(double));
	*/
	S = (double*) malloc(n * sizeof(double));
	new_dist = (double*) malloc(n*(n - 1)/2 * sizeof(double));
	otu_label = (int*) malloc(n * sizeof(int));
	DI = (double*) malloc((n - 2) * sizeof(double));
	if(S == NULL || new_dist == NULL || otu_label == NULL || DI == NULL){
		printf("Memory allocation fails!\n");
		exit(1);
	}

	for (i = 0; i < n; i++) otu_label[i] = i + 1;
	k = 0;

	while (n > 3) {

		for (i = 0; i < n; i++)
			S[i] = sum_dist_to_i(n, D, i + 1);

		ij = 0;
		smallest_S = 1e50;
		B = n - 2;
		for (i = 0; i < n - 1; i++) {
			for (j = i + 1; j < n; j++) {
				A = D[ij] - (S[i] + S[j])/B;
				if (A < smallest_S) {
					OTU1 = i + 1;
					OTU2 = j + 1;
					smallest_S = A;
					smallest = ij;
				}
				ij++;
			}
		}

		edge2[k] = otu_label[OTU1 - 1];
		edge2[k + 1] = otu_label[OTU2 - 1];
		edge1[k] = edge1[k + 1] = cur_nod;

		/* get the distances between all OTUs but the 2 selected ones
		   and the latter:
		   a) get the sum for both
		   b) compute the distances for the new OTU */
		A = B = ij = 0;
		for (i = 1; i <= n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			x = D[give_index(i, OTU1, n)]; /* dist between OTU1 and i */
 			y = D[give_index(i, OTU2, n)]; /* dist between OTU2 and i */
			new_dist[ij] = (x + y)/2;
			A += x;
			B += y;
			ij++;
		}
		/* compute the branch lengths */
		A /= n - 2;
		B /= n - 2;
		edge_length[k] = (D[smallest] + A - B)/2;
		edge_length[k + 1] = (D[smallest] + B - A)/2;
		DI[cur_nod - *N - 1] = D[smallest];

		/* update before the next loop */
		if (OTU1 > OTU2) { /* make sure that OTU1 < OTU2 */
			i = OTU1;
			OTU1 = OTU2;
			OTU2 = i;
		}
		if (OTU1 != 1)
			for (i = OTU1 - 1; i > 0; i--)
				otu_label[i] = otu_label[i - 1];
		if (OTU2 != n)
			for (i = OTU2; i < n; i++)
				otu_label[i - 1] = otu_label[i];
		otu_label[0] = cur_nod;

		for (i = 1; i < n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			for (j = i + 1; j <= n; j++) {
				if (j == OTU1 || j == OTU2) continue;
				new_dist[ij] = D[DINDEX(i, j)];
				ij++;
			}
		}

		n--;
		for (i = 0; i < n*(n - 1)/2; i++) D[i] = new_dist[i];

		cur_nod--;
		k = k + 2;
	}

	for (i = 0; i < 3; i++) {
		edge1[*N*2 - 4 - i] = cur_nod;
		edge2[*N*2 - 4 - i] = otu_label[i];
	}

	edge_length[*N*2 - 4] = (D[0] + D[1] - D[2])/2;
	edge_length[*N*2 - 5] = (D[0] + D[2] - D[1])/2;
	edge_length[*N*2 - 6] = (D[2] + D[1] - D[0])/2;

	for (i = 0; i < *N*2 - 3; i++) {
		if (edge2[i] <= *N) continue;
		/* In case there are zero branch lengths: */
		if (DI[edge2[i] - *N - 1] == 0) continue;
		edge_length[i] -= DI[edge2[i] - *N - 1]/2;
	}

	free(S);
	free(new_dist);
	free(otu_label);
	free(DI);
} /* End of ape_nj(). */


/* This file contains a function R_edist() called by R wraps edist_matrix() in
 * "R/f_edist.r", and this function calls the relative functions
 * edist() in "src/phyclust_qmatrix.c".
 *
 * Writen: Wei-Chen Chen on 2009/08/20. */

#include <R.h>
#include <Rinternals.h>
#include "phyclust/phyclust.h"

void R_edist_matrix(int edist_model, int N_X, int L, int **X, double *ret){
	int i, j, I = N_X - 1, total;
	double (*edist_D)(int, int*, int*) = get_edist_D(edist_model);

	total = 0;
	for(i = 0; i < I; i++){
		for(j = 0; j < (I - i); j++){
			ret[total + j] = edist_D(L, X[i], X[i + j + 1]);
		}
		total += (I - i);
	}
} /* End of R_edist_matrix(). */


/* This function calls edist_matrix() in "src/phyclust_qmatrix.c" and is
 * called by edist() using .Call() in "R/f_edist.r".
 * Input:
 *   R_edist_model: SEXP[1], index of edist model.
 *   R_N_X: SEXP[1], number of sequences.
 *   R_L: SEXP[1], length of sequences.
 *   R_X: SEXP[1], sequences.
 * Output:
 *   ret: SEXP[N_X * (N_X - 1) / 2], an array contains distance. */
SEXP R_phyclust_edist(SEXP R_edist_model, SEXP R_N_X, SEXP R_L, SEXP R_X){
	/* Declare variables for calling C. */
	int *C_edist_model, *C_N_X, *C_L, **C_X;

	/* Declare variables for R's returning. */
	SEXP ret;
	double *ret_ptr;

	/* Declare variables for processing. */
	int i, *tmp_ptr;
  
	/* Set initial values. */
	C_edist_model = INTEGER(R_edist_model);
	C_N_X = INTEGER(R_N_X);
	C_L = INTEGER(R_L);

	/* Assign data. */
	C_X = allocate_int_2D_AP(*C_N_X);
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X; i++){
		C_X[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}

	/* For return. */
	PROTECT(ret = allocVector(REALSXP, *C_N_X * (*C_N_X - 1) / 2));
	ret_ptr = REAL(ret);

	/* Compute. */
	R_edist_matrix(*C_edist_model, *C_N_X, *C_L, C_X, ret_ptr);

	/* Free memory and release protectation. */
	UNPROTECT(1);

	return(ret);
} /* End of SEXP R_phyclust_edist(). */


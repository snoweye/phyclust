/* This file contains a function R_Pt() called by R wraps Pt() in
 * "R/f_Pt.r", and this function calls the relative functions
 * in "src/phyclust_qmatrix.c".
 *
 * Writen: Wei-Chen Chen on 2009/10/04. */

#include <R.h>
#include <Rinternals.h>
#include "phyclust/phyclust.h"

Q_matrix* R_initialize_Q_matrix(int code_type, int substitution_model);
void R_free_Q_matrix(Q_matrix *Q);

/* This function calls functions in "src/phyclust_qmatrix.c" and is
 * called by Pt() using .Call() in "R/f_Pt.r".
 * Input:
 *   R_pi: SEXP[1], equilibirium distribution, pi's.
 *   R_kappa: SEXP[1], selection parameters, kappa.
 *   R_Tt: SEXP[1], total evolved time.
 *   R_substitution_mode: SEXP[1], index for substitution model.
 * Output:
 *   ret: a list contains everythings returned from Update_log_Pt in C. */
SEXP R_phyclust_logPt(SEXP R_pi, SEXP R_kappa, SEXP R_Tt,
		SEXP R_code_type, SEXP R_substitution_model){
	/* Declare variables for calling C. */
	int *C_code_type, *C_substitution_model;
	Q_matrix *Q;

	/* Declare variables for R's returning. */
	SEXP ret, ret_names, Pt, log_Pt, H;
	double *tmp_ptr;

	/* Declare variables for processing. */
	int i, ncode;
  
	/* Set initial values. */
	C_code_type = INTEGER(R_code_type);
	C_substitution_model = INTEGER(R_substitution_model);
	ncode = NCODE[*C_code_type];

	/* Assign data. */
	Q = R_initialize_Q_matrix(*C_code_type, *C_substitution_model);	
	Q->pi = REAL(R_pi);
	Q->kappa = REAL(R_kappa);
	Q->Tt = REAL(R_Tt);

	/* For return. */
	PROTECT(ret = allocVector(VECSXP, 3));
  	PROTECT(ret_names = allocVector(STRSXP, 3));
	PROTECT(Pt = allocVector(REALSXP, ncode * ncode));
	PROTECT(log_Pt = allocVector(REALSXP, ncode * ncode));
	PROTECT(H = allocVector(REALSXP, ncode));

	/* Set the elements and names. */
	SET_VECTOR_ELT(ret, 0, Pt);
	SET_VECTOR_ELT(ret, 1, log_Pt);
	SET_VECTOR_ELT(ret, 2, H);
	SET_STRING_ELT(ret_names, 0, mkChar("Pt"));
	SET_STRING_ELT(ret_names, 1, mkChar("log.Pt"));
	SET_STRING_ELT(ret_names, 2, mkChar("H"));
	setAttrib(ret, R_NamesSymbol, ret_names);

	/* Assign ret. */
	tmp_ptr = REAL(Pt);
	for(i = 0; i < ncode; i++){
		Q->Pt[i] = tmp_ptr;
		tmp_ptr += ncode;
	}
	tmp_ptr = REAL(log_Pt);
	for(i = 0; i < ncode; i++){
		Q->log_Pt[i] = tmp_ptr;
		tmp_ptr += ncode;
	}
	Q->H = REAL(H);

	/* Compute. */
	Q->Update_log_Pt(Q);

	/* Free memory and release protectation. */
	R_free_Q_matrix(Q);
	UNPROTECT(5);

	return(ret);
} /* End of SEXP R_phyclust_edist(). */


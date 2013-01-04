/* This file contains a function R_phyclust_logL() called by R wraps
 * phyclust.logL() in "R/f_phyclust_logL.r", and this function calls
 * the relative functions in "src/phyclust/" and "src/".
 *
 * Writen: Wei-Chen Chen on 2009/11/07. */

#include <R.h>
#include <Rinternals.h>
#include "R_phyclust.h"


/* This function calls LogL_observed() in "src/phyclust/phyclust_em_tool.c"
 * and is called by phyclust.logL() using .Call() in "R/f_phyclust_logL.r".
 * Input:
 *   R_N_X_org: SEXP[1], number of sequences.
 *   R_L: SEXP[1], length of sequences.
 *   R_X: SEXP[1], sequences.
 *   R_K: SEXP[1], number of clusters.
 *   R_Eta: SEXP[1], Eta.
 *   R_Mu: SEXP[1], Mu.
 *   R_vect: SEXP[1], vect contains pi, kappa, and Tt.
 *   R_substitution_model: SEXP[1], substitution model.
 *   R_identifier: SEXP[1], identifier.
 *   R_code_type: SEXP[1], code_type.
 *   R_label: SEXP[1], labels.
 * Output:
 *   ret: logL. */
SEXP R_phyclust_logL(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_K;
	double *C_vect;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *QA;
	em_phyclust_struct *empcs;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	SEXP ret_logL;
	double *tmp_logL;

	/* Declare variables for processing. */
	int i, *tmp_ptr;


	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_K = INTEGER(R_K);
	C_vect = REAL(R_vect);

	/* Assign controler. */
	EMC = initialize_em_control();
	EMC->substitution_model = INTEGER(R_substitution_model)[0];
	EMC->identifier = INTEGER(R_identifier)[0];
	EMC->code_type = INTEGER(R_code_type)[0];
	update_em_control(EMC);

	/* Assign data. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}
	tmp_ptr = INTEGER(R_Mu);
	for(i = 0; i < *C_K; i++){
		pcs->Mu[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}
	pcs->Eta = REAL(R_Eta);
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign empcs. */
	empcs = initialize_em_phyclust_struct(pcs);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);
	QA->Convert_vect_to_Q_matrix_array(C_vect, QA);
	QA->Update_log_Pt(QA);

	/* Assign returns. */
	PROTECT(ret_logL = allocVector(REALSXP, 1));
	tmp_logL = REAL(ret_logL);

	/* Compute. */
	*tmp_logL = EMFP->LogL_observed(empcs, QA);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(QA);
	free_em_phyclust_struct(empcs);

	UNPROTECT(1);
	return(ret_logL);
} /* End of SEXP R_phyclust_logL(). */


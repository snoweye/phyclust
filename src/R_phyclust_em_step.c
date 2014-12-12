/* This file contains a function R_phyclust_*_step() called by R wraps
 * phyclust.*.step() in "R/f_phyclust_em_step.r", and this function calls
 * the relative functions in "src/phyclust/" and "src/".
 *
 * Writen: Wei-Chen Chen on 2009/11/12. */

#include "R_phyclust.h"

/* This function calls functions in "src/phyclust/phyclust_em_step.c"
 * and is called by phyclust.em.step() using .Call()
 * in "R/f_phyclust_em_step.r".
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
 *   R_label: SEXP[1], labeles.
 * Output:
 *   ret: a list contains everythings returned from phyclust in C. */
SEXP R_phyclust_em_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_K;
	double *C_vect;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *new_QA, *org_QA, *tmp_QA = NULL;
	em_phyclust_struct *empcs, *tmp_empcs = NULL;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	EMPTR emptr = allocate_emptr();
	SEXP emobj;
	int C_protect_length;
	
	/* Declare variables for processing. */
	int i, j, *tmp_ptr;
	double *tmp_ptr_double;

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

	/* Assign data, read only. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	emobj = initialize_emptr(emptr, pcs);			/* !! Don't move this. */
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;			/* Assign poiners. */
		tmp_ptr += *C_L;
	}

	/* Assign parameters. Updates are required, so make a copy. */
	tmp_ptr = INTEGER(R_Mu);				/* Read only. */
	for(i = 0; i < *C_K; i++){
		for(j = 0; j < *C_L; j++){
			pcs->Mu[i][j] = *tmp_ptr;		/* Copy from the original input. */
			tmp_ptr++;
		}
	}
	tmp_ptr_double = REAL(R_Eta);				/* Read only. */
	for(i = 0; i < *C_K; i++){
		pcs->Eta[i] = *tmp_ptr_double;			/* Copy from the original input. */
		tmp_ptr_double++;
	}
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign empcs. */
	empcs = initialize_em_phyclust_struct(pcs);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	org_QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);
	org_QA->Convert_vect_to_Q_matrix_array(C_vect, org_QA);	/* Copy from the original input. */
	org_QA->Update_log_Pt(org_QA);
	new_QA = duplicate_Q_matrix_array(org_QA);

	/* Compute. */
	E_step_simple(empcs, new_QA, EMFP);
	M_step_simple(empcs, new_QA, org_QA, EMC, EMFP, tmp_empcs, tmp_QA);
	empcs->logL_observed = EMFP->LogL_observed(empcs, new_QA);
	EMFP->Copy_empcs_to_pcs(empcs, pcs);

	/* For return. */
	copy_all_to_emptr(pcs, new_QA, EMC, emptr);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(new_QA);
	free_Q_matrix_array(org_QA);
	free_em_phyclust_struct(empcs);
	C_protect_length = emptr->C_protect_length;
	free(emptr);

	UNPROTECT(C_protect_length);
	return(emobj);
} /* End of SEXP R_phyclust_em_step(). */


/* This function calls E_step_simple() in "src/phyclust/phyclust_em_step.c"
 * and is called by phyclust.e.step() using .Call()
 * in "R/f_phyclust_em_step.r".
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
 *   R_Z_state: SEXP[1], Z_state.
 *   R_label: SEXP[1], labels.
 * Output:
 *   ret: Z_normalized. */
SEXP R_phyclust_e_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_Z_state, SEXP R_label){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_K, *C_Z_state;
	double *C_vect;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *QA;
	em_phyclust_struct *empcs;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	SEXP ret_Z_normalized;

	/* Declare variables for processing. */
	int i, k, *tmp_ptr;
	double *tmp_ptr_double;

	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_K = INTEGER(R_K);
	C_Z_state = INTEGER(R_Z_state);
	C_vect = REAL(R_vect);

	/* Assign controler. */
	EMC = initialize_em_control();
	EMC->substitution_model = INTEGER(R_substitution_model)[0];
	EMC->identifier = INTEGER(R_identifier)[0];
	EMC->code_type = INTEGER(R_code_type)[0];
	update_em_control(EMC);

	/* Assign data, read only. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;			/* Assign pointers. */
		tmp_ptr += *C_L;
	}

	/* Assign parameters. No updates are required, so make a link. */
	tmp_ptr = INTEGER(R_Mu);				/* Read only. */
	for(i = 0; i < *C_K; i++){
		pcs->Mu[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}
	pcs->Eta = REAL(R_Eta);					/* Read only. */
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign empcs. */
	empcs = initialize_em_phyclust_struct(pcs);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);
	QA->Convert_vect_to_Q_matrix_array(C_vect, QA);		/* Copy from the original input. */
	QA->Update_log_Pt(QA);

	/* Assign returns. */
	PROTECT(ret_Z_normalized = allocVector(REALSXP, *C_N_X_org * *C_K));
	tmp_ptr_double = REAL(ret_Z_normalized);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->Z_normalized[i] = tmp_ptr_double;
		tmp_ptr_double += *C_K;
	}

	/* Compute. */
	if(*C_Z_state == 1){
		E_step_simple(empcs, QA, EMFP);
	} else{
		EMFP->Update_Z_modified(empcs, QA);
		/* Manually make a fake copy since Z_modified does not copy
		   back to R. Only Z_normalized did. */
		for(i = 0; i < empcs->N_X; i++){
			for(k = 0; k < empcs->K; k++){
				empcs->Z_normalized[i][k] = empcs->Z_modified[i][k];
			}
		}
	}
	EMFP->Copy_empcs_to_pcs(empcs, pcs);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(QA);
	free_em_phyclust_struct(empcs);

	UNPROTECT(1);
	return(ret_Z_normalized);
} /* End of SEXP R_phyclust_e_step(). */


/* This function calls functions in "src/phyclust/phyclust_em_step.c"
   and is called by phyclust.m.step() using .Call()
   in "R/f_phyclust_em_step.r".
   Input:
     R_N_X_org: SEXP[1], number of sequences.
     R_L: SEXP[1], length of sequences.
     R_X: SEXP[1], sequences.
     R_K: SEXP[1], number of clusters.
     R_vect: SEXP[1], vect contains pi, kappa, and Tt. (posterior)
     R_Z_normalized: SEXP[1], Z_normalized. (posterior)
     R_substitution_model: SEXP[1], substitution model.
     R_identifier: SEXP[1], identifier.
     R_code_type: SEXP[1], code_type.
     R_label: SEXP[1], labels.
   Output:
     ret: a list contains everythings returned from phyclust in C.
*/
SEXP R_phyclust_m_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_vect, SEXP R_Z_normalized,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_K;
	double *C_vect, *C_Z_normalized;
	em_control *EMC;
	phyclust_struct *pcs;
	Q_matrix_array *new_QA, *org_QA, *tmp_QA = NULL;
	em_phyclust_struct *empcs, *tmp_empcs = NULL;
	em_fp *EMFP;

	/* Declare variables for R's returning. */
	EMPTR emptr = allocate_emptr();
	SEXP emobj;
	int C_protect_length;
	
	/* Declare variables for processing. */
	int i, k, *tmp_ptr;

	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_K = INTEGER(R_K);
	C_vect = REAL(R_vect);
	C_Z_normalized = REAL(R_Z_normalized);

	/* Assign controler. */
	EMC = initialize_em_control();
	EMC->substitution_model = INTEGER(R_substitution_model)[0];
	EMC->identifier = INTEGER(R_identifier)[0];
	EMC->code_type = INTEGER(R_code_type)[0];
	update_em_control(EMC);

	/* Assign data, read only. */
	pcs = R_initialize_phyclust_struct(EMC->code_type, *C_N_X_org, *C_L, *C_K);
	emobj = initialize_emptr(emptr, pcs);
	tmp_ptr = INTEGER(R_X);
	for(i = 0; i < *C_N_X_org; i++){
		pcs->X_org[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}
	for(i = 0; i < *C_N_X_org; i++){
		for(k = 0; k < *C_K; k++){
			pcs->Z_normalized[i][k] = *C_Z_normalized++;
		}
	}
	assign_class(pcs);
	assign_Mu_by_class(pcs->N_X_org, pcs->K, pcs->L, pcs->ncode, pcs->gap_index, pcs->class_id, pcs->X_org, pcs->Mu);
	update_phyclust_struct(pcs);

	/* Assign labels. */
	R_update_phyclust_label(pcs, R_label);

	/* Assign empcs. */
	empcs = initialize_em_phyclust_struct(pcs);

	/* Assign function pointers. */
	EMFP = initialize_em_fp(EMC, pcs);

	/* Assign QA. */
	org_QA = initialize_Q_matrix_array(EMC->code_type, *C_K, EMC->substitution_model, EMC->identifier);
	org_QA->Convert_vect_to_Q_matrix_array(C_vect, org_QA);
	org_QA->Update_log_Pt(org_QA);
	new_QA = duplicate_Q_matrix_array(org_QA);

	/* Compute. */
	EMFP->Copy_pcs_to_empcs(pcs, empcs);
	M_step_simple(empcs, new_QA, org_QA, EMC, EMFP, tmp_empcs, tmp_QA);
	empcs->logL_observed = EMFP->LogL_observed(empcs, new_QA);
	EMFP->Copy_empcs_to_pcs(empcs, pcs);

	/* For return. */
	copy_all_to_emptr(pcs, new_QA, EMC, emptr);

	/* Free memory and release protectation. */
	free_em_control(EMC);
	R_free_phyclust_struct(pcs);
	free_em_fp(EMFP);
	free_Q_matrix_array(new_QA);
	free_Q_matrix_array(org_QA);
	free_em_phyclust_struct(empcs);
	C_protect_length = emptr->C_protect_length;
	free(emptr);

	UNPROTECT(C_protect_length);
	return(emobj);
} /* End of SEXP R_phyclust_m_step(). */

/* This file contains updating functions for Q matrix arraies. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_qmatrix.h"
#include "phyclust_qmatrix_array.h"
#include "phyclust_tool.h"


/* Allocate a Q matrix array. */
Q_matrix_array* initialize_Q_matrix_array(int code_type, int K, int substitution_model, int identifier){
	int i, k;
	Q_matrix_array *QA;

	QA = (Q_matrix_array*) malloc(sizeof(Q_matrix_array));
	QA->code_type = code_type;
	QA->ncode = NCODE[code_type];
	QA->K = K;
	QA->identifier = identifier;
	QA->total_n_param = 0;
	QA->substitution_model = substitution_model;
	QA->n_param = 0;
	QA->check_param = 1;
	QA->lower_bound = 1e-16;
	QA->upper_bound = 1 - QA->lower_bound;
	QA->Q = (Q_matrix**) malloc(K * sizeof(Q_matrix*));

	/* Allocate QA as a template. */
	QA->Q[0] = (Q_matrix*) malloc(sizeof(Q_matrix));
	QA->Q[0]->check_param = allocate_int_1D(1);
	QA->Q[0]->Pt = allocate_double_SQ(QA->ncode);
	QA->Q[0]->log_Pt = allocate_double_SQ(QA->ncode);
	QA->Q[0]->H = allocate_double_1D(QA->ncode);
	QA->Q[0]->pi = allocate_double_1D(QA->ncode);
	QA->Q[0]->kappa = allocate_double_1D(1);
	QA->Q[0]->Tt = allocate_double_1D(1);

	/* Assign initialized value. */
	QA->Q[0]->code_type = &QA->code_type;
	QA->Q[0]->ncode = &QA->ncode;
	QA->Q[0]->substitution_model = &QA->substitution_model;
	QA->Q[0]->n_param = &QA->n_param;
	assign_FP_to_Q_matrix(substitution_model, QA->Q[0]);
	QA->Q[0]->lower_bound = &QA->lower_bound;
	QA->Q[0]->upper_bound = &QA->upper_bound;
	for(i = 0; i < QA->ncode; i++){
		QA->Q[0]->pi[i] = 1 / (double) QA->ncode;
	}
	*QA->Q[0]->kappa = 1.0;
	*QA->Q[0]->Tt = 1.0;
	*QA->Q[0]->check_param = 1;

	/* EE has a common storage for Pt, log_Pt, but EV, VE, VV have their own. */
	switch(identifier){
		case EE:	/* Q's are the same, Tt's are the same. */
			QA->total_n_param = *QA->Q[0]->n_param;
			QA->Update_log_Pt = &Update_log_Pt_common;
			QA->Check_param = &Check_param_common;
			QA->Convert_vect_to_Q_matrix_array = &Convert_vect_to_Q_matrix_array_EE;
			QA->Convert_Q_matrix_array_to_vect = &Convert_Q_matrix_array_to_vect_EE;
			QA->Copy_Q_matrix_array = &Copy_Q_matrix_array_EE;
			for(k = 1; k < K; k++){
				QA->Q[k] = repoint_Q_matrix(QA->Q[0]);
			}
			break;
		case EV:	/* Q's are the same, Tt's are different. */
			QA->total_n_param = *QA->Q[0]->n_param + K - 1;
			QA->Update_log_Pt = &Update_log_Pt_split;
			QA->Check_param = &Check_param_split;
			QA->Convert_vect_to_Q_matrix_array = &Convert_vect_to_Q_matrix_array_EV;
			QA->Convert_Q_matrix_array_to_vect = &Convert_Q_matrix_array_to_vect_EV;
			QA->Copy_Q_matrix_array = &Copy_Q_matrix_array_EV;
			for(k = 1; k < K; k++){
				QA->Q[k] = repoint_Q_matrix(QA->Q[0]);
				QA->Q[k]->Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->log_Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->H = allocate_double_1D(QA->ncode);
				QA->Q[k]->Tt = allocate_double_1D(1);
				QA->Q[k]->check_param = allocate_int_1D(1);
				*QA->Q[k]->Tt = 1.0;
				*QA->Q[k]->check_param = 1;
			}
			break;
		case VE:	/* Q's are different, Tt's are the same. */
			QA->total_n_param = (*QA->Q[0]->n_param - 1) * K + 1;
			QA->Update_log_Pt = &Update_log_Pt_split;
			QA->Check_param = &Check_param_split;
			QA->Convert_vect_to_Q_matrix_array = &Convert_vect_to_Q_matrix_array_VE;
			QA->Convert_Q_matrix_array_to_vect = &Convert_Q_matrix_array_to_vect_VE;
			QA->Copy_Q_matrix_array = &Copy_Q_matrix_array_VE;
			for(k = 1; k < K; k++){
				QA->Q[k] = repoint_Q_matrix(QA->Q[0]);
				QA->Q[k]->Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->log_Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->H = allocate_double_1D(QA->ncode);
				QA->Q[k]->pi = allocate_double_1D(QA->ncode);
				QA->Q[k]->kappa = allocate_double_1D(1);
				QA->Q[k]->check_param = allocate_int_1D(1);
				for(i = 0; i < QA->ncode; i++){
					QA->Q[0]->pi[i] = 1 / (double) QA->ncode;
				}
				*QA->Q[k]->kappa = 1.0;
				*QA->Q[k]->check_param = 1;
			}
			break;
		case VV:	/* Q's are different, Tt's are different. */
			QA->total_n_param = *QA->Q[0]->n_param * K;
			QA->Update_log_Pt = &Update_log_Pt_split;
			QA->Check_param = &Check_param_split;
			QA->Convert_vect_to_Q_matrix_array = &Convert_vect_to_Q_matrix_array_VV;
			QA->Convert_Q_matrix_array_to_vect = &Convert_Q_matrix_array_to_vect_VV;
			QA->Copy_Q_matrix_array = &Copy_Q_matrix_array_VV;
			for(k = 1; k < K; k++){
				QA->Q[k] = repoint_Q_matrix(QA->Q[0]);
				QA->Q[k]->Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->log_Pt = allocate_double_SQ(*QA->Q[0]->ncode);
				QA->Q[k]->H = allocate_double_1D(QA->ncode);
				QA->Q[k]->pi = allocate_double_1D(QA->ncode);
				QA->Q[k]->kappa = allocate_double_1D(1);
				QA->Q[k]->Tt = allocate_double_1D(1);
				QA->Q[k]->check_param = allocate_int_1D(1);
				for(i = 0; i < QA->ncode; i++){
					QA->Q[0]->pi[i] = 1 / (double) QA->ncode;
				}
				*QA->Q[k]->kappa = 1.0;
				*QA->Q[k]->Tt = 1.0;
				*QA->Q[k]->check_param = 1;
			}
			break;
		default:
			fprintf_stderr("PE: Identifier is not found.\n");
			exit(1);
	}

	QA->tmp_vect = allocate_double_1D(QA->n_param);

	QA->Update_log_Pt(QA);
	return(QA);
} /* End of inititalize_Q_matrix_array(). */

void free_Q_matrix_array(Q_matrix_array *QA){
	int k, K = QA->K;

	switch(QA->identifier){
		case EE:	/* Q's are the same, Tt's are the same. */
			free_double_RT(QA->ncode, QA->Q[0]->Pt);
			free_double_RT(QA->ncode, QA->Q[0]->log_Pt);
			free(QA->Q[0]->H);
			free(QA->Q[0]->pi);
			free(QA->Q[0]->kappa);
			free(QA->Q[0]->Tt);
			free(QA->Q[0]->check_param);
			free(QA->Q[0]);
			for(k = 1; k < K; k++){
				free(QA->Q[k]);
			}
			break;
		case EV:	/* Q's are the same, Tt's are different. */
			free(QA->Q[0]->pi);
			free(QA->Q[0]->kappa);
			for(k = 0; k < K; k++){
				free_double_RT(QA->ncode, QA->Q[k]->Pt);
				free_double_RT(QA->ncode, QA->Q[k]->log_Pt);
				free(QA->Q[k]->H);
				free(QA->Q[k]->Tt);
				free(QA->Q[k]->check_param);
				free(QA->Q[k]);
			}
			break;
		case VE:	/* Q's are different, Tt's are the same. */
			free(QA->Q[0]->Tt);
			for(k = 0; k < K; k++){
				free_double_RT(QA->ncode, QA->Q[k]->Pt);
				free_double_RT(QA->ncode, QA->Q[k]->log_Pt);
				free(QA->Q[k]->H);
				free(QA->Q[k]->pi);
				free(QA->Q[k]->kappa);
				free(QA->Q[k]->check_param);
				free(QA->Q[k]);
			}
			break;
		case VV:	/* Q's are different, Tt's are different. */
			for(k = 0; k < K; k++){
				free_double_RT(QA->ncode, QA->Q[k]->Pt);
				free_double_RT(QA->ncode, QA->Q[k]->log_Pt);
				free(QA->Q[k]->H);
				free(QA->Q[k]->pi);
				free(QA->Q[k]->kappa);
				free(QA->Q[k]->Tt);
				free(QA->Q[k]->check_param);
				free(QA->Q[k]);
			}
			break;
		default:
			fprintf_stderr("PE: Identifier is not found.\n");
			exit(1);
	}
	free(QA->Q);
	free(QA->tmp_vect);
	free(QA);
} /* End of free_Q_matrix_array_EE(). */

Q_matrix_array* duplicate_Q_matrix_array(Q_matrix_array *org_QA){
	Q_matrix_array *new_QA;

	new_QA = initialize_Q_matrix_array(org_QA->code_type, org_QA->K, org_QA->substitution_model, org_QA->identifier);
	new_QA->total_n_param = org_QA->total_n_param;
	new_QA->n_param = org_QA->n_param;
	new_QA->check_param = org_QA->check_param;
	new_QA->lower_bound = org_QA->lower_bound;
	new_QA->upper_bound = org_QA->upper_bound;
	new_QA->Copy_Q_matrix_array(org_QA, new_QA);
	copy_double_1D(org_QA->n_param, new_QA->tmp_vect, org_QA->tmp_vect);
	return(new_QA);
} /* End of duplicate_Q_matrix_array(). */




/* These functions update log(P(t)). */
void Update_log_Pt_common(Q_matrix_array *QA){		/* For EE. */
	QA->Q[0]->Update_log_Pt(QA->Q[0]);
} /* Update_log_Pt_common(). */

void Update_log_Pt_split(Q_matrix_array *QA){		/* For EV, VE, VV. */
	int k;

	for(k = 0; k < QA->K; k++){
		QA->Q[k]->Update_log_Pt(QA->Q[k]);
	}
} /* Update_log_Pt_split(). */




/* These functions update log(P(t)). */
void Check_param_common(Q_matrix_array *QA){		/* For EE. */
	QA->check_param = *QA->Q[0]->check_param;
} /* Chark_param_common(). */

void Check_param_split(Q_matrix_array *QA){		/* For EV, VE, VV. */
	int k;

	QA->check_param = 1;
	for(k = 0; k < QA->K; k++){
		QA->check_param &= *QA->Q[k]->check_param;
	}
} /* Check_param_split(). */




/* Storage of vect for EE:
 * 	0 to (n_param - 1):				for k = 0, 1, 2,..., K-1.
 * Storage of vect for EV:
 * 	0 to (n_param - 2):				for k = 0, 1, 2,..., K-1.
 *	(n_param - 1) to (n_param - 1 + K):		for Tt[0] to Tt[K-1].
 * Storage of vect for VE:
 * 	0 to (n_param - 2):				for k = 0,
 *	(n_param - 1) to (2*(n_param - 1) - 1):		for k = 1,
 *	(2*(n_param - 1)) to (3*(n_param - 1) - 1):	for k = 2,
 *	... until k = K-1.
 *	((K-1)*(n_param - 1))				for Tt.
 * Storage of vect for VV:
 *	0 to (n_param - 1):				for k = 0,
 *	(n_param) to (2*n_param - 1):			for k = 1,
 *	(2*n_param) to (3*n_param - 1):			for k = 2,
 *	... until k = K-1.
 *
 * All of these can be escaped by initializing vect in QA at very beginning,
 * but that will be not easy for possible extensions, e.g. codon model.
 * I keep in this original approach for easy mantainizes and extensions.
 *
 * These copy actions may cause longer computing time, but the parameters
 * are relative fewer in profile function. Even vect is inititalized in QA
 * at very beginning, it still has few vect needed to copy back and forth,
 * e.g. pi[A] to pi[C] in vect have to be converted to Q->pi. However, this
 * can be ignored by reparameterizing Q and changing the whole computation
 * of log_Pt for Q.
 */

/* These functions convert a vector to parameters. */
void Convert_vect_to_Q_matrix_array_EE(double *vect, Q_matrix_array *QA){
	QA->Q[0]->Convert_vect_to_Q_matrix(vect, QA->Q[0]);
	QA->Check_param(QA);
} /* Convert_vect_to_Q_matrix_array_EE(). */

void Convert_vect_to_Q_matrix_array_EV(double *vect, Q_matrix_array *QA){
	int i, k, tmp_n_param = QA->n_param - 1;
	double *tmp_ptr = vect;

	for(i = 0; i < tmp_n_param; i++){
		QA->tmp_vect[i] = vect[i];
	}
	tmp_ptr += tmp_n_param;
	for(k = 0; k < QA->K; k++){
		QA->tmp_vect[tmp_n_param] = *tmp_ptr;
		tmp_ptr++;
		QA->Q[k]->Convert_vect_to_Q_matrix(QA->tmp_vect, QA->Q[k]);
	}
	QA->Check_param(QA);
} /* Convert_vect_to_Q_matrix_array_EV(). */

void Convert_vect_to_Q_matrix_array_VE(double *vect, Q_matrix_array *QA){
	int i, k, tmp_n_param = QA->n_param - 1;
	double *tmp_ptr = vect;

	QA->tmp_vect[tmp_n_param] = vect[QA->total_n_param - 1];
	for(k = 0; k < QA->K; k++){
		for(i = 0; i < tmp_n_param; i++){
			QA->tmp_vect[i] = tmp_ptr[i];
		}
		tmp_ptr += tmp_n_param;
		QA->Q[k]->Convert_vect_to_Q_matrix(QA->tmp_vect, QA->Q[k]);
	}
	QA->Check_param(QA);
} /* Convert_vect_to_Q_matrix_array_VE(). */

void Convert_vect_to_Q_matrix_array_VV(double *vect, Q_matrix_array *QA){
	int k;
	double *tmp_ptr= vect;

	for(k = 0; k < QA->K; k++){
		QA->Q[k]->Convert_vect_to_Q_matrix(tmp_ptr, QA->Q[k]);
		tmp_ptr += QA->n_param;
	}
	QA->Check_param(QA);
} /* Convert_vect_to_Q_matrix_array_VV(). */




/* These functions convert parameters to a vector. */
void Convert_Q_matrix_array_to_vect_EE(Q_matrix_array *QA, double *vect){
	QA->Q[0]->Convert_Q_matrix_to_vect(QA->Q[0], vect);
} /* Convert_Q_matrix_array_to_vect_EE(). */

void Convert_Q_matrix_array_to_vect_EV(Q_matrix_array *QA, double *vect){
	int i, k, tmp_n_param = QA->n_param - 1;
	double *tmp_ptr = vect;

	QA->Q[0]->Convert_Q_matrix_to_vect(QA->Q[0], QA->tmp_vect);
	for(i = 0; i < tmp_n_param; i++){
		vect[i] = QA->tmp_vect[i];
	}
	tmp_ptr += tmp_n_param;
	*tmp_ptr = QA->tmp_vect[tmp_n_param];
	tmp_ptr++;
	for(k = 1; k < QA->K; k++){
		QA->Q[k]->Convert_Q_matrix_to_vect(QA->Q[k], QA->tmp_vect);
		*tmp_ptr = QA->tmp_vect[tmp_n_param];
		tmp_ptr++;
	}
} /* Convert_Q_matrix_array_to_vect_EV(). */

void Convert_Q_matrix_array_to_vect_VE(Q_matrix_array *QA, double *vect){
	int i, k, tmp_n_param = QA->n_param - 1;
	double *tmp_ptr = vect;

	for(k = 0; k < QA->K; k++){
		QA->Q[k]->Convert_Q_matrix_to_vect(QA->Q[k], QA->tmp_vect);
		for(i = 0; i < tmp_n_param; i++){
			 tmp_ptr[i] = QA->tmp_vect[i];
		}
		tmp_ptr += tmp_n_param;
	}
	vect[QA->total_n_param - 1] = QA->tmp_vect[tmp_n_param];
} /* Convert_Q_matrix_array_to_vect_VE(). */

void Convert_Q_matrix_array_to_vect_VV(Q_matrix_array *QA, double *vect){
	int k;
	double *tmp_ptr= vect;

	for(k = 0; k < QA->K; k++){
		QA->Q[k]->Convert_Q_matrix_to_vect(QA->Q[k], tmp_ptr);
		tmp_ptr += QA->n_param;
	}
} /* Convert_Q_matrix_array_to_vect_VV(). */




/* For copy. */
/* Q's are the same, Tt's are the same. */
void Copy_Q_matrix_array_EE(Q_matrix_array *QA_from, Q_matrix_array *QA_to){
	QA_to->check_param = QA_from->check_param;
	copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[0]->log_Pt, QA_to->Q[0]->log_Pt);
	copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[0]->Pt, QA_to->Q[0]->Pt);
	copy_double_1D(QA_from->ncode, QA_from->Q[0]->H, QA_to->Q[0]->H);
	copy_double_1D(QA_from->ncode, QA_from->Q[0]->pi, QA_to->Q[0]->pi);
	*QA_to->Q[0]->kappa = *QA_from->Q[0]->kappa;
	*QA_to->Q[0]->Tt = *QA_from->Q[0]->Tt;
	*QA_to->Q[0]->check_param = *QA_from->Q[0]->check_param;
} /* End of Copy_Q_matrix_array_EE(). */

/* Q's are the same, Tt's are different. */
void Copy_Q_matrix_array_EV(Q_matrix_array *QA_from, Q_matrix_array *QA_to){
	int k;

	QA_to->check_param = QA_from->check_param;
	copy_double_1D(QA_from->ncode, QA_from->Q[0]->pi, QA_to->Q[0]->pi);
	*QA_to->Q[0]->kappa = *QA_from->Q[0]->kappa;
	for(k = 0; k < QA_from->K; k++){
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->Pt, QA_to->Q[k]->Pt);
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->log_Pt, QA_to->Q[k]->log_Pt);
		copy_double_1D(QA_from->ncode, QA_from->Q[k]->H, QA_to->Q[k]->H);
		*QA_to->Q[k]->Tt = *QA_from->Q[k]->Tt;
		*QA_to->Q[k]->check_param = *QA_from->Q[k]->check_param;
	}
} /* End of Copy_Q_matrix_array_EV(). */

/* Q's are different, Tt's are the same. */
void Copy_Q_matrix_array_VE(Q_matrix_array *QA_from, Q_matrix_array *QA_to){
	int k;

	QA_to->check_param = QA_from->check_param;
	for(k = 0; k < QA_from->K; k++){
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->Pt, QA_to->Q[k]->Pt);
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->log_Pt, QA_to->Q[k]->log_Pt);
		copy_double_1D(QA_from->ncode, QA_from->Q[k]->H, QA_to->Q[k]->H);
		copy_double_1D(QA_from->ncode, QA_from->Q[k]->pi, QA_to->Q[k]->pi);
		*QA_to->Q[k]->kappa = *QA_from->Q[k]->kappa;
		*QA_to->Q[k]->check_param = *QA_from->Q[k]->check_param;
	}
	*QA_to->Q[0]->Tt = *QA_from->Q[0]->Tt;
} /* End of Copy_Q_matrix_array_VE(). */

/* Q's are different, Tt's are different. */
void Copy_Q_matrix_array_VV(Q_matrix_array *QA_from, Q_matrix_array *QA_to){
	int k;

	QA_to->check_param = QA_from->check_param;
	for(k = 0; k < QA_from->K; k++){
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->Pt, QA_to->Q[k]->Pt);
		copy_double_RT(QA_from->ncode, QA_from->ncode, QA_from->Q[k]->log_Pt, QA_to->Q[k]->log_Pt);
		copy_double_1D(QA_from->ncode, QA_from->Q[k]->H, QA_to->Q[k]->H);
		copy_double_1D(QA_from->ncode, QA_from->Q[k]->pi, QA_to->Q[k]->pi);
		*QA_to->Q[k]->kappa = *QA_from->Q[k]->kappa;
		*QA_to->Q[k]->Tt = *QA_from->Q[k]->Tt;
		*QA_to->Q[k]->check_param = *QA_from->Q[k]->check_param;
	}
} /* End of Copy_Q_matrix_array_VV(). */


void reset_Q_matrix_array(Q_matrix_array *QA){
	Q_matrix_array *org_QA = initialize_Q_matrix_array(QA->code_type, QA->K, QA->substitution_model, QA->identifier);
	org_QA->Copy_Q_matrix_array(org_QA, QA);
	free_Q_matrix_array(org_QA);
} /* End of reset_Q_matrix_array(). */




/* ----- For debug. ----- */
void print_QA(Q_matrix_array *QA){
	int k;

	printf("identifier: %s, total_n_param = %d\n", IDENTIFIER[QA->identifier], QA->total_n_param);
	for(k = 0; k < QA->K; k++){
		printf("k = %d\n", k);
		QA->Q[k]->Print_Q_matrix(QA->Q[k]);
	}
} /* End of print_QA(). */


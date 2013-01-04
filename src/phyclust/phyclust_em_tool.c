/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_tool.h"
#include "phyclust_em_tool.h"


/* ----- Independent summary tool. ----- */
/* This is a original version of LogL_observed, and has some numerical problems such as overflow or underflow. */
double LogL_observed_original(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;
	double logL_observed, total_sum, a_Z_modified;

	logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			total_sum += exp(a_Z_modified);
		}
		if(empcs->replication_X[n_X] == 1){
			logL_observed += log(total_sum);
		} else{
			logL_observed += log(total_sum) * empcs->replication_X[n_X];
		}
	}
	
	return(logL_observed);
} /* End of LogL_observed_original(). */

double LogL_observed(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;
	int K = empcs->K, flag_out_range;
	double logL_observed, a_Z_normalized[empcs->K], total_sum, scale_exp;

	logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < K; k++){
			a_Z_normalized[k] = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_normalized[k] += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
		}

		/* This function takes log of component densities "a_z_normalized",
		 * stores a total density in "total_sum",
		 * stores a scale for exponent in "scal_exp", and
		 * turn on "flag_out_range" if exponent is unstable. */
		e_step_with_stable_exp(&K, a_Z_normalized, &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			logL_observed += total_sum;
		} else{
			logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
	
	return(logL_observed);
} /* End of LogL_observed(). */

double LogL_observed_label_semi(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;
	int K = empcs->K, flag_out_range;
	double logL_observed, a_Z_normalized[empcs->K], total_sum, scale_exp;

	logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X_unlabeled; n_X++){
		for(k = 0; k < K; k++){
			a_Z_normalized[k] = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_normalized[k] += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
		}

		e_step_with_stable_exp(&K, a_Z_normalized, &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			logL_observed += total_sum;
		} else{
			logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}

	/* For labeled part, semi_unique stores the k-th cluster which the sequence belongs to.
	 * The logL is equal to the log complete-data likelihood. */
	for(n_X = 0; n_X < empcs->N_X_labeled; n_X++){
		k = empcs->label_semi[n_X];
		total_sum = empcs->log_Eta[k];
		for(s_from = 0; s_from < empcs->ncode; s_from++){
			for(s_to = 0; s_to < empcs->ncode; s_to++){
				total_sum += QA->Q[k]->log_Pt[s_from][s_to] *
					empcs->count_Mu_X[n_X][k][s_from][s_to];
			}
		}

		/* Update logL_observed. */
		if(empcs->replication_X[n_X] == 1){
			logL_observed += total_sum;
		} else{
			logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
	
	return(logL_observed);
} /* End of LogL_observed_label_semi(). */

double LogL_observed_label_general(em_phyclust_struct *empcs, Q_matrix_array *QA){
	/* Need to be implemented with solving a system of linear equations. */
	return(0.0);
} /* End of LogL_observed_label_general(). */




/* log complete-data likelihood, for debuging only. */
double LogL_complete(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){		/* QA_H == QA */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_complete(). */

double LogL_complete_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){	/* QA_H != QA */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				a_Z_modified += QA_H->Q[k]->H[s_from] * empcs->count_Mu_X_gap[n_X][k][s_from];
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_complete_gap(). */




/* log complete-data profile likelihood, used in M-steps. */
double LogL_profile(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){		/* QA_H == QA */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_profile(). */

double LogL_profile_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){	/* QA_H != QA */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				a_Z_modified += QA_H->Q[k]->H[s_from] * empcs->count_Mu_X_gap[n_X][k][s_from];
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_profile_gap(). */




/* ----- Initialization tool. ----- */
/* This function is called by
 * initialize_em_phyclust_struct() in "phyclust_em.c" and
 * update_init_<method>() in "phyclust_init_method.c".
 * "ALWAYS and ONLY" unique sequences (decided by empcs->N_X) are counted.
 * For non-unique data, empcs->replication_X will involve or multiply to the targets,
 * such as computing likelihood, finding Mu. */
void initialize_count_Mu_X_and_gap(em_phyclust_struct *empcs){
	int s_from, s_to, n_X, k, l;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					empcs->count_Mu_X[n_X][k][s_from][s_to] = 0;
				}
				if(empcs->gap_flag){
					empcs->count_Mu_X_gap[n_X][k][s_from] = 0;
				}
			}
			for(l = 0; l < empcs->L; l++){
				if(empcs->X[n_X][l] >= 0 && empcs->X[n_X][l] < empcs->gap_index){
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				} else if(empcs->X[n_X][l] == empcs->gap_index){
					/* For gaps. */
					empcs->count_Mu_X_gap[n_X][k][empcs->Mu[k][l]]++;
				} else{
					/* For missing, skip and do nothing. */
				}
			}
		}
	}
} /* End of initialize_count_Mu_X_and_gap(). */

void reset_Mu_non_seg_site(em_phyclust_struct *empcs){
	int i = 0, k, l;
	int K = empcs->K, L = empcs->L, N_seg_site = empcs->N_seg_site;

	for(l = 0; l < L; l++){
		if(i < N_seg_site && l == empcs->seg_site_id[i]){
			i++;
		} else{
			for(k = 0; k < K; k++){
				empcs->Mu[k][l] = empcs->X[0][l];
			}
		}
	}
} /* End of reset_Mu_non_seg_site(). */




/* ----- Checking tool. ----- */
int is_finite(double x){
	return(x < Inf && x > -Inf);
} /* End of is_finite(). */




/* ----- For copy. ----- */
void copy_EMC(em_control *EMC_from, em_control *EMC_to){
	EMC_to->converge_eps = EMC_from->converge_eps;
	EMC_to->converge_error = EMC_from->converge_error;
	EMC_to->converge_flag = EMC_from->converge_flag;
	EMC_to->converge_iter = EMC_from->converge_iter;
	EMC_to->converge_inner_iter = EMC_from->converge_inner_iter;
	EMC_to->converge_cm_iter = EMC_from->converge_cm_iter;
	EMC_to->update_flag = EMC_from->update_flag;
} /* End of copy_EMC(). */

void reassign_label_pointer(em_phyclust_struct *empcs){
	int n_X, n_X_labeled = 0, n_X_unlabeled = 0, N_X_labeled = empcs->N_X_labeled;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		if(n_X_labeled >= N_X_labeled || empcs->label_index[n_X_labeled] != n_X){
			empcs->X_unlabeled[n_X_unlabeled] = empcs->X[n_X];
			empcs->Z_modified_unlabeled[n_X_unlabeled] = empcs->Z_modified[n_X];
			empcs->Z_normalized_unlabeled[n_X_unlabeled] = empcs->Z_normalized[n_X];
			n_X_unlabeled++;
		} else{
			empcs->X_labeled[n_X_labeled] = empcs->X[n_X];
			empcs->Z_modified_labeled[n_X_labeled] = empcs->Z_modified[n_X];
			empcs->Z_normalized_labeled[n_X_labeled] = empcs->Z_normalized[n_X];
			n_X_labeled++;
		}
	}
} /* End of reassign_label_pointer(). */

void Copy_empcs(em_phyclust_struct *empcs_from, em_phyclust_struct *empcs_to){
	int N_X = empcs_from->N_X, L = empcs_from->L, K = empcs_from->K;

	/* For Em. */
	copy_int_1D(K, empcs_from->n_class, empcs_to->n_class);
	copy_int_RT(K, L, empcs_from->Mu, empcs_to->Mu);
	copy_double_RT(N_X, K, empcs_from->Z_modified, empcs_to->Z_modified);
	copy_double_RT(N_X, K, empcs_from->Z_normalized, empcs_to->Z_normalized);
	copy_double_1D(K, empcs_from->Eta, empcs_to->Eta);
	copy_double_1D(K, empcs_from->log_Eta, empcs_to->log_Eta);
	empcs_to->logL_observed = empcs_from->logL_observed;
	copy_int_RT_4D(empcs_from->N_X, empcs_from->K, empcs_from->ncode, empcs_from->ncode, empcs_from->count_Mu_X,
			empcs_to->count_Mu_X);
	if(empcs_to->gap_flag){
		copy_int_RT_3D(empcs_from->N_X, empcs_from->K, empcs_from->ncode, empcs_from->count_Mu_X_gap,
				empcs_to->count_Mu_X_gap);
	}
} /* End of Copy_empcs(). */


void Copy_empcs_to_pcs(em_phyclust_struct *empcs, phyclust_struct *pcs){
	int n_X_org, n_X, k, L = empcs->L, K = empcs->K;
	double a_entropy;

	copy_int_RT(K, L, empcs->Mu, pcs->Mu);
	copy_double_1D(K, empcs->Eta, pcs->Eta);
	pcs->logL_entropy = empcs->logL_observed;
	for(n_X_org = 0; n_X_org < empcs->N_X_org; n_X_org++){
		n_X = empcs->map_X_org_to_X[n_X_org];
		a_entropy = 0.0;
		for(k = 0; k < K; k++){
			pcs->Z_normalized[n_X_org][k] = empcs->Z_normalized[n_X][k];
			if(empcs->Z_normalized[n_X][k] != 0.0){
				a_entropy += empcs->Z_normalized[n_X][k] * log(empcs->Z_normalized[n_X][k]);
			}
		}
		pcs->logL_entropy += a_entropy * empcs->replication_X[n_X];
	}
	pcs->logL_observed = empcs->logL_observed;
} /* End of Copy_empcs_to_pcs(). */


/* For M-step lonely, or semi-supervised. */
void Copy_pcs_to_empcs(phyclust_struct *pcs, em_phyclust_struct *empcs){
	int n_X_org, n_X, k, K = empcs->K;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		n_X_org = empcs->map_X_to_X_org[n_X];
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = pcs->Z_normalized[n_X_org][k];
		}
	}
} /* End of Copy_pcs_to_empcs(). */

void Copy_pcs_to_empcs_label(phyclust_struct *pcs, em_phyclust_struct *empcs){
	int n_X, n_X_org, n_X_labeled = 0, k, K = empcs->K, N_X_labeled = empcs->N_X_labeled;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		if(n_X_labeled >= N_X_labeled && empcs->label_index[n_X_labeled] != n_X){
			n_X_org = empcs->map_X_to_X_org[n_X];
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = pcs->Z_normalized[n_X_org][k];
			}
		} else{
/* Copy pre-assigned information (supervised) from pcl to empcl.
 * This part should be done by initialize_em_phyclust_label() called by
 * initialize_em_phyclust_struct() when emcps is created.
 * So, there is no need to do again, and copy unlabeled information is sufficient.
			for(k = 0; k < K; k++){
				empcs->Z_normalized_labeled[n_X_labeled][k] = pcs->label->prob_unique[n_X_labeled][k];
			}
 * */
			n_X_labeled++;
		}
	}
} /* End of Copy_pcs_to_empcs_label(). */




/* ----- For debug. ----- */
void print_empcs(em_phyclust_struct *empcs){
	int k;

	printf("em_phyclust_struct:\n");
	if(is_finite(empcs->logL_observed)){
		printf("  logL_obs: %.8f\n", empcs->logL_observed);
	} else{
		printf("  logL_obs: %.4e\n", empcs->logL_observed);
	}
	printf("  Eta:");
	for(k = 0; k < empcs->K; k++){
		printf(" %.8f", empcs->Eta[k]);
	}
	printf("\n");
	printf("  n_class:");
	for(k = 0; k < empcs->K; k++){
		printf(" %d", empcs->n_class[k]);
	}
	printf("\n");
} /* End of print_empcs(). */

void print_EMC(em_control *EMC){
	printf("em_control:\n");
	printf("  code type: %s, em method: %s.\n", CODE_TYPE[EMC->code_type], EM_METHOD[EMC->em_method]);
	printf("  init procedure: %s, method: %s\n", INIT_PROCEDURE[EMC->init_procedure], INIT_METHOD[EMC->init_method]);
	printf("  model substitution: %s, distance: %s\n", SUBSTITUTION_MODEL[EMC->substitution_model],
			EDISTANCE_MODEL[EMC->edist_model]);
	printf("  exhaust iter: %d\n", EMC->exhaust_iter);
	printf("  short iter: %d, eps: %.4e\n", EMC->short_iter, EMC->short_eps);
	printf("  EM iter: %d, eps: %.4e\n", EMC->EM_iter, EMC->EM_eps);
	printf("  CM reltol: %.4e, maxit: %d\n", EMC->cm_reltol, EMC->cm_maxit);
	printf("  NM_Mu_given_QA abstol: %.4e, reltol: %.4e, maxit: %d\n",
		EMC->nm_abstol_Mu_given_QA, EMC->nm_reltol_Mu_given_QA, EMC->nm_maxit_Mu_given_QA);
	printf("  NM_QA_given_Mu abstol: %.4e, reltol: %.4e, maxit: %d\n",
		EMC->nm_abstol_QA_given_Mu, EMC->nm_reltol_QA_given_Mu, EMC->nm_maxit_QA_given_Mu);
	printf("  est_non_seg_site: %d\n", EMC->est_non_seg_site);
	if(EMC->converge_flag < 3){
		printf("iter: %d %d %d, convergence: %d, eps: %.4e.\n",
			EMC->converge_iter, EMC->converge_inner_iter, EMC->converge_cm_iter, EMC->converge_flag,
			EMC->converge_eps);
	} else {
		printf("iter: %d %d %d, convergence: %d,\n",
			EMC->converge_iter, EMC->converge_inner_iter, EMC->converge_cm_iter, EMC->converge_flag);
		printf("  eps: %.4e, error: %.4e\n", EMC->converge_eps, EMC->converge_error);
	}
} /* End of print_EMC(). */

void print_rich_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC){
	print_result(pcs, QA, EMC);
	print_Mu(pcs);
	print_class_id(pcs);
} /* End of print_rich_result(). */

void print_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC){
	int k;

	printf("Phyclust Results:\n");
	printf("code type: %s, em method: %s.\n", CODE_TYPE[EMC->code_type], EM_METHOD[EMC->em_method]);
	printf("init procedure: %s, method: %s.\n",
		       	INIT_PROCEDURE[EMC->init_procedure], INIT_METHOD[EMC->init_method]);
	printf("model substitution: %s, distance: %s.\n",
		       	SUBSTITUTION_MODEL[EMC->substitution_model], EDISTANCE_MODEL[EMC->edist_model]);
        printf("label method: %s.\n", LABEL_METHOD[pcs->label->label_method]);
	if(EMC->converge_flag < 3){
		printf("iter: %d %d %d, convergence: %d, check_param: %d, eps: %.4e.\n",
			EMC->converge_iter, EMC->converge_inner_iter, EMC->converge_cm_iter, EMC->converge_flag,
			QA->check_param, EMC->converge_eps);
	} else {
		printf("iter: %d %d %d, convergence: %d, check_param: %d.\n",
			EMC->converge_iter, EMC->converge_inner_iter, EMC->converge_cm_iter, EMC->converge_flag,
			QA->check_param);
		printf("eps: %.4e, error: %.4e.\n", EMC->converge_eps, EMC->converge_error);
	}
	printf("N_X_org: %d, N_X: %d, L: %d, K: %d, p: %d, N_seg_site: %d.\n",
		       	pcs->N_X_org, pcs->N_X, pcs->L, pcs->K, pcs->n_param + QA->total_n_param,
			pcs->N_seg_site);
	if(is_finite(pcs->logL_observed)){
		printf("logL_obs: %.8f, BIC: %.8f, AIC: %.8f, ICL: %.8f.\n",
			       	pcs->logL_observed, pcs->bic, pcs->aic, pcs->icl);
	} else{
		printf("logL_obs: %.4e, BIC: %.4e, AIC: %.4e, ICL: %.4e.\n",
			       	pcs->logL_observed, pcs->bic, pcs->aic, pcs->icl);
	}
	printf("  Eta:");
	for(k = 0; k < pcs->K; k++){
		printf(" %.8f", pcs->Eta[k]);
	}
	printf(".\n");
	printf("  n_class:");
	for(k = 0; k < pcs->K; k++){
		printf(" %d", pcs->n_class[k]);
	}
	printf(".\n");
	print_QA(QA);

	if(pcs->se_type == SE_YES){
		pcs->SE_P->Print_f_err(pcs->SE_P);
	}
} /* End of print_result(). */

void print_Z_modified(em_phyclust_struct *empcs){
	int n_X, k;

	printf("Z_modified:\n");
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		printf("  %d ", n_X);
		for(k = 0; k < empcs->K; k++){
			printf("%.8e ", empcs->Z_modified[n_X][k]);
		}
		printf("\n");
	}
} /* End of print_Z_modified(). */

void print_Z_normalized(em_phyclust_struct *empcs){
	int n_X, k;

	printf("Z_normalized:\n");
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		printf("  %d ", n_X);
		for(k = 0; k < empcs->K; k++){
			printf("%.8e ", empcs->Z_normalized[n_X][k]);
		}
		printf("\n");
	}
} /* End of print_Z_normalized(). */

void print_Eta(em_phyclust_struct *empcs){
	int k;

	printf("Eta: %.8e\n", 1.0 / (double) empcs->N_X_org);
	for(k = 0; k < empcs->K; k++){
		printf(" %.8e", empcs->Eta[k]);
	}
	printf("\n");
} /* End of print_Eta(). */

void print_empcs_Mu(em_phyclust_struct *empcs){
	int k, l;

	printf("Mu:\n");
	for(k = 0; k < empcs->K; k++){
		printf("    ");
		for(l = 0; l < empcs->L; l++){
		#if PRINT_CODE_TYPE == 0
			if(empcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_CODE[empcs->Mu[k][l]]);
			} else if(empcs->code_type == SNP){
				printf("%c ", SNP_CODE[empcs->Mu[k][l]]);
			}
		#else
			if(empcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_ID[empcs->Mu[k][l]]);
			} else if(empcs->code_type == SNP){
				printf("%c ", SNP_ID[empcs->Mu[k][l]]);
			}
		#endif
		}
		printf("\n");
	}
} /* End of print_empcs_Mu(). */

void print_empcs_Mu_seg_site(em_phyclust_struct *empcs){
	int k, l;

	printf("Mu:\n");
	for(k = 0; k < empcs->K; k++){
		printf("    ");
		for(l = 0; l < empcs->N_seg_site; l++){
		#if PRINT_CODE_TYPE == 0
			if(empcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_CODE[empcs->Mu[k][empcs->seg_site_id[l]]]);
			} else if(empcs->code_type == SNP){
				printf("%c ", SNP_CODE[empcs->Mu[k][empcs->seg_site_id[l]]]);
			}
		#else
			if(empcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_ID[empcs->Mu[k][empcs->seg_site_id[l]]]);
			} else if(empcs->code_type == SNP){
				printf("%c ", SNP_ID[empcs->Mu[k][empcs->seg_site_id[l]]]);
			}
		#endif
		}
		printf("\n");
	}
} /* End of print_empcs_Mu_seg_site(). */

void print_count_Mu_X(em_phyclust_struct *empcs, int n_X, int k){
	int i, j, total = 0;

	printf("n=%d, k=%d:", n_X, k);
	for(i = 0; i < empcs->ncode; i++){
		for(j = 0; j < empcs->ncode; j++){
			printf(" %d", empcs->count_Mu_X[n_X][k][i][j]);
			total += empcs->count_Mu_X[n_X][k][i][j];
		}
		printf(" ");
	}
	printf(" total: %d\n", total);
} /* End of print_count_Mu_X(). */

void print_count_Mu_X_gap(em_phyclust_struct *empcs, int n_X, int k){
	int i, total = 0;

	printf("n=%d, k=%d:", n_X, k);
	for(i = 0; i < empcs->ncode; i++){
		printf(" %d", empcs->count_Mu_X_gap[n_X][k][i]);
		total += empcs->count_Mu_X_gap[n_X][k][i];
		printf(" ");
	}
	printf(" total: %d\n", total);
} /* End of print_count_Mu_X_gap(). */


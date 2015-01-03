/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* #include <errno.h> */
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"
#include "phyclust_tool.h"
#include "phyclust_init_method.h"


/* Gamma[n][k]
   = pi_k * L_k(X_n) / sum_i(pi_i * L_i(X_n))
   = 1 / sum_i(pi_i * L_i(X_n) / (pi_k * L_k))
   = 1 / sum_i(exp(log(pi_i) + log(L_i(X_n)) - log(pi_k) - log(L_k(X_n))))
   This function return stable exponential values with a scale exponential value and flag.
   If flag = 1, scale_exp will be used to adjust results everywhere. Otherwise scale_exp = 0.
   *total_sum is for logL_observed.
*/
void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum, double *scale_exp, int *flag_out_range){
	int k;
	double tmp_exp, max_exp, tmp_exp_K, K_double;

	*total_sum = 0.0;
	*scale_exp = 0.0;
	*flag_out_range = 0;
	max_exp = a_Z_normalized[0];
	for(k = 1; k < *K; k++){
		if(a_Z_normalized[k] > max_exp){
			max_exp = a_Z_normalized[k];
		}
	}

	/* tmp_exp = HUGE_VAL for overflow and 0 for underflow.
	 *   e.g. max_exp is large when parameters are near the boundary and is tiny when too many products.
	 *        errno = ERANGE, only when tmp_exp is too huge (HUGE_VAL).
	 *        errno = 0, when tmp_exp is too small. "Avoid to use ERANGE to check overflow or underflow!"
	 * Scale max_exp by 2 such that close to +0 or -0 until errno is not ERANGE or tmp_exp is not HUGE_VAL. */
/* BUG!
	errno = 0;
	tmp_exp = exp(max_exp);
	if(tmp_exp == HUGE_VAL){
		*flag_out_range = 1;
		*scale_exp = (tmp_exp == HUGE_VAL) ? max_exp : -max_exp;
		do{
			errno = 0;
			*scale_exp *= 0.5;
			tmp_exp = exp(*scale_exp);
		} while(errno == ERANGE);
		*scale_exp = max_exp - *scale_exp;
	}
*/
	tmp_exp = exp(max_exp);
	K_double = (double) *K;
	tmp_exp_K = tmp_exp * K_double;	/* Worst case of *total_sum. */
	if(tmp_exp == HUGE_VAL || tmp_exp == 0.0 || tmp_exp_K == HUGE_VAL){
		*flag_out_range = 1;
		*scale_exp = (tmp_exp == HUGE_VAL) ? max_exp : -max_exp;
		do{
			*scale_exp *= 0.5;
			tmp_exp = exp(*scale_exp);
			tmp_exp_K = tmp_exp * K_double;
		} while(tmp_exp == HUGE_VAL || tmp_exp_K == HUGE_VAL);
		*scale_exp = max_exp - *scale_exp;
		/* The *scale_exp is always greater than 0.
		 * If max_exp > 0 and too large, then shift all to left and computable.
		 *   c = max_exp - *scale_exp > 0, and all should minus c.
		 * If max_exp < 0 and too small, then shift all to right and positive computable.
		 *   c = max_exp - *scale_exp < max_exp < 0, and all should minus c. */
	}
	if(*flag_out_range){
		for(k = 0; k < *K; k++){
			a_Z_normalized[k] -= *scale_exp;
		}
	}

	*total_sum = 0.0;
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = exp(a_Z_normalized[k]);
		*total_sum += a_Z_normalized[k];
	}
	for(k = 0; k < *K; k++){
		a_Z_normalized[k] = a_Z_normalized[k] / *total_sum;
	}
} /* End of e_step_with_stable_exp(); */




/* This is a original version of E-step, and has some numerical problems such as overflow or underflow. */
void E_step_original(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP){
	int n_X, k, K = empcs->K;
	double total_sum;

	EMFP->Update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized. */
		total_sum = 0.0;
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = exp(empcs->Z_modified[n_X][k] + empcs->log_Eta[k]);
			total_sum += empcs->Z_normalized[n_X][k];
		}

		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_normalized[n_X][k] / total_sum;
		}
	}
} /* End of E_step_original(). */

/* E-steps simple verion. */
void E_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	EMFP->Update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_modified[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized[n_X], &total_sum, &scale_exp, &flag_out_range);
	}
} /* End of E_step_simple(). */

/* E-steps with none labels. */
void E_step_logL_observed(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	EMFP->Update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	empcs->logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){	/* Update Z_normalized and logL_observed. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = empcs->Z_modified[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized[n_X], &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			empcs->logL_observed += total_sum;
		} else{
			empcs->logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
} /* End of E_step_logL_observed(). */

/* E-steps with simple labels. */
void E_step_logL_observed_label_semi(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP){
	int n_X, k, K = empcs->K, flag_out_range;
	double total_sum, scale_exp;

	EMFP->Update_Z_modified(empcs, QA);	/* Update empcs->Z_modified (log and unnormalized). */
	empcs->logL_observed = 0.0;

	/* For unlabeled part. */
	for(n_X = 0; n_X < empcs->N_X_unlabeled; n_X++){	/* Update Z_normalized and logL_observed. */
		for(k = 0; k < K; k++){
			empcs->Z_normalized_unlabeled[n_X][k] = empcs->Z_modified_unlabeled[n_X][k] + empcs->log_Eta[k];
		}

		e_step_with_stable_exp(&K, empcs->Z_normalized_unlabeled[n_X], &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			empcs->logL_observed += total_sum;
		} else{
			empcs->logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}

	/* For labeled part, semi_unique stores the k-th cluster which the sequence belongs to.
	 * The logL is equal to the log complete-data likelihood. */
	for(n_X = 0; n_X < empcs->N_X_labeled; n_X++){
		k = empcs->label_semi[n_X];
		total_sum = empcs->Z_modified_labeled[n_X][k] + empcs->log_Eta[k];

		/* Update logL_observed. */
		if(empcs->replication_X[n_X] == 1){
			empcs->logL_observed += total_sum;
		} else{
			empcs->logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
} /* End of E_step_logL_observed_label_semi(). */

/* E-steps with general labels. */
void E_step_logL_observed_label_general(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP){
} /* End of E_step_logL_observed_label_general(). */




/* The original M-step:
 * Given Eta, for each fixed QA and Tt find the Mu to maximize the R, and use NM to search the best QA and Tt.
 * By default, set EMC->update_flag = 1 in Em_step().
 * The maximize_logpL() in "phyclust_logpL.c",
 * update_flag = 0 for update Mu given QA. (for profile logL)
 *             = 1 for update QA given Mu. (for profile logL and ECM/AECM) */
int M_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0;

	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC);	/* Find Eta given Z_normalized. */
	if(ret_stop > 0){
		return(ret_stop);
	}
	EMFP->Maximize_logpL(empcs, QA, QA_H, EMC, EMFP);	/* Find QA, Tt, and Mu. */
	return(ret_stop);
} /* End of M_step_simple(). */


/* The conditional M-step:
 * Update Eta given Z_normalized. Iteratively Update Mu and QA given Eta until the R function is converged.
 * Theoretically, this will converge to the same values of M_step_simple().
 * All of these CM steps require to set EMC->update_flag = 1 in Em_step().
 * The maximize_logpL() in "phyclust_logpL.c",
 * update_flag = 0 for update Mu given QA. (for profile logL)
 *             = 1 for update QA given Mu. (for profile logL and ECM/AECM) */
int M_step_CM(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0, cm_iter = 1, flag = 1;
	double cm_reltol = 0.0, R = 0.0, tmp_R = 0.0;

	EMFP->Copy_empcs(empcs, tmp_empcs);
	tmp_QA->Copy_Q_matrix_array(QA, tmp_QA);

	ret_stop = EMFP->Update_Eta_given_Z(tmp_empcs, EMC);	/* Update Eta given Z_normalized. */
	if(ret_stop){
		return(ret_stop);
	}

	EMFP->Update_Mu_given_QA(tmp_empcs, tmp_QA, QA_H);	/* Update Mu given Eta and QA. */ 
	EMFP->Update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
	EMFP->Maximize_logpL(tmp_empcs, QA, QA_H, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
	tmp_QA->Update_log_Pt(tmp_QA);
	EMFP->Update_Z_modified(tmp_empcs, tmp_QA);		/* Update empcs->Z_modified (log and unnormalized). */
	tmp_R = EMFP->Compute_R(tmp_empcs, tmp_QA, QA_H);	/* Compute R(Mu, QA, Tt). */
	do{
		EMFP->Copy_empcs(tmp_empcs, empcs);
		tmp_QA->Copy_Q_matrix_array(tmp_QA, QA);
		R = tmp_R;

		EMFP->Update_Mu_given_QA(tmp_empcs, tmp_QA, QA_H);		/* Update Mu given Eta and QA. */ 
		EMFP->Update_Z_modified(tmp_empcs, tmp_QA);			/* Update empcs->Z_modified (log and unnormalized). */
		EMFP->Maximize_logpL(tmp_empcs, tmp_QA, QA_H, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
		tmp_QA->Update_log_Pt(tmp_QA);
		EMFP->Update_Z_modified(tmp_empcs, tmp_QA);			/* Update empcs->Z_modified (log and unnormalized). */
		tmp_R = EMFP->Compute_R(tmp_empcs, tmp_QA, QA_H);		/* Compute R(Mu, QA, Tt). */

		if(tmp_R < R){
			flag = 0;
			break;
		}

		cm_reltol = fabs(tmp_R / R - 1.0);
		cm_iter++;
	} while((cm_reltol > EMC->cm_reltol) && (cm_iter < EMC->cm_maxit));

	if(flag){
		EMFP->Copy_empcs(tmp_empcs, empcs);
		tmp_QA->Copy_Q_matrix_array(tmp_QA, QA);
	}
	EMC->converge_cm_iter += cm_iter;
	return(ret_stop);
} /* End of M_step_CM(). */

/* The alternative conditional M-step:
 * Update Eta given Z_normalized. Update Mu given Eta and QA. Update QA given Eta and MU.
 * These 3 updates form a cycle of a M-step. The complete logL may not be optimized in a M-step,
 * but it will increase observed logL in each iterations.
 * All of these ACM steps require to set EMC->update_flag = 1 in Em_step().
 * The maximize_logpL() in "phyclust_logpL.c",
 * update_flag = 0 for update Mu given QA. (for profile logL)
 *             = 1 for update QA given Mu. (for profile logL and ECM/AECM) */
int M_step_ACM(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA){
	int ret_stop = 0;

	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC);	/* Update Eta given Z_normalized. */
	if(ret_stop){
		return(ret_stop);
	}
	E_step_simple(empcs, QA, EMFP);
	EMFP->Update_Mu_given_QA(empcs, QA, QA_H);		/* Update Mu given Eta and QA. */ 
	E_step_simple(empcs, QA, EMFP);
	EMFP->Maximize_logpL(empcs, QA, QA_H, EMC, EMFP);	/* Update QA given Eta and Mu. */ 
	QA->Update_log_Pt(QA);
	EMC->converge_cm_iter++;
	return(ret_stop);
} /* End of M_step_ACM(). */


int Update_Eta_given_Z_ADJUST(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k, check_Eta[empcs->K], flag_out_range = 0;
	double total_sum = 0.0, tmp_sum_in_range = 0.0, tmp_sum_out_range = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			if(empcs->replication_X[n_X] == 1){
				empcs->Eta[k] += empcs->Z_normalized[n_X][k];
			} else{
				empcs->Eta[k] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
			}
		}
		total_sum += empcs->Eta[k];
	}

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;

		/* Double check if out range occurs and set to boundary if any. */
		if(empcs->Eta[k] < EMC->Eta_lower_bound){
			empcs->Eta[k] = EMC->Eta_lower_bound;
			tmp_sum_out_range += empcs->Eta[k];
			check_Eta[k] = 1;
			flag_out_range |= 1;
		} else if(empcs->Eta[k] > EMC->Eta_upper_bound){
			empcs->Eta[k] = EMC->Eta_upper_bound;
			tmp_sum_out_range += empcs->Eta[k];
			check_Eta[k] = 1;
			flag_out_range |= 1;
		} else{
			tmp_sum_in_range += empcs->Eta[k];
			check_Eta[k] = 0;
		}
	}

	/* Normalize in range part to constrained total if out range occurs. */
	if(flag_out_range == 1){
		for(k = 0; k < empcs->K; k++){
			if(check_Eta[k] == 0){
				empcs->Eta[k] *= (1.0 - tmp_sum_out_range) / tmp_sum_in_range;
			}
		}
	}

	/* Update log_Eta before return. */
	for(k = 0; k < empcs->K; k++){
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}

	#if (EMDEBUG & 4) == 4
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_ADJUST(). */

/* Return 1 for stopping em algorithm since Eta is out of range. */
int Update_Eta_given_Z_IGNORE(em_phyclust_struct *empcs, em_control *EMC){
	int n_X, k;
	double total_sum = 0.0;

	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] = 0.0;
		for(n_X = 0; n_X < empcs->N_X; n_X++){
			if(empcs->replication_X[n_X] == 1){
				empcs->Eta[k] += empcs->Z_normalized[n_X][k];
			} else{
				empcs->Eta[k] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
			}
		}
		total_sum += empcs->Eta[k];
	}
	for(k = 0; k < empcs->K; k++){
		empcs->Eta[k] /= total_sum;
		empcs->log_Eta[k] = log(empcs->Eta[k]);
	}
	for(k = 0; k < empcs->K; k++){
		if(empcs->Eta[k] < EMC->Eta_lower_bound || empcs->Eta[k] > EMC->Eta_upper_bound){
			#if (EMDEBUG & 4) == 4
				printf("  Eta[%d]=%f is out of (%f, %f)\n", k, empcs->Eta[k], EMC->Eta_lower_bound,
						EMC->Eta_upper_bound);
			#endif
			return(1);
		}
	}

	#if (EMDEBUG & 4) == 4
		printf("    Eta:");
		for(k = 0; k < empcs->K; k++){
			printf(" %f", empcs->Eta[k]);
		}
		printf("\n");
	#endif

	return(0);
} /* End of Update_Eta_given_Z_IGNORE(). */


/* For each sequence, compute log_Pt for all (n, k) cells and save in empcs->Z_modified (log and unnormalized). */
void Update_Z_modified(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_modified[n_X][k] = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					empcs->Z_modified[n_X][k] += QA->Q[k]->log_Pt[s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
					/* For gap.
					 * This can be skipped for Z_normalized, but can NOT be skipped
					 * for the complete logL, the R function and the M steps. */
				}
			}
		}
	}
} /* End of Update_Z_modified(). */


void print_status(em_phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, int verbosity){
	int k;
	if(!verbosity) return;
	if(verbosity==1){
		printf(".");
		/* fflush(stdout); */
		return;
	}
	if(verbosity==2){
		printf("%5d %12.3f\n", EMC->converge_iter, pcs->logL_observed);
		return;
	}
	if(verbosity==3) {
		printf("%5d eta", EMC->converge_iter);
		for(k = 0; k<pcs->K; k++){
			printf(" %6.4f", pcs->Eta[k]);
		}
		print_QA(QA);
		printf(" %12.3f\n", pcs->logL_observed);
		return;
	}
} /* End of print_status(). */


/* converge_flag = 0, successed to converge.
 *               = 1, out of iterations.
 *               = 2, maybe converged where logL_observed decreasing or Eta's go to 0.
 *               = 9, fail to converge.
 * update_flag = 0 for update Mu given QA. (for profile logL)
 *             = 1 for update QA given Mu. (for profile logL and ECM/AECM) */
int Check_convergence_em(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, Q_matrix_array *new_QA,
		Q_matrix_array *org_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	if(new_empcs->logL_observed < org_empcs->logL_observed){	/* logL decreasing. */
		if(EMC->update_flag == 0){				/* change to update Q, Tt given Mu. */
			EMC->update_flag = 1;
		} else{							/* otherwise stop. */
			EMC->converge_flag = 9;
			EMC->converge_error = new_empcs->logL_observed - org_empcs->logL_observed;
			ret_stop = 1;
		}

		/* Restore to orginal setting with heigher logL. */
		EMFP->Copy_empcs(org_empcs, new_empcs);			/* Repalce object with higher logL. */
		org_QA->Copy_Q_matrix_array(org_QA, new_QA);
		QA_H->Copy_Q_matrix_array(QA_H, org_QA);
	} else{
		if(EMC->update_flag == 1){			/* logL increasing and current is update QA, Tt given Mu*/
			EMC->update_flag = 0;			/* change to original plan. */
		}
	}
	return(ret_stop);
} /* End of Check_convergence_em(). */

int Check_convergence_org(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, Q_matrix_array *new_QA,
		Q_matrix_array *org_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	if(new_empcs->logL_observed < org_empcs->logL_observed){	/* logL decreasing. */
		EMC->converge_flag = 9;
		EMC->converge_error = new_empcs->logL_observed - org_empcs->logL_observed;
		ret_stop = 1;

		/* Restore to orginal setting with heigher logL. */
		EMFP->Copy_empcs(org_empcs, new_empcs);			/* Repalce object with higher logL. */
		org_QA->Copy_Q_matrix_array(org_QA, new_QA);
		QA_H->Copy_Q_matrix_array(QA_H, org_QA);
	}
	return(ret_stop);
} /* End of Check_convergence_org(). */


void Em_step(em_phyclust_struct *org_empcs, Q_matrix_array *org_QA, em_control *EMC, em_fp *EMFP){
	/* Initials are taken and Results are stored in org_empcs, org_QA and EMC. */
	int ret_stop = 0;
	em_phyclust_struct *new_empcs, *tmp_empcs;		/* new_empcs is for the new results, (i+1)-th step.
								 * tmp_empcs is for CM temporary storages. */
	Q_matrix_array *new_QA, *tmp_QA, *QA_H;			/* new_QA is for the new results, (i+1)-th step..
								 * tmp_QA is for the CM temporary storages.
								 * QA_H is for the entropy of E-steps, i-th step. */

	reset_em_control(EMC);
	new_empcs = duplicate_em_phyclust_struct(org_empcs);
	new_QA = duplicate_Q_matrix_array(org_QA);
	QA_H = duplicate_Q_matrix_array(org_QA);
	if(EMC->em_method == ECM){
		tmp_empcs = duplicate_em_phyclust_struct(org_empcs);
		tmp_QA = duplicate_Q_matrix_array(org_QA);
	} else{
		tmp_empcs = NULL;
		tmp_QA = NULL;
	}

	#if (EMDEBUG & 1) == 1
		int k;
		printf("Start: EM\n");
		printf("  Eta:");
		for(k = 0; k < new_empcs->K; k++){
			printf(" %.4f", new_empcs->Eta[k]);
		}
		printf("\n");
		print_QA(new_QA);
		printf("\n");
		if(EMC->se_type == SE_YES){
			new_empcs->SE_P->Print_f_err(new_empcs->SE_P);
		}
	#endif

	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
	EMFP->E_step_logL_observed(new_empcs, new_QA, EMFP);
	#if (EMDEBUG & 1) == 1
		double tmp_logL;
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter - 1, EMC->update_flag);
		tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
		if(is_finite(new_empcs->logL_observed)){
			printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
		} else{
			printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
		}
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
	#endif
	do{
		#if verbosity_em_step > 0
			print_status(new_empcs, new_QA, EMC, verbosity_em_step);
		#endif

		EMFP->Copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(org_QA, QA_H);
		new_QA->Copy_Q_matrix_array(new_QA, org_QA);

		#if (EMDEBUG & 2) == 2
			double tmp_R, tmp_Q, tmp_obs;
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA, org_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA, org_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("  M-step:\n");
			printf("    init.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif
		ret_stop = EMFP->M_step(new_empcs, new_QA, org_QA, EMC, EMFP, tmp_empcs, tmp_QA);
		if(ret_stop){
			EMC->converge_flag = 2;	/* Eta < 1/N or > 1 - 1/N. */
			break;
		}
		#if (EMDEBUG & 2) == 2
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA, org_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA, org_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("    conv.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif

		EMFP->E_step_logL_observed(new_empcs, new_QA, EMFP);

		EMC->converge_eps = fabs(new_empcs->logL_observed / org_empcs->logL_observed - 1.0);
		EMC->converge_iter++;
		#if (EMDEBUG & 1) == 1
			tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
			if(is_finite(new_empcs->logL_observed)){
				printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			} else{
				printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			}
			printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
		#endif

		ret_stop = EMFP->Check_convergence(new_empcs, org_empcs, new_QA, org_QA, QA_H, EMC, EMFP);
		if(ret_stop){
			break;	/* logL is decreasing and fails after switch. */
		}
	} while((EMC->converge_eps > EMC->EM_eps) && (EMC->converge_iter < EMC->EM_iter));

	#if verbosity_em_step > 1
		printf("\n");
	#endif

	if(EMC->converge_iter > EMC->EM_iter){
		EMC->converge_flag = 1;
	}

	if(EMC->converge_flag < 2){
		EMFP->Copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);
	}

	free_em_phyclust_struct(new_empcs);
	free_Q_matrix_array(new_QA);
	free_Q_matrix_array(QA_H);
	if(EMC->em_method == ECM){
		free_em_phyclust_struct(tmp_empcs);
		free_Q_matrix_array(tmp_QA);
	}
} /* End of Em_step(). */

void Short_em_step(em_phyclust_struct *org_empcs, Q_matrix_array *org_QA, em_control *EMC, em_fp *EMFP){
	/* Initials are taken and Results are stored in org_empcs, org_QA and EMC. */
	double logL_0;
	int ret_stop = 0;
	em_phyclust_struct *new_empcs, *tmp_empcs;		/* new_empcs is for the new results, (i+1)-th step.
								 * tmp_empcs is for CM temporary storages. */
	Q_matrix_array *new_QA, *tmp_QA, *QA_H;			/* new_QA is for the new results, (i+1)-th step..
								 * tmp_QA is for the CM temporary storages.
								 * QA_H is for the entropy of E-steps, i-th step. */

	reset_em_control(EMC);
	new_empcs = duplicate_em_phyclust_struct(org_empcs);
	new_QA = duplicate_Q_matrix_array(org_QA);
	QA_H = duplicate_Q_matrix_array(org_QA);
	if(EMC->em_method == ECM){
		tmp_empcs = duplicate_em_phyclust_struct(org_empcs);
		tmp_QA = duplicate_Q_matrix_array(org_QA);
	} else{
		tmp_empcs = NULL;
		tmp_QA = NULL;
	}

	#if (EMDEBUG & 1) == 1
		int k;
		double tmp_logL;
		printf("Start: EM\n");
		printf("  Eta:");
		for(k = 0; k < new_empcs->K; k++){
			printf(" %.4f", new_empcs->Eta[k]);
		}
		printf("\n");
		print_QA(new_QA);
	#endif

	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
	EMFP->E_step_logL_observed(new_empcs, new_QA, EMFP);
	logL_0 = new_empcs->logL_observed;
	#if (EMDEBUG & 1) == 1
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter - 1, EMC->update_flag);
		EMFP->Update_Z_modified(org_empcs, org_QA);
		tmp_logL = EMFP->LogL_observed(org_empcs, org_QA);
		if(is_finite(org_empcs->logL_observed)){
			printf("  logL: %.8f, %.8f [indep fcn]\n", org_empcs->logL_observed, tmp_logL);
		} else{
			printf("  logL: %.8e, %.8e [indep fcn]\n", org_empcs->logL_observed, tmp_logL);
		}
		printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
	#endif
	do{
		EMFP->Copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(org_QA, QA_H);
		new_QA->Copy_Q_matrix_array(new_QA, org_QA);

		#if (EMDEBUG & 2) == 2
			double tmp_R, tmp_Q, tmp_obs;
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA, org_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA, org_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("  M-step:\n");
			printf("    init.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif
		ret_stop = EMFP->M_step(new_empcs, new_QA, org_QA, EMC, EMFP, tmp_empcs, tmp_QA);
		if(ret_stop){
			EMC->converge_flag = 2;	/* Eta < 1/N or > 1 - 1/N. */
			break;
		}
		#if (EMDEBUG & 2) == 2
			tmp_R = EMFP->LogL_profile(new_empcs, new_QA, org_QA);
			tmp_Q = EMFP->LogL_complete(new_empcs, new_QA, org_QA);
			tmp_obs = EMFP->LogL_observed(new_empcs, new_QA);
			printf("    conv.: R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		#endif

		EMFP->E_step_logL_observed(new_empcs, new_QA, EMFP);

		EMC->converge_eps = (org_empcs->logL_observed - new_empcs->logL_observed) /
			(logL_0 - new_empcs->logL_observed);
		EMC->converge_iter++;
		#if (EMDEBUG & 1) == 1
			tmp_logL = EMFP->LogL_observed(new_empcs, new_QA);
			if(is_finite(new_empcs->logL_observed)){
				printf("  logL: %.8f, %.8f [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			} else{
				printf("  logL: %.8e, %.8e [indep fcn]\n", new_empcs->logL_observed, tmp_logL);
			}
			printf("iter: %d, update_flag: %d\n", EMC->converge_iter, EMC->update_flag);
		#endif

		ret_stop = EMFP->Check_convergence(new_empcs, org_empcs, new_QA, org_QA, QA_H, EMC, EMFP);
		if(ret_stop){
			break;	/* logL is decreasing and fails after switch. */
		}
	} while((EMC->converge_eps > EMC->short_eps) && (EMC->converge_iter < EMC->short_iter));

	if(EMC->converge_iter > EMC->short_iter){
		EMC->converge_flag = 1;
	}

	if(EMC->converge_flag < 2){
		EMFP->Copy_empcs(new_empcs, org_empcs);
		org_QA->Copy_Q_matrix_array(new_QA, org_QA);
	}

	free_em_phyclust_struct(new_empcs);
	free_Q_matrix_array(new_QA);
	free_Q_matrix_array(QA_H);
	if(EMC->em_method == ECM){
		free_em_phyclust_struct(tmp_empcs);
		free_Q_matrix_array(tmp_QA);
	}
} /* End of Short_em_step(). */


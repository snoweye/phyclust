/* This file contains all fucntions for initialization. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_edist.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"
#include "phyclust_tool.h"
#include "phyclust_ape_nj.h"


/* Define double rdunif(). */
#ifdef RMATH_H
	#ifdef MATHLIB_STANDALONE
		/* Require set_seed(SEED1, SEED2) & get_seed(SEED1, SEED2) to call this. */
		#include <Rmath.h>
		int rdunif(int n){
			return((int) floor(n * unif_rand()));
		}
	#else
		#include <R.h>
		#include <Rmath.h>
		int rdunif(int n){
			int ret = 0;
			GetRNGstate();
			ret = (int) floor(n * runif(0, 1));
			PutRNGstate();
			return(ret);
		}
	#endif
#else
	/* Require srand(seed) to call this. */
	#include <math.h>
	int rdunif(int n){
		return((int) floor(n * (double) rand() / ((double) RAND_MAX + 1.0)));
	}
#endif


/* Simple random sample withor replacement. Choose k from 0 to n-1.
 * Modified from Dr. Maitra's code. */
void srswor(int n, int k, int *x){
	int i, j;
	int *tmp_x = allocate_int_1D(n);

	for(i = 0; i < n; i++){
		tmp_x[i] = i;
	}

	for(i = 0; i < k; i++){
		j = rdunif(n);
		x[i] = tmp_x[j];
		tmp_x[j] = tmp_x[--n];
	}

	free(tmp_x);
}


/* Copy from phyclust_em.c and modified for initialization method. */
int init_m_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	Q_matrix_array *QA_H;

	initialize_count_Mu_X_and_gap(empcs);		/* Initialize count_Mu_X and count_Mu_X_gap. */
	ret_stop = EMFP->Update_Eta_given_Z(empcs, EMC);	/* Find Eta. */
	if(ret_stop > 0){
		return(ret_stop);
	}

	/* Find QA. */
	EMC->update_flag = 1;		/* For update QA, given Mu. */
	QA_H = duplicate_Q_matrix_array(QA);
	ret_stop = EMFP->Maximize_logpL(empcs, QA, QA_H, EMC, EMFP);
	QA->Update_log_Pt(QA);
	EMC->update_flag = 0;		/* Reset to 0 for update Mu, given QA. */

	free_Q_matrix_array(QA_H);
	return(ret_stop);
} /* End of init_m_step(). */

int check_all_min_n_class(int K, int *n_class, int min_n_class){
	int k, ret = 1;

	for(k = 0; k < K; k++){
		ret &= (n_class[k] >= min_n_class);
	}

	return(ret);
} /* End of check_n_class(). */

void assign_Mu_by_class(int N_X_org, int K, int L, int ncode, int gap_index, int *class_id, int **X_org, int **Mu){
	int i, n_X_org, k, l;
	int count_K_N[K][ncode], count_N[ncode], tmp_count;

	for(l = 0; l < L; l++){
		/* Find the most common nucleotide in all sequences for replacing when tie in each k. */
		for(i = 0; i < ncode; i++){
			count_N[i] = 0;
			for(k = 0; k < K; k++){
				count_K_N[k][i] = 0;
			}
		}
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			if(X_org[n_X_org][l] == gap_index || X_org[n_X_org][l] == MISSING_ALLELE){
				continue;
			}
			count_N[X_org[n_X_org][l]]++;
			count_K_N[class_id[n_X_org]][X_org[n_X_org][l]]++;
		}

		/* Find the most common nucleotide for each k. */
		for(k = 0; k < K; k++){
			tmp_count = -1;
			for(i = 0; i < ncode; i++){
				if(count_K_N[k][i] > tmp_count){
					tmp_count = count_K_N[k][i];
					Mu[k][l] = i;
				} else if(count_K_N[k][i] == tmp_count && count_N[i] > count_N[Mu[k][l]]){
					Mu[k][l] = i;	/* Return the most common nucleotide if tie. */
				}
			}
		}
	}
} /* End of assign_Mu_by_class(). */

void find_consensus_Mu(int N_X_org, int L, int ncode, int gap_index, int **X_org, int *consensus_Mu){
	int i, n_X_org, l, flag, max_i;
	int count_L_N[L][ncode], count_N[ncode], tmp_count;

	/* Summarize overall. */
	for(i = 0; i < ncode; i++){
		count_N[i] = 0;
		for(l = 0; l < L; l++){
			count_L_N[l][i] = 0;
		}
	}
	for(l = 0; l < L; l++){
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			if(X_org[n_X_org][l] == gap_index || X_org[n_X_org][l] == MISSING_ALLELE){
				continue;
			}
			count_L_N[l][X_org[n_X_org][l]]++;
			count_N[X_org[n_X_org][l]]++;
		}
	}

	max_i = 0;
	tmp_count = count_N[0];
	for(i = 1; i < ncode; i++){
		if(count_N[i] > tmp_count){
			max_i = i;
			tmp_count = count_N[i];
		}
	}

	/* Find the consensus. */
	for(l = 0; l < L; l++){
		flag = 0;
		for(i = 0; i < ncode; i++){
			if(count_L_N[l][i] > 0){
				flag = 1;
				break;
			}
		}

		if(flag){
			tmp_count = -1;
			for(i = 0; i < ncode; i++){
				if(count_L_N[l][i] > tmp_count){
					tmp_count = count_L_N[l][i];
					consensus_Mu[l] = i;
				} else if(count_L_N[l][i] == tmp_count && count_N[i] > count_N[consensus_Mu[l]]){
					/* tie, compare the overall. */
					consensus_Mu[l] = i;
				}
			}
		} else{	/* All are gap, replace by the max. */
			consensus_Mu[l] = max_i;
		}
	}
} /* End of find_consensus_Mu(). */

void find_consensus_Mu_gap(int N_X_org, int L, int ncode_wigap, int gap_index, int **X_org, int *consensus_Mu){
	int i, n_X_org, l, flag, max_i;
	int count_L_N[L][ncode_wigap], count_N[ncode_wigap], tmp_count;

	/* Summarize overall. */
	for(i = 0; i < ncode_wigap; i++){
		count_N[i] = 0;
		for(l = 0; l < L; l++){
			count_L_N[l][i] = 0;
		}
	}
	for(l = 0; l < L; l++){
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			if(X_org[n_X_org][l] == MISSING_ALLELE){
				continue;
			}
			count_L_N[l][X_org[n_X_org][l]]++;
			count_N[X_org[n_X_org][l]]++;
		}
	}

	max_i = 0;
	tmp_count = count_N[0];
	for(i = 1; i < ncode_wigap; i++){
		if(count_N[i] > tmp_count){
			max_i = i;
			tmp_count = count_N[i];
		}
	}

	/* Find the consensus. */
	for(l = 0; l < L; l++){
		flag = 0;
		for(i = 0; i < ncode_wigap; i++){
			if(count_L_N[l][i] > 0){
				flag = 1;
				break;
			}
		}

		if(flag){
			tmp_count = -1;
			for(i = 0; i < ncode_wigap; i++){
				if(count_L_N[l][i] > tmp_count){
					tmp_count = count_L_N[l][i];
					consensus_Mu[l] = i;
				} else if(count_L_N[l][i] == tmp_count && count_N[i] > count_N[consensus_Mu[l]]){
					/* tie, compare the overall. */
					consensus_Mu[l] = i;
				}
			}
		} else{	/* All are gap, replace by the max. */
			consensus_Mu[l] = max_i;
		}
	}
} /* End of find_consensus_Mu_gap(). */




int Update_init_manually(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int n_X, k, ret_stop = 0;

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][empcs->class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
	}

	reset_Q_matrix_array(QA);
	if(EMC->se_type == SE_YES){
		reset_SE_P_matrix(empcs->SE_P);
	}
	assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->gap_index,
				empcs->class_id, empcs->X_org, empcs->Mu);
	ret_stop = init_m_step(empcs, QA, EMC, EMFP);
	if(ret_stop > 0){
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: Initialization error.\n");
		#endif
		return(ret_stop);
	}
	if(!is_finite(EMFP->LogL_observed(empcs, QA))){
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: manual initialization leads to non-finite observed log likelihood\n");
		#endif
		return(1);
	}

	return(ret_stop);
} /* End of Update_init_manually(). */




/* These functions with "_unique" pick centers from unique sequences
 * unlike other methods, then maps id back to usual "X_org". */

/* Randomly pick Mu from X. */
int Update_init_random_Mu_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X, k, l, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], tmp_id;
	int consensus_Mu[L];
	double tmp, tmp_min, init_logL_observed = 0.0;
	edist_struct *eds;

	find_consensus_Mu(empcs->N_X_org, L, empcs->ncode, empcs->gap_index, empcs->X_org, consensus_Mu);
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		if(EMC->se_type == SE_YES){
			reset_SE_P_matrix(empcs->SE_P);
		}

		srswor(N_X, K, center_id);
		
		/* Assign Mu by centers. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
		for(n_X = 0; n_X < N_X; n_X++){
			tmp_min = eds->get_pair_edist(eds, n_X, center_id[0]);
			tmp_id = 0;
			for(k = 1; k < K; k++){
				tmp = eds->get_pair_edist(eds, n_X, center_id[k]);
				if(tmp < tmp_min){
					tmp_min = tmp;
					tmp_id = k;
				}
			}
	
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][tmp_id] = 1.0;
			empcs->n_class[tmp_id] += empcs->replication_X[n_X];
		}

		/* Replace gaps by the concensus. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				if(empcs->Mu[k][l] == empcs->gap_index || empcs->Mu[k][l] == MISSING_ALLELE){
					empcs->Mu[k][l] = consensus_Mu[l];
				}
			}
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop > 0){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(ret_stop);
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(1);
		}
	}

	free_edist_struct(eds);
	return(ret_stop);
} /* End of Update_init_random_Mu_unique(). */

int Update_init_random_Mu_unique_label(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X, k, l, N_X = empcs->N_X, K = empcs->K, L = empcs->L,
		N_X_unlabeled = empcs->N_X_unlabeled, K_labeled = empcs->K_labeled, K_unlabeled = K - K_labeled,
		tmp_n, tmp_k;
	int center_id[K], tmp_center_id[K], tmp_id;
	int consensus_Mu[L];
	double tmp, tmp_min, init_logL_observed = 0.0;
	edist_struct *eds;

	find_consensus_Mu(empcs->N_X_org, L, empcs->ncode, empcs->gap_index, empcs->X_org, consensus_Mu);
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		// reset_SE_P_matrix(empcs->SE_P);

		/* Randomly pick centers from X. */
		for(k = 0; k < K_labeled; k++){
			tmp_n = 0;
			for(n_X = 0; n_X < empcs->N_X_labeled; n_X++){
				if(empcs->label_semi[n_X] == k){
					tmp_n++;
				}
			}
			srswor(tmp_n, 1, tmp_center_id);

			tmp_n = -1;
			for(n_X = 0; n_X < empcs->N_X_labeled; n_X++){
				if(empcs->label_semi[n_X] == k){
					tmp_n++;
					if(tmp_n == tmp_center_id[0]){
						break;
					}
				}
			}
			tmp_n = n_X;

			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->X[n_X] == empcs->X_labeled[tmp_n]){
					center_id[k] = n_X;
					break;
				}
			}
		}

		if(K_unlabeled > 0){
			srswor(N_X_unlabeled, K_unlabeled, tmp_center_id);

			for(k = 0; k < K_unlabeled; k++){
				tmp_n = tmp_center_id[k];
				tmp_k = k + K_labeled;
				for(n_X = 0; n_X < N_X; n_X++){
					if(empcs->X[n_X] == empcs->X_unlabeled[tmp_n]){
						center_id[tmp_k] = n_X;
						break;
					}
				}
			}
		}

		/* Assign Mu by centers. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
		for(n_X = 0; n_X < N_X; n_X++){
			tmp_min = eds->get_pair_edist(eds, n_X, center_id[0]);
			tmp_id = 0;
			for(k = 1; k < K; k++){
				tmp = eds->get_pair_edist(eds, n_X, center_id[k]);
				if(tmp < tmp_min){
					tmp_min = tmp;
					tmp_id = k;
				}
			}
	
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][tmp_id] = 1.0;
			empcs->n_class[tmp_id] += empcs->replication_X[n_X];
		}

		/* Replace gaps by the concensus. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				if(empcs->Mu[k][l] == empcs->gap_index || empcs->Mu[k][l] == MISSING_ALLELE){
					empcs->Mu[k][l] = consensus_Mu[l];
				}
			}
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop > 0){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(ret_stop);
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(1);
		}
	}

	free_edist_struct(eds);
	return(ret_stop);
} /* End of Update_init_random_Mu_unique(). */




/* Pick clusters by cutting the longest K internal branches of neighbor-joining tree.
 * There is no randomness for this method, so the following settings may be suggested.
 * EMC->init_procedure = exhaustEM;
 * EMC->exhaust_iter = 1; */
int Update_init_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int n_X_org, n_X, k, ret_stop = 0;
	int N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int largest_branch_id[N_X - 3], class_id[N_X];
	double init_logL_observed;
	edist_struct *eds;
	nj_struct *njs;
	
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);
	njs = initialize_nj_struct(N_X);
	njs->D = eds->EDM[0];
	phyclust_ape_nj(njs);
	if(! check_njs(njs)){
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: NJ may be not valid!\n");
		#endif
		print_njs(njs->n_edge, njs);
		free_edist_struct(eds);
		free_nj_struct(njs);
		return(1);
	}
	#if INITDEBUG > 0
		print_njs(njs->n_edge, njs);
	#endif

	search_largest_branch(njs, largest_branch_id);
	ret_stop = assign_class_by_njs_branch(K, njs, largest_branch_id, class_id);
	if(ret_stop != 0){
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: Class assignment fails.\n");
		#endif
		free_edist_struct(eds);
		free_nj_struct(njs);
		return(1);
	}

	/* Assign Mu and recreate Z_normalized. */
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		empcs->class_id[n_X_org] = class_id[empcs->map_X_org_to_X[n_X_org]];
	}
	assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->gap_index,
				empcs->class_id, empcs->X_org, empcs->Mu);
	for(k = 0; k < K; k++){
		empcs->n_class[k] = 0;
	}
	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][class_id[empcs->map_X_org_to_X[empcs->map_X_to_X_org[n_X]]]] = 1.0;
	}
	for(n_X = 0; n_X < N_X; n_X++){
		empcs->n_class[class_id[n_X]] += empcs->replication_X[n_X];
	}

	if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			free_nj_struct(njs);
			return(ret_stop);
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			free_nj_struct(njs);
			return(1);
		}
	} else{
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class,
					INIT_METHOD[EMC->init_method]);
		#endif
		free_edist_struct(eds);
		free_nj_struct(njs);
		return(1);
	}

	free_edist_struct(eds);
	free_nj_struct(njs);
	return(0);
} /* End of Update_init_nj_unique(). */


/* Pick clusters by randomly cutting the longest K internal branches of neighbor-joining tree. */
int Update_init_random_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X_org, n_X, k;
	int N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int random_branch_id[N_X - 3], class_id[N_X];
	double init_logL_observed;
	edist_struct *eds;
	nj_struct *njs;
	
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);
	njs = initialize_nj_struct(N_X);
	njs->D = eds->EDM[0];
	phyclust_ape_nj(njs);
	if(! check_njs(njs)){
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: NJ may be not valid!\n");
		#endif
		print_njs(njs->n_edge, njs);
		free_edist_struct(eds);
		free_nj_struct(njs);
		return(1);
	}
	#if INITDEBUG > 0
		print_njs(njs->n_edge, njs);
	#endif

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		if(EMC->se_type == SE_YES){
			reset_SE_P_matrix(empcs->SE_P);
		}

		/* Randomly pick mu from X. */
		random_branch(njs, random_branch_id);
		ret_stop = assign_class_by_njs_branch(K, njs, random_branch_id, class_id);
		if(ret_stop != 0){
			free_edist_struct(eds);
			free_nj_struct(njs);
			continue;
		}

		/* Assign Mu and recreate Z_normalized. */
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			empcs->class_id[n_X_org] = class_id[empcs->map_X_org_to_X[n_X_org]];
		}
		assign_Mu_by_class(empcs->N_X_org, empcs->K, empcs->L, empcs->ncode, empcs->gap_index,
					empcs->class_id, empcs->X_org, empcs->Mu);
		for(k = 0; k < K; k++){
			empcs->n_class[k] = 0;
		}
		for(n_X = 0; n_X < N_X; n_X++){
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][class_id[empcs->map_X_org_to_X[empcs->map_X_to_X_org[n_X]]]]
				= 1.0;
		}
		for(n_X = 0; n_X < N_X; n_X++){
			empcs->n_class[class_id[n_X]] += empcs->replication_X[n_X];
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop > 0){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			free_nj_struct(njs);
			return(ret_stop);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			free_nj_struct(njs);
			return(1);
		}
	}

	free_edist_struct(eds);
	free_nj_struct(njs);
	return(ret_stop);
} /* End of Update_init_random_nj_unique(). */




int Update_init_k_medoids(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0, ret_stop = 0;
	int n_X_org, n_X, k, l, N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X_org];
	int consensus_Mu[L];
	double init_logL_observed;
	edist_struct *eds;

	find_consensus_Mu(empcs->N_X_org, L, empcs->ncode, empcs->gap_index, empcs->X_org, consensus_Mu);
	eds = initialize_edist_struct_UT(EMC->edist_model, N_X_org, L, empcs->X_org);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		if(EMC->se_type == SE_YES){
			reset_SE_P_matrix(empcs->SE_P);
		}

		assign_class_unique_by_k_medoids(N_X_org, K, eds->EDM, N_X, empcs->map_X_to_X_org,
							center_id, class_id);

		/* Pick mu from X. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X_org[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
		for(n_X = 0; n_X < N_X; n_X++){
			for(k = 0; k < K; k++){
				empcs->Z_normalized[n_X][k] = 0.0;
			}
			empcs->Z_normalized[n_X][class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
		}
		for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
			empcs->n_class[class_id[n_X_org]]++;
		}

		/* Replace gaps by the concensus. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				if(empcs->Mu[k][l] == empcs->gap_index || empcs->Mu[k][l] == MISSING_ALLELE){
					empcs->Mu[k][l] = consensus_Mu[l];
				}
			}
		}

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			ret_stop = init_m_step(empcs, QA, EMC, EMFP);
			if(ret_stop > 0){
				continue;
			}
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter >= EMC->max_init_iter){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(ret_stop);	
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(1);
		}
	}

	free_edist_struct(eds);
	return(ret_stop);
} /* End of Update_init_k_medoids(). */




/* Pick clusters by PAM.
 * There is no randomness for this method, so the following settings may be suggested.
 * EMC->init_procedure = exhaustEM;
 * EMC->exhaust_iter = 1; */
int Update_init_pam(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int ret_stop = 0;
	int n_X_org, n_X, k, l, N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X_org];
	int consensus_Mu[L];
	double init_logL_observed;
	edist_struct *eds;
	
	find_consensus_Mu(empcs->N_X_org, L, empcs->ncode, empcs->gap_index, empcs->X_org, consensus_Mu);
	eds = initialize_edist_struct_LT_pam(EMC->edist_model, N_X_org, L, empcs->X_org);
	assign_class_by_pam(N_X_org, K, eds->EDM, center_id, class_id);

	/* Pick mu from X. */
	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			empcs->Mu[k][l] = empcs->X_org[center_id[k]][l];
		}
		empcs->n_class[k] = 0;
	}

	/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			empcs->Z_normalized[n_X][k] = 0.0;
		}
		empcs->Z_normalized[n_X][class_id[empcs->map_X_to_X_org[n_X]]] = 1.0;
	}
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		empcs->n_class[class_id[n_X_org]]++;
	}

	/* Replace gaps by the concensus. */
	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			if(empcs->Mu[k][l] == empcs->gap_index || empcs->Mu[k][l] == MISSING_ALLELE){
				empcs->Mu[k][l] = consensus_Mu[l];
			}
		}
	}

	if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
		ret_stop = init_m_step(empcs, QA, EMC, EMFP);
		if(ret_stop > 0){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initialization error. (%s)\n", INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(ret_stop);
		}
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			#if PRINT_ERROR > 0
				fprintf_stderr("PE: Initial logL_observed is not finit. (%s)\n",
						INIT_METHOD[EMC->init_method]);
			#endif
			free_edist_struct(eds);
			return(1);
		}
	} else{
		#if PRINT_ERROR > 0
			fprintf_stderr("PE: Initialization is not valid for min_n_class = %d. (%s)\n",
					EMC->min_n_class, INIT_METHOD[EMC->init_method]);
		#endif
		free_edist_struct(eds);
		return(1);
	}

	free_edist_struct(eds);
	return(ret_stop);
} /* End of Update_init_pam(). */




/* Sample alleles for each locus by the frequence of alleles. */
int Update_init_sampleL_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	/* WCC: TBD. */

	return(0);
} /* End of Update_init_sampleL_unique(). */




/* ----- For debug. ----- */
void print_consensus_Mu(em_phyclust_struct *empcs,  int *consensus_Mu){
	int l;

	for(l = 0; l < empcs->L; l++){
	#if PRINT_CODE_TYPE == 0
		if(empcs->code_type == NUCLEOTIDE){
			printf("%c ", NUCLEOTIDE_CODE[consensus_Mu[l]]);
		} else if(empcs->code_type == SNP){
			printf("%c ", SNP_CODE[consensus_Mu[l]]);
		}
	#else
		if(empcs->code_type == NUCLEOTIDE){
			printf("%c ", NUCLEOTIDE_ID[consensus_Mu[l]]);
		} else if(empcs->code_type == SNP){
			printf("%c ", SNP_ID[consensus_Mu[l]]);
		}
	#endif
	}
} /* End of print_consensus_Mu(). */


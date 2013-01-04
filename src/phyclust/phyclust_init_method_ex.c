/* This file contains all fucntions for initialization with execeptions. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_edist.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"
#include "phyclust_tool.h"
#include "phyclust_ape_nj.h"


void Update_init_k_medoids_by_X_org(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0;
	int n_X_org, n_X, k, l, N_X_org = empcs->N_X_org, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X_org];
	double init_logL_observed;
	edist_struct *eds;

	eds = initialize_edist_struct_UT(EMC->edist_model, N_X_org, L, empcs->X_org);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		if(EMC->se_type == SE_YES){
			reset_SE_P_matrix(empcs->SE_P);
		}

		assign_class_by_k_medoids(N_X_org, K, eds->EDM, center_id, class_id);

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

		if(check_all_min_n_class(K, empcs->n_class, EMC->min_n_class)){
			init_m_step(empcs, QA, EMC, EMFP);
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter > EMC->max_init_iter){
		printf("Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class,
				INIT_METHOD[EMC->init_method]);
		printf("Reach the maximum initial iterations. (%s)\n", INIT_METHOD[EMC->init_method]);
		init_m_step(empcs, QA, EMC, EMFP);
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			printf("Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	}

	free_edist_struct(eds);
} /* End of Update_init_k_medoids_by_X_org(). */

void Update_init_k_medoids_by_X(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int init_iter = 0;
	int n_X, k, l, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X];
	double init_logL_observed;
	edist_struct *eds;

	eds = initialize_edist_struct_UT(EMC->edist_model, N_X, L, empcs->X);

	while(init_iter < EMC->max_init_iter){
		init_iter++;
		reset_Q_matrix_array(QA);
		if(EMC->se_type == SE_YES){
			reset_SE_P_matrix(empcs->SE_P);
		}

		assign_class_by_k_medoids(N_X, K, eds->EDM, center_id, class_id);

		/* Pick mu from X. */
		for(k = 0; k < K; k++){
			for(l = 0; l < L; l++){
				empcs->Mu[k][l] = empcs->X[center_id[k]][l];
			}
			empcs->n_class[k] = 0;
		}

		/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
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
			init_m_step(empcs, QA, EMC, EMFP);
			init_logL_observed = EMFP->LogL_observed(empcs, QA);
			if(is_finite(init_logL_observed)){
				break;
			}
		}
	}

	if(init_iter > EMC->max_init_iter){
		printf("Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class,
				INIT_METHOD[EMC->init_method]);
		printf("Reach the maximum initial iterations. (%s)\n", INIT_METHOD[EMC->init_method]);
		init_m_step(empcs, QA, EMC, EMFP);
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			printf("Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	}

	free_edist_struct(eds);
} /* End of Update_init_k_medoids_by_X(). */


/* Pick clusters by PAM.
 * There is no randomness for this method, so the following settings may be suggested.
 * EMC->init_procedure = exhaustEM;
 * EMC->exhaust_iter = 1; */
void Update_init_pam_by_X(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	int n_X, k, l, N_X = empcs->N_X, K = empcs->K, L = empcs->L;
	int center_id[K], class_id[N_X];
	double init_logL_observed;
	edist_struct *eds;
	
	eds = initialize_edist_struct_LT_pam(EMC->edist_model, N_X, L, empcs->X);
	assign_class_by_pam(N_X, K, eds->EDM, center_id, class_id);

	/* Pick mu from X. */
	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			empcs->Mu[k][l] = empcs->X[center_id[k]][l];
		}
		empcs->n_class[k] = 0;
	}

	/* Assign X to the nearest mu by distance, and recreate Z_normalized. */
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
		init_m_step(empcs, QA, EMC, EMFP);
		init_logL_observed = EMFP->LogL_observed(empcs, QA);
		if(!is_finite(init_logL_observed)){
			printf("Initial logL_observed is not finit. (%s)\n", INIT_METHOD[EMC->init_method]);
			exit(1);
		}
	} else{
		printf("Initialization is not valid for min_n_class = %d. (%s)\n", EMC->min_n_class,
				INIT_METHOD[EMC->init_method]);
		exit(1);
	}

	free_edist_struct(eds);
} /* End of Update_init_pam_by_X(). */


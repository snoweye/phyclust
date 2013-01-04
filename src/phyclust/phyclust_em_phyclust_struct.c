/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"
#include "phyclust_se_em.h"


/* Initial a em_phyclust structure by a phyclust structure and a QA. */
em_phyclust_struct* initialize_em_phyclust_struct(phyclust_struct *pcs){
	int i, j, N_X_org = pcs->N_X_org, N_X = pcs->N_X, L = pcs->L, K = pcs->K;
	em_phyclust_struct *empcs;

	empcs = (em_phyclust_struct*) malloc(sizeof(em_phyclust_struct));
	empcs->code_type = pcs->code_type;
	empcs->ncode = pcs->ncode;
	empcs->gap_index = pcs->gap_index;
	empcs->gap_flag = pcs->gap_flag;
	empcs->N_X_org = N_X_org;
	empcs->N_X = N_X;
	empcs->N_seg_site = pcs->N_seg_site;
	empcs->L = L;
	empcs->K = K;
	empcs->X_org = pcs->X_org;
	empcs->X = pcs->X;
	empcs->map_X_org_to_X = pcs->map_X_org_to_X;
	empcs->map_X_to_X_org = pcs->map_X_to_X_org;
	empcs->replication_X = pcs->replication_X;
	empcs->seg_site_id = pcs->seg_site_id;
	empcs->class_id = pcs->class_id;

	/* For EM. */
	empcs->n_class = allocate_int_1D(K);
	empcs->Mu = allocate_int_RT(K, L); 
	for(i = 0; i < K; i++){
		/* This loop will be replaced by initialization functions.
		 * pcs->MU have 0 as default, and are assigend for exhausted EM. */
		for(j = 0; j < L; j++){
			empcs->Mu[i][j] = pcs->Mu[i][j];
		}
	}
	empcs->Z_modified = allocate_double_RT(N_X, K); 
	empcs->Z_normalized = allocate_double_RT(N_X, K); 
	empcs->Eta = allocate_double_1D(K); 
	empcs->log_Eta = allocate_double_1D(K); 
	for(i = 0; i < K; i++){
		empcs->Eta[i] = pcs->Eta[i];
		empcs->log_Eta[i] = log(pcs->Eta[i]);
	}
	empcs->logL_observed = 0.0;
	empcs->count_Mu_X = allocate_int_RT_4D(N_X, pcs->K, pcs->ncode, pcs->ncode);
	if(empcs->gap_flag){
		empcs->count_Mu_X_gap = allocate_int_RT_3D(N_X, pcs->K, pcs->ncode);
	}

	reset_Mu_non_seg_site(empcs);	/* Reset Mu for non-segregating sites. */
	initialize_count_Mu_X_and_gap(empcs);

	/* For labels, only owned pointers memory and need to be freed.
	 * This can be an independent object, em_phyclust_label/empcl. */
	initialize_em_phyclust_label(empcs, pcs);

	/* For sequencing error model. */
	empcs->se_type = pcs->se_type;
	initialize_em_phyclust_struct_se(empcs, pcs);

	return(empcs);
} /* End of initialize_em_phyclust_struct(). */

void free_em_phyclust_struct(em_phyclust_struct *empcs){
	/* For sequencing error model. */
	free_em_phyclust_struct_se(empcs);

	/* For labels, only owned pointers memory and need to be freed.
	 * This can be an independent object, em_phyclust_label/empcl. */
	free_em_phyclust_label(empcs);

	/* For EM. */
	free(empcs->n_class);
	free_int_RT(empcs->K, empcs->Mu);
	free_double_RT(empcs->N_X, empcs->Z_modified);
	free_double_RT(empcs->N_X, empcs->Z_normalized);
	free(empcs->Eta);
	free(empcs->log_Eta);
	free_int_RT_4D(empcs->N_X, empcs->K, empcs->ncode, empcs->count_Mu_X);
	if(empcs->gap_flag){
		free_int_RT_3D(empcs->N_X, empcs->K, empcs->count_Mu_X_gap);
	}

	free(empcs);
} /* End of free_em_phyclust_struct(). */

em_phyclust_struct* duplicate_em_phyclust_struct(em_phyclust_struct *org_empcs){
	em_phyclust_struct *new_empcs;
	int N_X = org_empcs->N_X, L = org_empcs->L, K = org_empcs->K;

	new_empcs = (em_phyclust_struct*) malloc(sizeof(em_phyclust_struct));
	new_empcs->code_type = org_empcs->code_type;
	new_empcs->ncode = org_empcs->ncode;
	new_empcs->gap_index = org_empcs->gap_index;
	new_empcs->gap_flag = org_empcs->gap_flag;
	new_empcs->N_X_org = org_empcs->N_X_org;
	new_empcs->N_X = N_X;
	new_empcs->N_seg_site = org_empcs->N_seg_site;
	new_empcs->L = L;
	new_empcs->K = K;
	new_empcs->X_org = org_empcs->X_org;
	new_empcs->X = org_empcs->X;
	new_empcs->map_X_org_to_X = org_empcs->map_X_org_to_X;
	new_empcs->map_X_to_X_org = org_empcs->map_X_to_X_org;
	new_empcs->replication_X = org_empcs->replication_X;
	new_empcs->seg_site_id = org_empcs->seg_site_id;
	new_empcs->class_id = org_empcs->class_id;

	/* For EM. */
	new_empcs->n_class = allocate_int_1D(K);
	new_empcs->Mu = allocate_int_RT(K, L); 
	new_empcs->Z_modified = allocate_double_RT(N_X, K); 
	new_empcs->Z_normalized = allocate_double_RT(N_X, K); 
	new_empcs->Eta = allocate_double_1D(K);
	new_empcs->log_Eta = allocate_double_1D(K);
	new_empcs->count_Mu_X = allocate_int_RT_4D(org_empcs->N_X, org_empcs->K, org_empcs->ncode, org_empcs->ncode);
	if(org_empcs->gap_flag){
		new_empcs->count_Mu_X_gap = allocate_int_RT_3D(org_empcs->N_X, org_empcs->K, org_empcs->ncode);
	}

	/* For labels. */
	duplicate_em_phyclust_label(org_empcs, new_empcs);

	/* For sequencing error models. */
	new_empcs->se_type = org_empcs->se_type;
	duplicate_em_phyclust_struct_se(org_empcs, new_empcs);

	Copy_empcs(org_empcs, new_empcs);
	return(new_empcs);
} /* End of duplicate_em_phyclust_struct(). */


/* Only SEMI and GENERAL should be dealed differently. */
void initialize_em_phyclust_label(em_phyclust_struct *empcs, phyclust_struct *pcs){
	int n_X, k, K = empcs->K;

	if(pcs->label->label_method == NONE){
		empcs->K_labeled = 0;
		empcs->N_X_labeled = 0;
		empcs->N_X_unlabeled = empcs->N_X;
		empcs->X_labeled = NULL;
		empcs->X_unlabeled = NULL;
		empcs->label_semi = NULL;
		empcs->label_index = NULL;
		empcs->Z_modified_labeled = NULL;
		empcs->Z_modified_unlabeled = NULL;
		empcs->Z_normalized_labeled = NULL;
		empcs->Z_normalized_unlabeled = NULL;
	} else{
		empcs->K_labeled = 0;
		empcs->N_X_labeled = pcs->label->N_index;
		empcs->N_X_unlabeled = empcs->N_X - empcs->N_X_labeled;
		empcs->X_labeled = allocate_int_2D_AP(empcs->N_X_labeled);
		empcs->X_unlabeled = allocate_int_2D_AP(empcs->N_X_unlabeled);
		empcs->label_semi = pcs->label->semi;
		empcs->label_index = pcs->label->index;
		empcs->Z_modified_labeled = allocate_double_2D_AP(empcs->N_X_labeled);
		empcs->Z_modified_unlabeled = allocate_double_2D_AP(empcs->N_X_unlabeled);
		empcs->Z_normalized_labeled = allocate_double_2D_AP(empcs->N_X_labeled);
		empcs->Z_normalized_unlabeled = allocate_double_2D_AP(empcs->N_X_unlabeled);
		reassign_label_pointer(empcs);

		for(n_X = 0; n_X < empcs->N_X_labeled; n_X++){
			/* Copy pre-assigned information (supervised) from pcl to empcl. */
			for(k = 0; k < K; k++){
				empcs->Z_normalized_labeled[n_X][k] = pcs->label->prob[n_X][k];
			}

			/* Find the total of labeled clusters. */
			if(empcs->label_semi[n_X] >= empcs->K_labeled){
				empcs->K_labeled = empcs->label_semi[n_X];
			}
		}
		empcs->K_labeled = empcs->K_labeled + 1;
	}
} /* End of initialize_em_phyclust_label(). */

void free_em_phyclust_label(em_phyclust_struct *empcs){
	if(empcs->N_X_labeled > 0){
		free(empcs->X_labeled);
		free(empcs->X_unlabeled);
		free(empcs->Z_modified_labeled);
		free(empcs->Z_modified_unlabeled);
		free(empcs->Z_normalized_labeled);
		free(empcs->Z_normalized_unlabeled);
	}
} /* End of free_em_phyclust_label(). */

void duplicate_em_phyclust_label(em_phyclust_struct *org_empcs, em_phyclust_struct *new_empcs){
	new_empcs->K_labeled = org_empcs->K_labeled;
	new_empcs->N_X_labeled = org_empcs->N_X_labeled;
	new_empcs->N_X_unlabeled = org_empcs->N_X_unlabeled;
	new_empcs->label_semi = org_empcs->label_semi;
	new_empcs->label_index = org_empcs->label_index;

	if(org_empcs->N_X_labeled == 0){
		new_empcs->X_labeled = NULL;
		new_empcs->X_unlabeled = NULL;
		new_empcs->Z_modified_labeled = NULL;
		new_empcs->Z_modified_unlabeled = NULL;
		new_empcs->Z_normalized_labeled = NULL;
		new_empcs->Z_normalized_unlabeled = NULL;
	} else{
		new_empcs->X_labeled = allocate_int_2D_AP(org_empcs->N_X_labeled);
		new_empcs->X_unlabeled = allocate_int_2D_AP(org_empcs->N_X_unlabeled);
		new_empcs->Z_modified_labeled = allocate_double_2D_AP(org_empcs->N_X_labeled);
		new_empcs->Z_modified_unlabeled = allocate_double_2D_AP(org_empcs->N_X_unlabeled);
		new_empcs->Z_normalized_labeled = allocate_double_2D_AP(org_empcs->N_X_labeled);
		new_empcs->Z_normalized_unlabeled = allocate_double_2D_AP(org_empcs->N_X_unlabeled);
		reassign_label_pointer(new_empcs);
	}
} /* End of duplicate_em_phyclust_label(). */


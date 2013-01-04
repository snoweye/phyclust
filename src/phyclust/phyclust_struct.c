/* This file contains functions for initialization and memory deallocation. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_struct.h"
#include "phyclust_em.h"
#include "phyclust_tool.h"
#include "phyclust_edist.h"
#include "phyclust_se_struct.h"


/* Initial a phyclust structure without assigning data X.
 * Assign X later by calling update_phyclust_struct() to update pcs.
 * Assign GAP_CODE by code_type.
 * Assign labeled data if any by calling update_phyclust_label().
 * Assign sequencing error if any by calling update_phyclust_struct_se().
 */
phyclust_struct* initialize_phyclust_struct(int code_type, int N_X_org, int L, int K){
	phyclust_struct *pcs;

	pcs = (phyclust_struct*) malloc(sizeof(phyclust_struct));
	pcs->code_type = code_type;
	pcs->ncode = NCODE[code_type];
	pcs->gap_index = GAP_INDEX[code_type];		/* For gaps. */
	pcs->gap_flag = 0;				/* Assigned by update_phyclust_struct(). */
	pcs->n_param = K - 1 + K * L;
	pcs->N_X_org = N_X_org;
	pcs->N_X = 0;					/* Assigned by update_phyclust_struct(). */
	pcs->N_seg_site = 0;				/* Assigned by update_phyclust_struct(). */
	pcs->L = L;
	pcs->K = K;
	pcs->X_org = allocate_int_2D_AP(N_X_org);	/* Assigned by ins or R. */
	pcs->X = NULL;					/* Assigned by update_phyclust_struct(). */
	pcs->map_X_org_to_X = NULL;			/* Assigned by update_phyclust_struct(). */
	pcs->map_X_to_X_org = NULL;			/* Assigned by update_phyclust_struct(). */
	pcs->replication_X = NULL;			/* Assigned by update_phyclust_struct(). */
	pcs->seg_site_id = NULL;			/* Assigned by update_phyclust_struct(). */
	pcs->Mu = allocate_int_RT(K, L); 
	pcs->Eta = allocate_double_1D(K); 
	/* For report, so use original dimensions. */
	pcs->Z_normalized = allocate_double_RT(N_X_org, K);
	
	pcs->logL_observed = 0.0;
	pcs->logL_entropy = 0.0;
	pcs->bic = 0.0;
	pcs->aic = 0.0;
	pcs->icl = 0.0;
	pcs->class_id = allocate_int_1D(N_X_org); 
	pcs->n_class = allocate_int_1D(K); 

	/* Labels. */
	pcs->label = NULL;				/* Point to update_phyclust_label(). */

	/* If code_type = NUCLEOTIDE & se_type == SE_YES, then this should be
	 * updated by updat_phyclust_struct_se(). */
	pcs->se_type = SE_NO;				/* Updated by outside function. */
	pcs->SE_P = NULL;				/* Assigned by update_phyclust_se_struct(). */

	return(pcs);
} /* End of initialize_phyclust_struct(). */

void free_phyclust_struct(phyclust_struct *pcs){
	free_phyclust_se_struct(pcs);
	free_phyclust_label(pcs->label);

	free(pcs->X_org);
	free(pcs->X);
	free(pcs->map_X_org_to_X);
	free(pcs->map_X_to_X_org);
	free(pcs->replication_X);
	free(pcs->seg_site_id);
	free_int_RT(pcs->K, pcs->Mu);
	free(pcs->Eta);
	free_double_RT(pcs->N_X_org, pcs->Z_normalized);
	free(pcs->class_id);
	free(pcs->n_class);
	free(pcs);
} /* End of free_phyclust_struct(). */

/* After assigning the data X_org, this function will extract the information
 * and update the pcs including: N_X, X, map_X_to_X_org, replication_X,
 * seg_site_id, and N_seg_site. */
void update_phyclust_struct(phyclust_struct *pcs){
	int i, n_X, l, flag, flag_gap, N_X_org = pcs->N_X_org, N_X, L = pcs->L;
	int map_X_to_X_org[N_X_org], replication_X[N_X_org], seg_site_id[L], N_seg_site = 0;

	pcs->map_X_org_to_X = allocate_int_1D(N_X_org);

	/* Assign N_X and map_X_to_X_org. */
	for(i = 0; i < N_X_org; i++){
		replication_X[i] = 0;
	}
	pcs->map_X_org_to_X[0] = 0;
	map_X_to_X_org[0] = 0;
	replication_X[0] = 1;
	N_X = 1;
	for(n_X = 1; n_X < N_X_org; n_X++){
		flag = 1;
		for(i = 0; i < N_X; i++){
			if(!(edist_D_HAMMING(L, pcs->X_org[n_X], pcs->X_org[map_X_to_X_org[i]]) > 0)){
				flag = 0;
				break;
			}
		}
		pcs->map_X_org_to_X[n_X] = i;
		if(flag){
			map_X_to_X_org[N_X++] = n_X;
		}
		replication_X[i]++;
	}
	/* Copy to pcs. */
	pcs->X = allocate_int_2D_AP(N_X);
	pcs->map_X_to_X_org = allocate_int_1D(N_X);
	pcs->replication_X = allocate_int_1D(N_X);
	for(i = 0; i < N_X; i++){
		pcs->X[i] = pcs->X_org[map_X_to_X_org[i]];
		pcs->map_X_to_X_org[i] = map_X_to_X_org[i];
		pcs->replication_X[i] = replication_X[i];
	}
	pcs->N_X = N_X;


	/* Assign gap_flag. */
	flag_gap = 0;
	for(n_X = 0; n_X < pcs->N_X; n_X++){
		for(l = 0; l < L; l++){
			if(pcs->X[n_X][l] == pcs->gap_index){
				flag_gap = 1;
				break;
			}
		}
		if(flag_gap){
			break;
		}
	}
	/* Copy to pcs. */
	pcs->gap_flag = flag_gap;


	/* Assign seg_site_id, and N_seg_site. If all are gap, set as seg_site. */
	for(l = 0; l < L; l++){
		flag = 0;
		flag_gap = 0;

		if(pcs->X[0][l] == pcs->gap_index || pcs->X[0][l] == MISSING_ALLELE){
			flag_gap = 1;
		}
		for(n_X = 1; n_X < pcs->N_X; n_X++){
			if(pcs->X[n_X][l] != pcs->X[0][l]){
				flag |= 1;
				break;
			}
			if(pcs->X[n_X][l] == pcs->gap_index || pcs->X[n_X][l] == MISSING_ALLELE){
				flag_gap++;
			}
		}

		if(flag || flag_gap == pcs->N_X){	/* All are gap, set as seg_site. */
			seg_site_id[N_seg_site++] = l;
		}
	}
	/* Copy to pcs. */
	pcs->seg_site_id = allocate_int_1D(N_seg_site);
	for(i = 0; i < N_seg_site; i++){
		pcs->seg_site_id[i] = seg_site_id[i];
	}
	pcs->N_seg_site = N_seg_site;

	/* Initialize pcl. */
	pcs->label = initialize_phyclust_label();
} /* End of update_phyclust_struct(). */




/* ----- For summary. ----- */
void assign_class(phyclust_struct *pcs){
	int n_X_org, k, tmp_k;
	double org_Z_normalized;

	for(k = 0; k < pcs->K; k++){
		pcs->n_class[k] = 0;
	}
	for(n_X_org = 0; n_X_org < pcs->N_X_org; n_X_org++){
		tmp_k = 0;
		org_Z_normalized = pcs->Z_normalized[n_X_org][0];
		for(k = 1; k < pcs->K; k++){
			if(org_Z_normalized < pcs->Z_normalized[n_X_org][k]){
				org_Z_normalized = pcs->Z_normalized[n_X_org][k];
				tmp_k = k;
			}
		}
		pcs->class_id[n_X_org] = tmp_k;
		pcs->n_class[tmp_k]++;
	}
} /* End of assign_class(). */

void update_ic(phyclust_struct *pcs, Q_matrix_array *QA){
	pcs->bic = -2 * pcs->logL_observed + (pcs->n_param + QA->total_n_param) * log(pcs->N_X_org);
	pcs->aic = -2 * pcs->logL_observed + 2 * (pcs->n_param + QA->total_n_param);
	pcs->icl = -2 * pcs->logL_entropy + (pcs->n_param + QA->total_n_param) * log(pcs->N_X_org);
} /* End of update_ic(). */




/* ----- For debug. ----- */
void print_X(phyclust_struct *pcs){
	int n_X, l;

	printf("X:\n");
	for(n_X = 0; n_X < pcs->N_X; n_X++){
		printf("    ");
		for(l = 0; l < pcs->L; l++){
		#if PRINT_CODE_TYPE == 0
			if(pcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_CODE[pcs->X[n_X][l]]);
			} else if(pcs->code_type == SNP){
				printf("%c ", SNP_CODE[pcs->X[n_X][l]]);
			}
		#else
			if(pcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_ID[pcs->X[n_X][l]]);
			} else if(pcs->code_type == SNP){
				printf("%c ", SNP_ID[pcs->X[n_X][l]]);
			}
		#endif
		}
		printf("\n");
	}
} /* End of print_X(). */

void print_Mu(phyclust_struct *pcs){
	int k, l;

	printf("Mu:\n");
	for(k = 0; k < pcs->K; k++){
		printf("    ");
		for(l = 0; l < pcs->L; l++){
		#if PRINT_CODE_TYPE == 0
			if(pcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_CODE[pcs->Mu[k][l]]);
			} else if(pcs->code_type == SNP){
				printf("%c ", SNP_CODE[pcs->Mu[k][l]]);
			}
		#else
			if(pcs->code_type == NUCLEOTIDE){
				printf("%c ", NUCLEOTIDE_ID[pcs->Mu[k][l]]);
			} else if(pcs->code_type == T_SNP){
				printf("%c ", SNP_ID[pcs->Mu[k][l]]);
			}
		#endif
		}
		printf("\n");
	}
} /* End of print_Mu(). */

void print_class_id(phyclust_struct *pcs){
	int n_X_org;
	printf("Class id:");
	for(n_X_org = 0; n_X_org < pcs->N_X_org; n_X_org++){
		printf(" %d", pcs->class_id[n_X_org]);
	}
	printf("\n");
} /* End of print_class_id(). */


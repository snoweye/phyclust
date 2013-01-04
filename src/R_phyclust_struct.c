/* This file contains functions for initialization and memory deallocation.
 * Copy from "phyclust_struct.c" to avoid coping some memory. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust/phyclust.h"


/* Initial a phyclust structure without assigning data X.
 * Assign X later and call update_phyclust_struct(). */
phyclust_struct* R_initialize_phyclust_struct(int code_type, int N_X_org, int L, int K){
	phyclust_struct *pcs = NULL;

	pcs = (phyclust_struct*) malloc(sizeof(phyclust_struct));
	pcs->code_type = code_type;
	pcs->ncode = NCODE[code_type];
	pcs->gap_index = GAP_INDEX[code_type];
	pcs->gap_flag = 0;
	pcs->n_param = K - 1 + K * L;
	pcs->N_X_org = N_X_org;
	pcs->N_X = 0;
	pcs->N_seg_site = 0;
	pcs->L = L;
	pcs->K = K;
	pcs->X_org = allocate_int_2D_AP(N_X_org);
	pcs->X = NULL;
	pcs->map_X_org_to_X = NULL;
	pcs->map_X_to_X_org = NULL;
	pcs->replication_X = NULL;
	pcs->seg_site_id = NULL;

	pcs->Mu = allocate_int_2D_AP(K);			/* Assigned by R. */
	pcs->Eta = NULL;					/* Assigned by R. */
	pcs->Z_normalized = allocate_double_2D_AP(N_X_org);	/* Assigned by R. */

	pcs->logL_observed = 0.0;
	pcs->logL_entropy = 0.0;
	pcs->bic = 0.0;
	pcs->aic = 0.0;
	pcs->icl = 0.0;
	pcs->class_id = NULL; 
	pcs->n_class = NULL; 

	/* Labels. */
	pcs->label = NULL;					/* Assigned by R_update_phyclust_label(pcs, R_label) */

	/* If code_type = NUCLEOTIDE & se_type == SE_YES, then this should be
	 * updated by updat_phyclust_struct_se(). */
	pcs->se_type = SE_NO;					/* Updated by outside function. */
	pcs->SE_P = NULL;					/* Assigned by update_phyclust_se_struct(). */

	return(pcs);
} /* End of R_initialize_phyclust_struct(). */


/* This function only frees pointers. The memory pointed by the pointer is
 * allocated in R and should not be free by R. */
void R_free_phyclust_struct(phyclust_struct *pcs){
	free(pcs->X_org);
	free(pcs->seg_site_id);
	free_phyclust_label(pcs->label);
	free(pcs->Mu);
	free(pcs->Z_normalized);
	free(pcs);
} /* End of free_phyclust_struct(). */


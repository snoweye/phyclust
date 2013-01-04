/* This file contains functions for initialization and memory deallocation.
 * Copy from "phyclust_label.c" to avoid coping some memory.
 *
 * Writen: Wei-Chen Chen on 2010/09/08. */


#include "R_phyclust.h"

/* Initial a phyclust label.
 *
 * label$semi or label (vector) will be pre-processed in R,
 * label$index: the indices (0, ..., N_X_org - 1) of labeled sequences,
 * label$prob: the probabilities for each sequence, dim = N_label * K.
 *
 * pcl should be initialized nomatter un- or semi-supervised
 * clustering, point to NULL if no labels.
 *
 * N_label, label_index and label_prob should be read only. */
void R_update_phyclust_label(phyclust_struct *pcs, SEXP R_label){
	/* Variables for preprocess R objects. */
	int label_method = NONE, N_label = 0, *label_semi = NULL, *label_index = NULL;
	double *tmp_ptr = NULL;
	 
	/* Update initial values. */
	if(R_label != R_NilValue){
		label_method = INTEGER(getListElement(R_label, "label.method"))[0];
		if(label_method != NONE){
			N_label = length(getListElement(R_label, "index"));
			label_semi = INTEGER(getListElement(R_label, "semi"));
			label_index = INTEGER(getListElement(R_label, "index"));
			tmp_ptr = REAL(getListElement(R_label, "prob"));
	       		/* For SEMI and GENERAL. */
			update_phyclust_label(pcs->label, label_method, N_label, label_semi, label_index, tmp_ptr,
				pcs->map_X_org_to_X, pcs->K);
		}
	}
} /* End of R_update_phyclust_label(). */


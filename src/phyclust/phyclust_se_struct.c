/* For sequencing error models. */

/* This file contains functions for initialization and memory deallocation. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_struct.h"
#include "phyclust_tool.h"
#include "phyclust_se_struct.h"
#include "phyclust_se_pmatrix.h"

void free_phyclust_se_struct(phyclust_struct *pcs){
	if(pcs->code_type == NUCLEOTIDE && pcs->se_type == SE_YES){
		free_SE_P_matrix(pcs->SE_P);
	}
} /* End of free_phyclust_se_struct(). */

void update_phyclust_se_struct(phyclust_struct *pcs, em_control *EMC){
	pcs->se_type = EMC->se_type;
	if(pcs->code_type == NUCLEOTIDE && EMC->se_type == SE_YES){
		pcs->SE_P = initialize_SE_P_matrix(pcs->code_type, EMC->se_model, EMC->se_constant, pcs->gap_flag, pcs->K);
		pcs->n_param = pcs->n_param + pcs->SE_P->n_param;
	}
} /* End of update_phyclust_se_struct(). */


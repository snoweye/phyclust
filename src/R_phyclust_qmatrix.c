/* This file contains updating functions for log transition probabilities. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust/phyclust.h"

/* Initial a Q matrix by given substitution model. */
Q_matrix* R_initialize_Q_matrix(int code_type, int substitution_model){
	Q_matrix *Q;

	Q = (Q_matrix*) malloc(sizeof(Q_matrix));
	Q->code_type = allocate_int_1D(1);
	Q->ncode = allocate_int_1D(1);
	Q->substitution_model = allocate_int_1D(1);
	Q->n_param = allocate_int_1D(1);
	Q->lower_bound = NULL;
	Q->upper_bound = NULL;

	Q->Pt = allocate_double_2D_AP(NCODE[code_type]);
	Q->log_Pt = allocate_double_2D_AP(NCODE[code_type]);
	Q->H = NULL;
	
	Q->pi = NULL;
	Q->kappa = NULL;
	Q->Tt = NULL;
	Q->check_param = NULL;

	*Q->code_type = code_type;
	*Q->ncode = NCODE[code_type];
	*Q->substitution_model = substitution_model;
	assign_FP_to_Q_matrix(substitution_model, Q);

	return(Q);
} /* End of R_initialize_Q_matrix(). */

void R_free_Q_matrix(Q_matrix *Q){
	free(Q->code_type);
	free(Q->substitution_model);
	free(Q->n_param);
	free(Q->Pt);
	free(Q->log_Pt);
	free(Q->ncode);
	free(Q);
} /* End of R_free_Q_matrix(). */


/* For sequencing error models. */

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
#include "phyclust_se_em.h"
#include "phyclust_se_convolution_logpL.h"


/* For each sequence, compute log_Pt for all (n, k) cells and save in empcs->Z_modified (log and unnormalized). */
void Update_Z_modified_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_modified[n_X][k] = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					empcs->Z_modified[n_X][k] += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
		}
	}
} /* End of Update_Z_modified_se_convolution(). */

/* For each sequence, compute log_Pt for all (n, k) cells and save in empcs->Z_modified (log and unnormalized). */
void Update_Z_modified_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < empcs->K; k++){
			empcs->Z_modified[n_X][k] = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					empcs->Z_modified[n_X][k] += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			s_to = empcs->ncode;	/* Observed gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				empcs->Z_modified[n_X][k] += empcs->SE_P->log_conv[k][s_from][s_to] *
					empcs->count_Mu_X_gap[n_X][k][s_from];
			}
		}
	}
} /* End of Update_Z_modified_gap_se_convolution(). */


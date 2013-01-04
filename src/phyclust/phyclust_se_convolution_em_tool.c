/* For sequencing error models. */

/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"
#include "phyclust_se_em.h"
#include "phyclust_se_pmatrix.h"


/* Independent summary tool. */
double LogL_observed_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;
	int K = empcs->K, flag_out_range;
	double logL_observed, a_Z_normalized[empcs->K], total_sum, scale_exp;

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < K; k++){
			a_Z_normalized[k] = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_normalized[k] += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
		}

		/* This function takes log of component densities "a_z_normalized",
		 * stores a total density in "total_sum",
		 * stores a scale for exponent in "scal_exp", and
		 * turn on "flag_out_range" if exponent is unstable. */
		e_step_with_stable_exp(&K, a_Z_normalized, &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			logL_observed += total_sum;
		} else{
			logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
	
	return(logL_observed);
} /* End of LogL_observed_se_convolution(). */

double LogL_observed_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA){
	int s_from, s_to, n_X, k;
	int K = empcs->K, flag_out_range;
	int fix_s_to = empcs->SE_P->ncode;	/* Observed gap. */
	double logL_observed, a_Z_normalized[empcs->K], total_sum, scale_exp;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	logL_observed = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		for(k = 0; k < K; k++){
			a_Z_normalized[k] = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_normalized[k] += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				a_Z_normalized[k] += empcs->SE_P->log_conv[k][s_from][fix_s_to] *
					empcs->count_Mu_X_gap[n_X][k][s_from];
			}
		}

		/* This function takes log of component densities "a_z_normalized",
		 * stores a total density in "total_sum",
		 * stores a scale for exponent in "scal_exp", and
		 * turn on "flag_out_range" if exponent is unstable. */
		e_step_with_stable_exp(&K, a_Z_normalized, &total_sum, &scale_exp, &flag_out_range);

		/* Update logL_observed. */
		total_sum = log(total_sum);
		if(flag_out_range){
			total_sum += scale_exp;
		}
		if(empcs->replication_X[n_X] == 1){
			logL_observed += total_sum;
		} else{
			logL_observed += total_sum * empcs->replication_X[n_X];
		}
	}
	
	return(logL_observed);
} /* End of LogL_observed_gap_se_convolution(). */


/* log complete-data likelihood, for debuging only. */
double LogL_complete_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){		/* QA_H unused */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_complete_se_convolution(). */

double LogL_complete_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){	/* QA_H unused */
	int s_from, s_to, n_X, k;
	int fix_s_to = empcs->SE_P->ncode;	/* Observed gap. */
	double logL_complete, total_sum, a_Z_modified;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = empcs->log_Eta[k];
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				a_Z_modified += empcs->SE_P->log_conv[k][s_from][fix_s_to] *
					empcs->count_Mu_X_gap[n_X][k][s_from];
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_complete_gap_se_convolution(). */


/* log complete-data profile likelihood, used in M-steps. */
double LogL_profile_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){		/* QA_H unused */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_profile_se_convolution(). */

double LogL_profile_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){	/* QA_H unused */
	int s_from, s_to, n_X, k;
	double logL_complete, total_sum, a_Z_modified;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	logL_complete = 0.0;
	for(n_X = 0; n_X < empcs->N_X; n_X++){
		total_sum = 0.0;
		for(k = 0; k < empcs->K; k++){
			a_Z_modified = 0.0;
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				for(s_to = 0; s_to < empcs->ncode; s_to++){
					a_Z_modified += empcs->SE_P->log_conv[k][s_from][s_to] *
						empcs->count_Mu_X[n_X][k][s_from][s_to];
				}
			}
			/* For gap. */
			s_to = empcs->SE_P->ncode;	/* Observed gap. */
			for(s_from = 0; s_from < empcs->ncode; s_from++){
				a_Z_modified += empcs->SE_P->log_conv[k][s_from][s_to] *
					empcs->count_Mu_X_gap[n_X][k][s_from];
			}
			total_sum += a_Z_modified * empcs->Z_normalized[n_X][k];
		}
		if(empcs->replication_X[n_X] == 1){
			logL_complete += total_sum;
		} else{
			logL_complete += total_sum * empcs->replication_X[n_X];
		}
	}

	return(logL_complete);
} /* End of LogL_profile_gap_se_convolution(). */




/* Utility functions unique for the sequencing error convolution model. */

/* The ret_pt is a double pointer point to a double array which has
 * the same K as QA->Q[k], and
 * the same row and column as empcs->SE_P->f_err.
 *
 * i.e. double array[K][NCODE][NCODE or NCODE_WIGAP].
 *
 * Return the log of convolution of mutation (Pt) and error (f_err) probabilities.
 *
 * This is a silly allocation for EE, but it is OK for EV, VE, and VV.
 * This should be in two different functions "_common" and "_split" as QA and
 * those will be pointed by a function pointer.
 */
void update_convolution_Pt_f_err(Q_matrix_array *QA, SE_P_matrix *SE_P){
	int k, s_from, s_to, s_between;
	double tmp_sum;

	if(QA->identifier == EE){
		for(s_from = 0; s_from < SE_P->ncode; s_from++){
			for(s_to = 0; s_to < SE_P->ncode; s_to++){
				tmp_sum = 0.0;
				for(s_between = 0; s_between < SE_P->ncode; s_between++){
					tmp_sum += QA->Q[0]->Pt[s_from][s_between] *
						SE_P->f_err[s_between][s_to];
				}
				SE_P->log_conv[0][s_from][s_to] = log(tmp_sum);
			}
		}
		for(k = 1; k < QA->K; k++){
			for(s_from = 0; s_from < SE_P->ncode; s_from++){
				for(s_to = 0; s_to < SE_P->ncode; s_to++){
					SE_P->log_conv[k][s_from][s_to] = SE_P->log_conv[0][s_from][s_to];
				}
			}
		}
	} else{
		for(k = 0; k < QA->K; k++){
			for(s_from = 0; s_from < SE_P->ncode; s_from++){
				for(s_to = 0; s_to < SE_P->ncode; s_to++){
					tmp_sum = 0.0;
					for(s_between = 0; s_between < SE_P->ncode; s_between++){
						tmp_sum += QA->Q[k]->Pt[s_from][s_between] *
							SE_P->f_err[s_between][s_to];
					}
					SE_P->log_conv[k][s_from][s_to] = log(tmp_sum);
				}
			}
		}
	}

	#if (EMDEBUG & 16) == 16
		print_convolution_Pt_f_err(SE_P->log_conv, QA->K, SE_P->ncode, SE_P->ncode);
	#endif

	#if (EMDEBUG & 32) == 32
		print_SE_P(SE_P);
	#endif
} /* End of update_convolution_Pt_f_err(). */

void update_convolution_Pt_f_err_gap(Q_matrix_array *QA, SE_P_matrix *SE_P){
	int k, s_from, s_to, s_between;
	double tmp_sum;

	if(QA->identifier == EE){
		for(s_from = 0; s_from < SE_P->ncode; s_from++){
			for(s_to = 0; s_to < SE_P->ncode_wigap; s_to++){
				tmp_sum = 0.0;
				for(s_between = 0; s_between < SE_P->ncode; s_between++){
					tmp_sum += QA->Q[0]->Pt[s_from][s_between] *
						SE_P->f_err[s_between][s_to];
				}
				SE_P->log_conv[0][s_from][s_to] = log(tmp_sum);
			}
		}
		for(k = 1; k < QA->K; k++){
			for(s_from = 0; s_from < SE_P->ncode; s_from++){
				for(s_to = 0; s_to < SE_P->ncode_wigap; s_to++){
					SE_P->log_conv[k][s_from][s_to] = SE_P->log_conv[0][s_from][s_to];
				}
			}
		}
	} else{
		for(k = 0; k < QA->K; k++){
			for(s_from = 0; s_from < SE_P->ncode; s_from++){
				for(s_to = 0; s_to < SE_P->ncode_wigap; s_to++){
					tmp_sum = 0.0;
					for(s_between = 0; s_between < SE_P->ncode; s_between++){
						tmp_sum += QA->Q[k]->Pt[s_from][s_between] *
							SE_P->f_err[s_between][s_to];
					}
					SE_P->log_conv[k][s_from][s_to] = log(tmp_sum);
				}
			}
		}
	}

	#if (EMDEBUG & 16) == 16
		print_convolution_Pt_f_err(SE_P->log_conv, QA->K, SE_P->ncode, SE_P->ncode_wigap);
	#endif

	#if (EMDEBUG & 32) == 32
		print_SE_P(SE_P);
	#endif
} /* End of update_convolution_Pt_f_err_gap(). */


void print_convolution_Pt_f_err(double ***log_conv, int K, int nrow, int ncol){
	int k, i, j;

	for(k = 0; k < K; k++){
		printf("k = %d\n", k);
		for(i = 0; i < nrow; i++){
			printf("  %c:", NUCLEOTIDE_CODE[i]);
			for(j = 0; j < ncol; j++){
				printf("  %.8f", log_conv[k][i][j]);
			}
			printf("\n");
		}
	}
} /* End of print_convolution_Pt_f_err(). */


/* For sequencing error models. */

/* This file contains updating functions for log sequencing error probabilities.
 * Copy from "phyclust_qmatrix.c" since the similarity of structure.
 * Probability f_err is defined rather than rates.
 * Currently, this model is for nucleotide data only and is unsupervised.
 *
 * WARNING: This is not a well optimized version and not quite extensible.
 *
 */


/* Without gap (row sum = 1, * sum = se_constant)
 *   A G C T
 * A ~ * * *
 * G * ~ * *
 * C * * ~ *
 * T * * * ~
 *
 * With gap (row sum = 1, * sum = se_constant)
 *   A G C T -
 * A ~ * * * *
 * G * ~ * * *
 * C * * ~ * *
 * T * * * ~ *
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_se_pmatrix.h"
#include "phyclust_tool.h"


/* Initial a SE_P matrix by given sequencing error model. */
SE_P_matrix* initialize_SE_P_matrix(int code_type, int se_model, double se_constant, int gap_flag, int K){
	SE_P_matrix *SE_P;
	int tmp_ncode;

	if(code_type != NUCLEOTIDE){
		fprintf_stderr("PE: The code_type is not supported except NUCLEOTIDE.\n");
		exit(1);
	}

	SE_P = (SE_P_matrix*) malloc(sizeof(SE_P_matrix));
	SE_P->code_type = code_type;
	SE_P->ncode = NCODE[code_type];
	SE_P->ncode_wigap = NCODE_WIGAP[code_type];
	SE_P->gap_index = GAP_INDEX[code_type];		/* For gaps. */
	SE_P->gap_flag = gap_flag;
	SE_P->se_model = se_model;

	assign_FP_to_SE_P_matrix(SE_P);

	SE_P->se_constant = se_constant;
	SE_P->lower_bound = 1e-16;
	SE_P->upper_bound = (SE_CONSTANT < 1.0 ? SE_CONSTANT : 1.0) - 1e-16;
	SE_P->lower_bound_diag = 1e-16;
	SE_P->upper_bound_diag = 1.0 - 1e-16;

	/* For f_err. */
	initialize_f_err(SE_P);
	SE_P->check_param = 1;

	/* For temporary storage only.
	 * Waste a few memory for the EE model, but OK for EV, VE, and VV.
	 * Require more codes as Q and QA to reduce memory and computing.
	 * update_convolution_Pt_f_err_*(...) need to be changed, too.
	 *
	 * pcs and empcs both have SE_P.
	 * Initial function update_phyclust_se_struct(...) for pcs will
	 * set K = 0. */
	SE_P->K = K;
	if(K > 0){
		tmp_ncode = gap_flag ? SE_P->ncode_wigap : SE_P->ncode;
		SE_P->log_conv = allocate_double_RT_3D(K, SE_P->ncode, tmp_ncode);
	}

	return(SE_P);
} /* End of initialize_SE_P_matrix(). */

void free_SE_P_matrix(SE_P_matrix *SE_P){
	free_double_RT(SE_P->ncode, SE_P->f_err);
	if(SE_P->K > 0){
		free_double_RT_3D(SE_P->K, SE_P->ncode, SE_P->log_conv);
	}
	free(SE_P);
} /* End of free_SE_P_matrix(). */

SE_P_matrix* duplicate_SE_P_matrix(SE_P_matrix *org_SE_P){
	SE_P_matrix *new_SE_P;

	new_SE_P = initialize_SE_P_matrix(org_SE_P->code_type, org_SE_P->se_model, org_SE_P->se_constant, org_SE_P->gap_flag, org_SE_P->K);
	new_SE_P->se_constant = org_SE_P->se_constant;
	new_SE_P->lower_bound = org_SE_P->lower_bound;
	new_SE_P->upper_bound = org_SE_P->upper_bound;
	new_SE_P->lower_bound_diag = org_SE_P->lower_bound_diag;
	new_SE_P->upper_bound_diag = org_SE_P->upper_bound_diag;
	copy_SE_P_matrix(org_SE_P, new_SE_P);

	/* No need for copy log_conv since it is a temporary storage. */

	return(new_SE_P);
} /* End of duplicate_SE_P_matrix(). */

void assign_FP_to_SE_P_matrix(SE_P_matrix *SE_P){
	switch(SE_P->se_model){
		case SE_CONVOLUTION:
			if(SE_P->gap_flag){
				SE_P->n_param = 15;
				SE_P->Check_param = &Check_param_f_err_se_convolution_gap;
				SE_P->Print_f_err = &Print_f_err_common_gap;
				SE_P->Convert_vect_to_f_err = &Convert_vect_to_f_err_se_convolution_gap;
				SE_P->Convert_f_err_to_vect = &Convert_f_err_to_vect_se_convolution_gap;
				SE_P->Copy_f_err = &Copy_f_err_common_gap;
			} else{
				SE_P->n_param = 11;
				SE_P->Check_param = &Check_param_f_err_se_convolution;
				SE_P->Print_f_err = &Print_f_err_common;
				SE_P->Convert_vect_to_f_err = &Convert_vect_to_f_err_se_convolution;
				SE_P->Convert_f_err_to_vect = &Convert_f_err_to_vect_se_convolution;
				SE_P->Copy_f_err = &Copy_f_err_common;
			}
			break;
		default:
			fprintf_stderr("PE: The SE_P model is not found.\n");
			exit(1);
			break;
	}
} /* End of assign_FP_to_SE_P_matrix(). */


void initialize_f_err(SE_P_matrix *SE_P){
	int i, j;
	int tmp_ncode = SE_P->gap_flag ? SE_P->ncode_wigap : SE_P->ncode;
	double tmp_prob, tmp_prob_err;

	switch(SE_P->se_model){
		case SE_CONVOLUTION:
			SE_P->f_err = allocate_double_RT(SE_P->ncode, tmp_ncode);
			tmp_prob = 1.0 - SE_P->se_constant / ((double) SE_P->ncode);
			tmp_prob_err = SE_P->se_constant / ((double) SE_P->ncode) / ((double) tmp_ncode - 1.0);
			for(i = 0; i < SE_P->ncode; i++){
				for(j = 0; j < tmp_ncode; j++){
					if(i != j){
						SE_P->f_err[i][j] = tmp_prob_err;
					} else{
						SE_P->f_err[i][j] = tmp_prob;
					}
				}
			}
			break;
		default:
			fprintf_stderr("PE: The SE_P model is not found.\n");
			exit(1);
			break;
	}
} /* End of initialize_f_err(). */




/* Negative of profile log-likelihood for minimizing. */
void Check_param_f_err_se_convolution(SE_P_matrix *SE_P){
	int i, j, flag = 1;
	double tmp_error = 0.0;

	for(i = 0; i < SE_P->ncode; i++){
		for(j = 0; j < SE_P->ncode; j++){
			if(i != j){
				flag = flag &&
					(SE_P->f_err[i][j] > SE_P->lower_bound) &&
					(SE_P->f_err[i][j] < SE_P->upper_bound);
				tmp_error = tmp_error + SE_P->f_err[i][j];
			} else{
				flag = flag &&
					(SE_P->f_err[i][j] > SE_P->lower_bound_diag) &&
					(SE_P->f_err[i][j] < SE_P->upper_bound_diag);
			}
		}
	}

	SE_P->check_param = flag;
} /* End of Check_param_f_err_se_convolution(). */

void Check_param_f_err_se_convolution_gap(SE_P_matrix *SE_P){
	int i, j, flag = 1;
	double tmp_error = 0.0;

	for(i = 0; i < SE_P->ncode; i++){
		for(j = 0; j < SE_P->ncode_wigap; j++){
			if(i != j){
				flag = flag &&
					(SE_P->f_err[i][j] > SE_P->lower_bound) &&
					(SE_P->f_err[i][j] < SE_P->upper_bound);
				tmp_error = tmp_error + SE_P->f_err[i][j];
			} else{
				flag = flag &&
					(SE_P->f_err[i][j] > SE_P->lower_bound_diag) &&
					(SE_P->f_err[i][j] < SE_P->upper_bound_diag);
			}
		}
	}

	SE_P->check_param = flag;
} /* End of Check_param_f_err_se_convolution_gap(). */




/* These functions will convert vectors to parameters. */
void Convert_vect_to_f_err_se_convolution(double *vect, SE_P_matrix *SE_P){
	double *tmp_vect = vect, tmp_sum, tmp_error;
	int i, j;

	tmp_error = 0.0;
	for(i = 0; i < SE_P->ncode - 1; i++){
		tmp_sum = 0.0;
		for(j = 0; j < SE_P->ncode; j++){
			if(i != j){
				SE_P->f_err[i][j] = *tmp_vect;
				tmp_sum += *tmp_vect;
				tmp_vect++;
			}
		}
		SE_P->f_err[i][i] = 1.0 - tmp_sum;
		tmp_error += tmp_sum;
	}

	/* i == SE_P->ncode - 1 */
	tmp_sum = 0.0;
	for(j = 0; j < SE_P->ncode - 2; j++){
		SE_P->f_err[i][j] = *tmp_vect;
		tmp_sum += *tmp_vect;
		tmp_vect++;
	}

	/* j == SE_P->ncode - 1. */
	tmp_error += tmp_sum;
	SE_P->f_err[i][j] = SE_P->se_constant - tmp_error;
	tmp_sum += SE_P->f_err[i][j];
	SE_P->f_err[i][i] = 1.0 - tmp_sum;

	SE_P->Check_param(SE_P);
} /* End of Convert_vect_to_f_err_se_convolution(). */

void Convert_vect_to_f_err_se_convolution_gap(double *vect, SE_P_matrix *SE_P){
	double *tmp_vect = vect, tmp_sum, tmp_error;
	int i, j;

	tmp_error = 0.0;
	for(i = 0; i < SE_P->ncode - 1; i++){
		tmp_sum = 0.0;
		for(j = 0; j < SE_P->ncode_wigap; j++){
			if(i != j){
				SE_P->f_err[i][j] = *tmp_vect;
				tmp_sum += *tmp_vect;
				tmp_vect++;
			}
		}
		SE_P->f_err[i][i] = 1.0 - tmp_sum;
		tmp_error += tmp_sum;
	}

	/* i == SE_P->ncode - 1 */
	tmp_sum = 0.0;
	for(j = 0; j < SE_P->ncode_wigap - 2; j++){
		SE_P->f_err[i][j] = *tmp_vect;
		tmp_sum += *tmp_vect;
		tmp_vect++;
	}

	/* j == SE_P->ncode_wigap - 1. */
	tmp_error += tmp_sum;
	SE_P->f_err[i][j + 1] = SE_P->se_constant - tmp_error;
	tmp_sum += SE_P->f_err[i][j + 1];
	SE_P->f_err[i][i] = 1.0 - tmp_sum;

	SE_P->Check_param(SE_P);
} /* End of Convert_vect_to_f_err_se_convolution_gap(). */




/* These functions will convert parameters to vectors. */
void Convert_f_err_to_vect_se_convolution(SE_P_matrix *SE_P, double *vect){
	double *tmp_vect = vect;
	int i, j;

	for(i = 0; i < SE_P->ncode - 1; i++){
		for(j = 0; j < SE_P->ncode; j++){
			if(i == j){
				continue;
			}
			*tmp_vect = SE_P->f_err[i][j];
			tmp_vect++;
		}
	}

	/* i == SE_P->ncode - 1 */
	for(j = 0; j < SE_P->ncode - 2; j++){
		*tmp_vect = SE_P->f_err[i][j];
		tmp_vect++;
	}
} /* End of Convert_f_err_to_vect_se_convolution(). */

void Convert_f_err_to_vect_se_convolution_gap(SE_P_matrix *SE_P, double *vect){
	double *tmp_vect = vect;
	int i, j;

	for(i = 0; i < SE_P->ncode - 1; i++){
		for(j = 0; j < SE_P->ncode_wigap; j++){
			if(i == j){
				continue;
			}
			*tmp_vect = SE_P->f_err[i][j];
			tmp_vect++;
		}
	}

	/* i == SE_P->ncode - 1 */
	for(j = 0; j < SE_P->ncode_wigap - 2; j++){
		*tmp_vect = SE_P->f_err[i][j];
		tmp_vect++;
	}
} /* End of Convert_f_err_to_vect_se_convolution_gap(). */




/* ----- For printing. ----- */
void Print_f_err_common(SE_P_matrix *SE_P){
	int i, j;
	double tmp_sum, total_error = 0.0;

	printf("SE_model: %s, n_param: %d\n", SE_MODEL[SE_P->se_model], SE_P->n_param);

	for(i = 0; i < SE_P->ncode; i++){
		printf("  p(.|%c):", NUCLEOTIDE_CODE[i]);
		tmp_sum = 0.0;
		for(j = 0; j < SE_P->ncode; j++){
			printf(" %.8f", SE_P->f_err[i][j]);
			tmp_sum = tmp_sum + SE_P->f_err[i][j];
			if(i != j){
				total_error = total_error + SE_P->f_err[i][j];
			}
		}
		printf("  sum = %.4f", tmp_sum);
		printf("\n");
	}

	printf("  total error = %.16f\n", total_error);
} /* End of Print_f_err_common(). */

void Print_f_err_common_gap(SE_P_matrix *SE_P){
	int i, j;
	double tmp_sum, total_error = 0.0;

	printf("SE_model: %s, n_param: %d\n", SE_MODEL[SE_P->se_model], SE_P->n_param);

	for(i = 0; i < SE_P->ncode; i++){
		printf("  p(.|%c):", NUCLEOTIDE_CODE[i]);
		tmp_sum = 0.0;
		for(j = 0; j < SE_P->ncode_wigap; j++){
			printf(" %.8f", SE_P->f_err[i][j]);
			tmp_sum = tmp_sum + SE_P->f_err[i][j];
			if(i != j){
				total_error = total_error + SE_P->f_err[i][j];
			}
		}
		printf("  sum = %.4f", tmp_sum);
		printf("\n");
	}

	printf("  total error = %.16f\n", total_error);
} /* End of Print_f_err_common_gap(). */




/* ----- For copy. ----- */
void Copy_f_err_common(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to){
	copy_double_RT(SE_P_from->ncode, SE_P_from->ncode, SE_P_from->f_err, SE_P_to->f_err);
} /* End of Copy_f_err_common(). */

void Copy_f_err_common_gap(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to){
	copy_double_RT(SE_P_from->ncode, SE_P_from->ncode_wigap, SE_P_from->f_err, SE_P_to->f_err);
} /* End of Copy_f_err_common_gap(). */


void copy_SE_P_matrix(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to){
	SE_P_to->Copy_f_err(SE_P_from, SE_P_to);
	SE_P_to->check_param = SE_P_from->check_param;
} /* End of copy_SE_P_matrix(). */

void reset_SE_P_matrix(SE_P_matrix *SE_P){
	SE_P_matrix *org_SE_P = initialize_SE_P_matrix(SE_P->code_type, SE_P->se_model, SE_P->se_constant, SE_P->gap_flag, SE_P->K);

	copy_SE_P_matrix(org_SE_P, SE_P);
	free_SE_P_matrix(org_SE_P);
} /* End of reset_SE_P_matrix(). */




/* ----- For debug. ----- */
void print_SE_P(SE_P_matrix *SE_P){
	SE_P->Print_f_err(SE_P);
} /* End of print_SE_P(). */



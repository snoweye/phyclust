/* This file contains updating functions for log transition probabilities. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_qmatrix.h"
#include "phyclust_tool.h"

/* Initial a Q matrix by given substitution model. */
Q_matrix* initialize_Q_matrix(int code_type, int substitution_model){
	int i;
	Q_matrix *Q;

	Q = (Q_matrix*) malloc(sizeof(Q_matrix));
	Q->code_type = allocate_int_1D(1);
	Q->ncode = allocate_int_1D(1);
	Q->substitution_model = allocate_int_1D(1);
	Q->n_param = allocate_int_1D(1);
	Q->lower_bound = allocate_double_1D(1);
	Q->upper_bound = allocate_double_1D(1);

	Q->Pt = allocate_double_SQ(NCODE[code_type]);
	Q->log_Pt = allocate_double_SQ(NCODE[code_type]);
	Q->H = allocate_double_1D(NCODE[code_type]);

	Q->pi = allocate_double_1D(NCODE[code_type]);
	Q->kappa = allocate_double_1D(1);
	Q->Tt = allocate_double_1D(1);
	Q->check_param = allocate_int_1D(1);

	*Q->code_type = code_type;
	*Q->ncode = NCODE[code_type];
	*Q->substitution_model = substitution_model;
	assign_FP_to_Q_matrix(substitution_model, Q);
	*Q->lower_bound = 1e-16;
	*Q->upper_bound = 1 - *Q->lower_bound;

	for(i = 0; i < NCODE[code_type]; i++){
		Q->pi[i] = 1.0 / (double) NCODE[code_type];
	}
	*Q->kappa = 1.0;
	*Q->Tt = 1.0;
	*Q->check_param = 1;

	return(Q);
} /* End of initialize_Q_matrix(). */

void free_Q_matrix(Q_matrix *Q){
	free(Q->code_type);
	free(Q->substitution_model);
	free(Q->n_param);
	free(Q->lower_bound);
	free(Q->upper_bound);
	free_double_RT(*Q->ncode, Q->Pt);
	free_double_RT(*Q->ncode, Q->log_Pt);
	free(Q->H);
	free(Q->ncode);
	free(Q->pi);
	free(Q->kappa);
	free(Q->Tt);
	free(Q->check_param);
	free(Q);
} /* End of free_Q_matrix(). */

Q_matrix* duplicate_Q_matrix(Q_matrix *org_Q){
	Q_matrix *new_Q;

	new_Q = initialize_Q_matrix(*org_Q->code_type, *org_Q->substitution_model);
	*new_Q->code_type = *org_Q->code_type;
	*new_Q->ncode = *org_Q->ncode;
	*new_Q->substitution_model = *org_Q->substitution_model;
	assign_FP_to_Q_matrix(*new_Q->substitution_model, new_Q);
	*new_Q->lower_bound = *org_Q->lower_bound;
	*new_Q->upper_bound = *org_Q->upper_bound;
	copy_Q_matrix(org_Q, new_Q);
	return(new_Q);
} /* End of duplicate_Q_matrix(). */

void assign_FP_to_Q_matrix(int substitution_model, Q_matrix *Q){
	switch(substitution_model){
		case JC69:
			*Q->n_param = 1;
			Q->Update_log_Pt = &Update_log_Pt_JC69;
			Q->Check_param = &Check_param_JC69;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_JC69;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_JC69;
			Q->Print_Q_matrix = &Print_Q_matrix_JC69;
			break;
		case K80:
			*Q->n_param = 2;
			Q->Update_log_Pt = &Update_log_Pt_K80;
			Q->Check_param = &Check_param_K80;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_K80;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_K80;
			Q->Print_Q_matrix = &Print_Q_matrix_K80;
			break;
		case F81:
			*Q->n_param = 4;
			Q->Update_log_Pt = &Update_log_Pt_F81;
			Q->Check_param = &Check_param_F81;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_F81;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_F81;
			Q->Print_Q_matrix = &Print_Q_matrix_F81;
			break;
		case HKY85:
			*Q->n_param = 5;
			Q->Update_log_Pt = &Update_log_Pt_HKY85;
			Q->Check_param = &Check_param_HKY85;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_HKY85;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_HKY85;
			Q->Print_Q_matrix = &Print_Q_matrix_HKY85;
			break;
		case SNP_JC69:
			*Q->n_param = 1;
			Q->Update_log_Pt = &Update_log_Pt_SNP_JC69;
			Q->Check_param = &Check_param_SNP_JC69;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_SNP_JC69;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_SNP_JC69;
			Q->Print_Q_matrix = &Print_Q_matrix_SNP_JC69;
			break;
		case SNP_F81:
			*Q->n_param = 2;
			Q->Update_log_Pt = &Update_log_Pt_SNP_F81;
			Q->Check_param = &Check_param_SNP_F81;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_SNP_F81;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_SNP_F81;
			Q->Print_Q_matrix = &Print_Q_matrix_SNP_F81;
			break;
		case E_F81:
			*Q->n_param = 1;
			Q->Update_log_Pt = &Update_log_Pt_F81;
			Q->Check_param = &Check_param_E_F81;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_E_F81;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_E_F81;
			Q->Print_Q_matrix = &Print_Q_matrix_F81;
			break;
		case E_HKY85:
			*Q->n_param = 2;
			Q->Update_log_Pt = &Update_log_Pt_HKY85;
			Q->Check_param = &Check_param_E_HKY85;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_E_HKY85;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_E_HKY85;
			Q->Print_Q_matrix = &Print_Q_matrix_HKY85;
			break;
		case E_SNP_F81:
			*Q->n_param = 1;
			Q->Update_log_Pt = &Update_log_Pt_SNP_F81;
			Q->Check_param = &Check_param_E_SNP_F81;
			Q->Convert_vect_to_Q_matrix = &Convert_vect_to_Q_matrix_E_SNP_F81;
			Q->Convert_Q_matrix_to_vect = &Convert_Q_matrix_to_vect_E_SNP_F81;
			Q->Print_Q_matrix = &Print_Q_matrix_SNP_F81;
			break;
		default:
			fprintf_stderr("PE: The substitution model is not found.\n");
			exit(1);
	}
} /* End of assign_FP_to_Q_matrix(). */

Q_matrix* repoint_Q_matrix(Q_matrix *Q_from){
	Q_matrix *Q_to;

	Q_to = (Q_matrix*) malloc(sizeof(Q_matrix));
	Q_to->code_type = Q_from->code_type;
	Q_to->ncode = Q_from->ncode;
	Q_to->substitution_model = Q_from->substitution_model;
	Q_to->n_param = Q_from->n_param;
	assign_FP_to_Q_matrix(*Q_to->substitution_model, Q_to);
	Q_to->lower_bound = Q_from->lower_bound;
	Q_to->upper_bound = Q_from->upper_bound;
	Q_to->Pt = Q_from->Pt;
	Q_to->log_Pt = Q_from->log_Pt;
	Q_to->H = Q_from->H;
	Q_to->pi = Q_from->pi;
	Q_to->kappa = Q_from->kappa;
	Q_to->Tt = Q_from->Tt;
	return(Q_to);
} /* End of repoint_Q_matrix(). */




/*   A G C T        A G C T
 * A . a a a      A . 1 1 1
 * G a . a a  or  G 1 . 1 1
 * C a a . a      C 1 1 . 1
 * T a a a .      T 1 1 1 .
 *
 * log_Pt[0]: x = y.
 * log_Pt[1]: x != y.
 * */
void Update_log_Pt_JC69(Q_matrix *Q){
	double A, Pt[2], log_Pt[2];

	A = exp(-4 * *Q->Tt);
	Pt[0] = 0.25 + 0.75 * A;
	Pt[1] = 0.25 - 0.25 * A;
	log_Pt[0] = log(Pt[0]);
	log_Pt[1] = log(Pt[1]);

	/* Compute transition probabilities. */
	Q->Pt[0][0] = Pt[0];
	Q->Pt[0][1] = Pt[1];
	Q->Pt[0][2] = Pt[1];
	Q->Pt[0][3] = Pt[1];
	Q->Pt[1][0] = Pt[1];
	Q->Pt[1][1] = Pt[0];
	Q->Pt[1][2] = Pt[1];
	Q->Pt[1][3] = Pt[1];
	Q->Pt[2][0] = Pt[1];
	Q->Pt[2][1] = Pt[1];
	Q->Pt[2][2] = Pt[0];
	Q->Pt[2][3] = Pt[1];
	Q->Pt[3][0] = Pt[1];
	Q->Pt[3][1] = Pt[1];
	Q->Pt[3][2] = Pt[1];
	Q->Pt[3][3] = Pt[0];

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log_Pt[0];
	Q->log_Pt[0][1] = log_Pt[1];
	Q->log_Pt[0][2] = log_Pt[1];
	Q->log_Pt[0][3] = log_Pt[1];
	Q->log_Pt[1][0] = log_Pt[1];
	Q->log_Pt[1][1] = log_Pt[0];
	Q->log_Pt[1][2] = log_Pt[1];
	Q->log_Pt[1][3] = log_Pt[1];
	Q->log_Pt[2][0] = log_Pt[1];
	Q->log_Pt[2][1] = log_Pt[1];
	Q->log_Pt[2][2] = log_Pt[0];
	Q->log_Pt[2][3] = log_Pt[1];
	Q->log_Pt[3][0] = log_Pt[1];
	Q->log_Pt[3][1] = log_Pt[1];
	Q->log_Pt[3][2] = log_Pt[1];
	Q->log_Pt[3][3] = log_Pt[0];

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_JC69(). */

/*   A G C T        A G C T
 * A . a b b      A . 1 k k
 * G a . b b  or  G 1 . k k
 * C b b . a      C k k . 1
 * T b b a .      T k k 1 .
 *
 * log_Pt[0]: x = y.
 * log_Pt[1]: transition. (A <-> G, C <-> T)
 * log_Pt[2]: transversion. (A,G <-> C,T)
 *   A G C T
 * A 0 1 2 2
 * G 1 0 2 2
 * C 2 2 0 1
 * T 2 2 1 0 */
void Update_log_Pt_K80(Q_matrix *Q){
	double exp_A, exp_B, Pt[3], log_Pt[3];

	exp_A = exp(-4 * *Q->Tt);
	exp_B = 2 * exp(-2 * (*Q->kappa + 1) * *Q->Tt);
	Pt[0] = (1 + exp_A + exp_B) * 0.25;
	Pt[1] = (1 + exp_A - exp_B) * 0.25;
	Pt[2] = (1 - exp_A) * 0.25;
	log_Pt[0] = log((1 + exp_A + exp_B) * 0.25);
	log_Pt[1] = log((1 + exp_A - exp_B) * 0.25);
	log_Pt[2] = log((1 - exp_A) * 0.25);

	/* Compute transition probabilities. */
	Q->Pt[0][0] = Pt[0];
	Q->Pt[0][1] = Pt[1];
	Q->Pt[0][2] = Pt[2];
	Q->Pt[0][3] = Pt[2];
	Q->Pt[1][0] = Pt[1];
	Q->Pt[1][1] = Pt[0];
	Q->Pt[1][2] = Pt[2];
	Q->Pt[1][3] = Pt[2];
	Q->Pt[2][0] = Pt[2];
	Q->Pt[2][1] = Pt[2];
	Q->Pt[2][2] = Pt[0];
	Q->Pt[2][3] = Pt[1];
	Q->Pt[3][0] = Pt[2];
	Q->Pt[3][1] = Pt[2];
	Q->Pt[3][2] = Pt[1];
	Q->Pt[3][3] = Pt[0];

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log_Pt[0];
	Q->log_Pt[0][1] = log_Pt[1];
	Q->log_Pt[0][2] = log_Pt[2];
	Q->log_Pt[0][3] = log_Pt[2];
	Q->log_Pt[1][0] = log_Pt[1];
	Q->log_Pt[1][1] = log_Pt[0];
	Q->log_Pt[1][2] = log_Pt[2];
	Q->log_Pt[1][3] = log_Pt[2];
	Q->log_Pt[2][0] = log_Pt[2];
	Q->log_Pt[2][1] = log_Pt[2];
	Q->log_Pt[2][2] = log_Pt[0];
	Q->log_Pt[2][3] = log_Pt[1];
	Q->log_Pt[3][0] = log_Pt[2];
	Q->log_Pt[3][1] = log_Pt[2];
	Q->log_Pt[3][2] = log_Pt[1];
	Q->log_Pt[3][3] = log_Pt[0];

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_K80(). */

/*      A    G    C    T
 * A    . pi_G pi_C pi_T
 * G pi_A    . pi_C pi_T
 * C pi_A pi_G    . pi_T
 * T pi_A pi_G pi_C    .
 */
void Update_log_Pt_F81(Q_matrix *Q){
	double pi_AG, pi_CT, exp_lambda_1_t, exp_lambda_2_t, exp_lambda_3_t,
		delta_AG, delta_CT, Delta_AGCT, Delta_CTAG, Delta,
		pi_A_Delta_CTAG, pi_G_Delta_CTAG, pi_C_Delta_AGCT, pi_T_Delta_AGCT,
		pi_A_delta_AG_exp_lambda_3_t, pi_G_delta_AG_exp_lambda_3_t,
		pi_C_delta_CT_exp_lambda_4_t, pi_T_delta_CT_exp_lambda_4_t;

	pi_AG = Q->pi[0] + Q->pi[1];
	pi_CT = Q->pi[2] + Q->pi[3];

	exp_lambda_1_t = 1;
	exp_lambda_2_t = exp(-*Q->Tt);
	exp_lambda_3_t = exp(-(pi_CT + pi_AG) * *Q->Tt);
	
	delta_AG = 1 / pi_AG;
	delta_CT = 1 / pi_CT;
	Delta_AGCT = exp_lambda_1_t + pi_AG / pi_CT * exp_lambda_2_t;
	Delta_CTAG = exp_lambda_1_t + pi_CT / pi_AG * exp_lambda_2_t;
	Delta = exp_lambda_1_t - exp_lambda_2_t;
	
	pi_A_Delta_CTAG = Q->pi[0] * Delta_CTAG;
	pi_G_Delta_CTAG = Q->pi[1] * Delta_CTAG;
	pi_C_Delta_AGCT = Q->pi[2] * Delta_AGCT;
	pi_T_Delta_AGCT = Q->pi[3] * Delta_AGCT;

	pi_A_delta_AG_exp_lambda_3_t = Q->pi[0] * delta_AG * exp_lambda_3_t;
	pi_G_delta_AG_exp_lambda_3_t = Q->pi[1] * delta_AG * exp_lambda_3_t;
	pi_C_delta_CT_exp_lambda_4_t = Q->pi[2] * delta_CT * exp_lambda_3_t;
	pi_T_delta_CT_exp_lambda_4_t = Q->pi[3] * delta_CT * exp_lambda_3_t;

	/* Compute transition probabilities. */
	Q->Pt[0][0] = pi_A_Delta_CTAG + pi_G_delta_AG_exp_lambda_3_t;
	Q->Pt[0][1] = pi_G_Delta_CTAG - pi_G_delta_AG_exp_lambda_3_t;
	Q->Pt[0][2] = Q->pi[2] * Delta;
	Q->Pt[0][3] = Q->pi[3] * Delta;
	Q->Pt[1][0] = pi_A_Delta_CTAG - pi_A_delta_AG_exp_lambda_3_t;
	Q->Pt[1][1] = pi_G_Delta_CTAG + pi_A_delta_AG_exp_lambda_3_t;
	Q->Pt[1][2] = Q->Pt[0][2];
	Q->Pt[1][3] = Q->Pt[0][3];
	Q->Pt[2][0] = Q->pi[0] * Delta;
	Q->Pt[2][1] = Q->pi[1] * Delta;
	Q->Pt[2][2] = pi_C_Delta_AGCT + pi_T_delta_CT_exp_lambda_4_t;
	Q->Pt[2][3] = pi_T_Delta_AGCT - pi_T_delta_CT_exp_lambda_4_t;
	Q->Pt[3][0] = Q->Pt[2][0];
	Q->Pt[3][1] = Q->Pt[2][1];
	Q->Pt[3][2] = pi_C_Delta_AGCT - pi_C_delta_CT_exp_lambda_4_t;
	Q->Pt[3][3] = pi_T_Delta_AGCT + pi_C_delta_CT_exp_lambda_4_t;

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log(Q->Pt[0][0]);
	Q->log_Pt[0][1] = log(Q->Pt[0][1]);
	Q->log_Pt[0][2] = log(Q->Pt[0][2]);
	Q->log_Pt[0][3] = log(Q->Pt[0][3]);
	Q->log_Pt[1][0] = log(Q->Pt[1][0]);
	Q->log_Pt[1][1] = log(Q->Pt[1][1]);
	Q->log_Pt[1][2] = Q->log_Pt[0][2];
	Q->log_Pt[1][3] = Q->log_Pt[0][3];
	Q->log_Pt[2][0] = log(Q->Pt[2][0]);
	Q->log_Pt[2][1] = log(Q->Pt[2][1]);
	Q->log_Pt[2][2] = log(Q->Pt[2][2]);
	Q->log_Pt[2][3] = log(Q->Pt[2][3]);
	Q->log_Pt[3][0] = Q->log_Pt[2][0];
	Q->log_Pt[3][1] = Q->log_Pt[2][1];
	Q->log_Pt[3][2] = log(Q->Pt[3][2]);
	Q->log_Pt[3][3] = log(Q->Pt[3][3]);

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_F81(). */

/*        A      G      C      T             A      G      C      T
 * A      . a*pi_G b*pi_C b*pi_T      A      . k*pi_G   pi_C   pi_T
 * G a*pi_A      . b*pi_C b*pi_T  or  G k*pi_A      .   pi_C   pi_T
 * C b*pi_A b*pi_G      . a*pi_T      C   pi_A   pi_G      . k*pi_T
 * T b*pi_A b*pi_G a*pi_C      .      T   pi_A   pi_G k*pi_C      .
 */
void Update_log_Pt_HKY85(Q_matrix *Q){
	double pi_AG, pi_CT, exp_lambda_1_t, exp_lambda_2_t, exp_lambda_3_t, exp_lambda_4_t,
		delta_AG, delta_CT, Delta_AGCT, Delta_CTAG, Delta,
		pi_A_Delta_CTAG, pi_G_Delta_CTAG, pi_C_Delta_AGCT, pi_T_Delta_AGCT,
		pi_A_delta_AG_exp_lambda_3_t, pi_G_delta_AG_exp_lambda_3_t,
		pi_C_delta_CT_exp_lambda_4_t, pi_T_delta_CT_exp_lambda_4_t;

	pi_AG = Q->pi[0] + Q->pi[1];
	pi_CT = Q->pi[2] + Q->pi[3];

	exp_lambda_1_t = 1;
	exp_lambda_2_t = exp(-*Q->Tt);
	exp_lambda_3_t = exp(-(pi_CT + *Q->kappa * pi_AG) * *Q->Tt);
	exp_lambda_4_t = exp(-(pi_CT * *Q->kappa + pi_AG) * *Q->Tt);
	
	delta_AG = 1 / pi_AG;
	delta_CT = 1 / pi_CT;
	Delta_AGCT = exp_lambda_1_t + pi_AG / pi_CT * exp_lambda_2_t;
	Delta_CTAG = exp_lambda_1_t + pi_CT / pi_AG * exp_lambda_2_t;
	Delta = exp_lambda_1_t - exp_lambda_2_t;
	
	pi_A_Delta_CTAG = Q->pi[0] * Delta_CTAG;
	pi_G_Delta_CTAG = Q->pi[1] * Delta_CTAG;
	pi_C_Delta_AGCT = Q->pi[2] * Delta_AGCT;
	pi_T_Delta_AGCT = Q->pi[3] * Delta_AGCT;

	pi_A_delta_AG_exp_lambda_3_t = Q->pi[0] * delta_AG * exp_lambda_3_t;
	pi_G_delta_AG_exp_lambda_3_t = Q->pi[1] * delta_AG * exp_lambda_3_t;
	pi_C_delta_CT_exp_lambda_4_t = Q->pi[2] * delta_CT * exp_lambda_4_t;
	pi_T_delta_CT_exp_lambda_4_t = Q->pi[3] * delta_CT * exp_lambda_4_t;

	/* Compute transition probabilities. */
	Q->Pt[0][0] = pi_A_Delta_CTAG + pi_G_delta_AG_exp_lambda_3_t;
	Q->Pt[0][1] = pi_G_Delta_CTAG - pi_G_delta_AG_exp_lambda_3_t;
	Q->Pt[0][2] = Q->pi[2] * Delta;
	Q->Pt[0][3] = Q->pi[3] * Delta;
	Q->Pt[1][0] = pi_A_Delta_CTAG - pi_A_delta_AG_exp_lambda_3_t;
	Q->Pt[1][1] = pi_G_Delta_CTAG + pi_A_delta_AG_exp_lambda_3_t;
	Q->Pt[1][2] = Q->Pt[0][2];
	Q->Pt[1][3] = Q->Pt[0][3];
	Q->Pt[2][0] = Q->pi[0] * Delta;
	Q->Pt[2][1] = Q->pi[1] * Delta;
	Q->Pt[2][2] = pi_C_Delta_AGCT + pi_T_delta_CT_exp_lambda_4_t;
	Q->Pt[2][3] = pi_T_Delta_AGCT - pi_T_delta_CT_exp_lambda_4_t;
	Q->Pt[3][0] = Q->Pt[2][0];
	Q->Pt[3][1] = Q->Pt[2][1];
	Q->Pt[3][2] = pi_C_Delta_AGCT - pi_C_delta_CT_exp_lambda_4_t;
	Q->Pt[3][3] = pi_T_Delta_AGCT + pi_C_delta_CT_exp_lambda_4_t;

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log(Q->Pt[0][0]);
	Q->log_Pt[0][1] = log(Q->Pt[0][1]);
	Q->log_Pt[0][2] = log(Q->Pt[0][2]);
	Q->log_Pt[0][3] = log(Q->Pt[0][3]);
	Q->log_Pt[1][0] = log(Q->Pt[1][0]);
	Q->log_Pt[1][1] = log(Q->Pt[1][1]);
	Q->log_Pt[1][2] = Q->log_Pt[0][2];
	Q->log_Pt[1][3] = Q->log_Pt[0][3];
	Q->log_Pt[2][0] = log(Q->Pt[2][0]);
	Q->log_Pt[2][1] = log(Q->Pt[2][1]);
	Q->log_Pt[2][2] = log(Q->Pt[2][2]);
	Q->log_Pt[2][3] = log(Q->Pt[2][3]);
	Q->log_Pt[3][0] = Q->log_Pt[2][0];
	Q->log_Pt[3][1] = Q->log_Pt[2][1];
	Q->log_Pt[3][2] = log(Q->Pt[3][2]);
	Q->log_Pt[3][3] = log(Q->Pt[3][3]);

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_HKY85(). */

/*   0 1        0 1
 * 0 . a  or  0 . 1  =  (1  1) * (0  0) * (0.5  0.5)
 * 1 a .      1 1 .     (1 -1)   (0 -2)   (0.5 -0.5)
 *
 * P(t) = (0.5+0.5*e^{-2t} 0.5-0.5*e^{-2t})
 *        (0.5-0.5*e^{-2t} 0.5+0.5*e^{-2t})
 *
 * log_Pt[0]: x = y.
 * log_Pt[1]: x != y.
 * */
void Update_log_Pt_SNP_JC69(Q_matrix *Q){
	double A, Pt[2], log_Pt[2];

	A = exp(-2 * *Q->Tt);
	Pt[0] = 0.5 + 0.5 * A;
	Pt[1] = 0.5 - 0.5 * A;
	log_Pt[0] = log(Pt[0]);
	log_Pt[1] = log(Pt[1]);

	/* Compute transition probabilities. */
	Q->Pt[0][0] = Pt[0];
	Q->Pt[0][1] = Pt[1];
	Q->Pt[1][0] = Pt[1];
	Q->Pt[1][1] = Pt[0];

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log_Pt[0];
	Q->log_Pt[0][1] = log_Pt[1];
	Q->log_Pt[1][0] = log_Pt[1];
	Q->log_Pt[1][1] = log_Pt[0];

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_SNP_JC69(). */

/*      0    1           0    1
 * 0    . pi_1  or  0    . pi_1  = (1      1/pi_0) * (0  0) * (        pi_0       1-pi_0)
 * 1 pi_0    .      1 pi_0    .    (1 -1/(1-pi_0))   (0 -1)   (pi_0(1-pi_0) pi_0(pi_0-1))
 *
 * P(t) = (pi_0+(1-pi_0)*e^{-t} (1-pi_0)-(1-pi_0)*e^{-t})
 *        (    pi_0-pi_0*e^{-t}     (1-pi_0)+pi_0*e^{-t})
 */
void Update_log_Pt_SNP_F81(Q_matrix *Q){
	double A;

	A = exp(-*Q->Tt);

	/* Compute transition probabilities. */
	Q->Pt[0][0] = Q->pi[0] + Q->pi[1] * A;
	Q->Pt[0][1] = Q->pi[1] - Q->pi[1] * A;
	Q->Pt[1][0] = Q->pi[0] - Q->pi[0] * A;
	Q->Pt[1][1] = Q->pi[1] + Q->pi[0] * A;

	/* Compute log transition probabilities. */
	Q->log_Pt[0][0] = log(Q->Pt[0][0]);
	Q->log_Pt[0][1] = log(Q->Pt[0][1]);
	Q->log_Pt[1][0] = log(Q->Pt[1][0]);
	Q->log_Pt[1][1] = log(Q->Pt[1][1]);

	/* Compute negative entropies. */
	Update_H(Q);
} /* End of Update_log_Pt_SNP_F81(). */




/* Compute negative entropies. */
void Update_H(Q_matrix *Q){
	int i, j;
	for(i = 0; i < *Q->ncode; i++){
		Q->H[i] = 0.0;
		for(j = 0; j < *Q->ncode; j++){
			Q->H[i] = Q->H[i] + Q->Pt[i][j] * Q->log_Pt[i][j];
		}
	}
} /* End of Update_H(). */




/* Negative of profile log-likelihood for minimizing.
 * Bvec stores parameters: eta_1, ..., eta_K, Q(substitution_model), Tt.
 *   JC69: Tt.
 *   K80: kappa, Tt.
 *   HKY85: pi_A, pi_C, pi_G, kappa, Tt.
 * vect stores parameters: Q(substitution_model), Tt.
 */
void Check_param_JC69(double *vect, Q_matrix *Q){
	*Q->check_param = (vect[0] > *Q->lower_bound);	/* Tt. */
} /* End of Check_param_JC69(). */

void Check_param_K80(double *vect, Q_matrix *Q){
	*Q->check_param =
		((vect[0] > *Q->lower_bound) &&		/* kappa. */
		 (vect[1] > *Q->lower_bound));		/* Tt. */
} /* End of Check_param_K80(). */

void Check_param_F81(double *vect, Q_matrix *Q){
	double pi_last = 1 - (vect[0] + vect[1] + vect[2]);
	*Q->check_param =
	       	((vect[0] > *Q->lower_bound) &&		/* pi_A > lower_bound. */
	       	 (vect[0] < *Q->upper_bound) &&		/* pi_A < upper_bound. */
		 (vect[1] > *Q->lower_bound) &&		/* pi_G > lower_bound. */
	       	 (vect[1] < *Q->upper_bound) &&		/* pi_G < upper_bound. */
		 (vect[2] > *Q->lower_bound) &&		/* pi_C > lower_bound. */
	       	 (vect[2] < *Q->upper_bound) &&		/* pi_C < upper_bound. */
		 (pi_last > *Q->lower_bound) &&		/* pi_T > lower_bound. */
		 (pi_last < *Q->upper_bound) &&		/* pi_T < upper_bound. */
		 (vect[3] > *Q->lower_bound) );		/* Tt. */
} /* End of Check_param_F81(). */

void Check_param_HKY85(double *vect, Q_matrix *Q){
	double pi_last = 1 - (vect[0] + vect[1] + vect[2]);
	*Q->check_param =
	       	((vect[0] > *Q->lower_bound) &&		/* pi_A > lower_bound. */
	       	 (vect[0] < *Q->upper_bound) &&		/* pi_A < upper_bound. */
		 (vect[1] > *Q->lower_bound) &&		/* pi_G > lower_bound. */
	       	 (vect[1] < *Q->upper_bound) &&		/* pi_G < upper_bound. */
		 (vect[2] > *Q->lower_bound) &&		/* pi_C > lower_bound. */
	       	 (vect[2] < *Q->upper_bound) &&		/* pi_C < upper_bound. */
		 (pi_last > *Q->lower_bound) &&		/* pi_T > lower_bound. */
		 (pi_last < *Q->upper_bound) &&		/* pi_T < upper_bound. */
		 (vect[3] > *Q->lower_bound) &&		/* kappa. */
		 (vect[4] > *Q->lower_bound) );		/* Tt. */
} /* End of Check_param_HKY85(). */

void Check_param_SNP_JC69(double *vect, Q_matrix *Q){
	*Q->check_param = (vect[0] > *Q->lower_bound);	/* Tt. */
} /* End of Check_param_SNP_JC69(). */

void Check_param_SNP_F81(double *vect, Q_matrix *Q){
	double pi_last = 1 - vect[0];
	*Q->check_param =
	       	((vect[0] > *Q->lower_bound) &&		/* pi_0 > lower_bound. */
	       	 (vect[0] < *Q->upper_bound) &&		/* pi_0 < upper_bound. */
		 (pi_last > *Q->lower_bound) &&		/* pi_1 > lower_bound. */
		 (pi_last < *Q->upper_bound) &&		/* pi_1 < upper_bound. */
		 (vect[1] > *Q->lower_bound) );		/* Tt. */
} /* End of Check_param_SNP_F81(). */

void Check_param_E_F81(double *vect, Q_matrix *Q){
	*Q->check_param = (vect[0] > *Q->lower_bound);	/* Tt. */
} /* End of Check_param_E_F81(). */

void Check_param_E_HKY85(double *vect, Q_matrix *Q){
	*Q->check_param =
	       	((vect[0] > *Q->lower_bound) &&		/* kappa. */
		 (vect[1] > *Q->lower_bound) );		/* Tt. */
} /* End of Check_param_E_HKY85(). */

void Check_param_E_SNP_F81(double *vect, Q_matrix *Q){
	*Q->check_param = (vect[0] > *Q->lower_bound);	/* Tt. */
} /* End of Check_param_E_SNP_F81(). */




/* These functions will convert vectors to parameters. */
void Convert_vect_to_Q_matrix_JC69(double *vect, Q_matrix *Q){
	*Q->Tt = vect[0];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_JC69(). */

void Convert_vect_to_Q_matrix_K80(double *vect, Q_matrix *Q){
	*Q->kappa = vect[0];
	*Q->Tt = vect[1];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_K80(). */

void Convert_vect_to_Q_matrix_F81(double *vect, Q_matrix *Q){
	Q->pi[0] = vect[0];
	Q->pi[1] = vect[1];
	Q->pi[2] = vect[2];
	Q->pi[3] = 1 - (Q->pi[0] + Q->pi[1] + Q->pi[2]);
	*Q->Tt = vect[3];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_F81(). */

void Convert_vect_to_Q_matrix_HKY85(double *vect, Q_matrix *Q){
	Q->pi[0] = vect[0];
	Q->pi[1] = vect[1];
	Q->pi[2] = vect[2];
	Q->pi[3] = 1 - (Q->pi[0] + Q->pi[1] + Q->pi[2]);
	*Q->kappa = vect[3];
	*Q->Tt = vect[4];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_HKY85(). */

void Convert_vect_to_Q_matrix_SNP_JC69(double *vect, Q_matrix *Q){
	*Q->Tt = vect[0];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_SNP_JC69(). */

void Convert_vect_to_Q_matrix_SNP_F81(double *vect, Q_matrix *Q){
	Q->pi[0] = vect[0];
	Q->pi[1] = 1 - Q->pi[0];
	*Q->Tt = vect[1];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_SNP_F81(). */

void Convert_vect_to_Q_matrix_E_F81(double *vect, Q_matrix *Q){
	*Q->Tt = vect[0];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_E_F81(). */

void Convert_vect_to_Q_matrix_E_HKY85(double *vect, Q_matrix *Q){
	*Q->kappa = vect[0];
	*Q->Tt = vect[1];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_E_HKY85(). */

void Convert_vect_to_Q_matrix_E_SNP_F81(double *vect, Q_matrix *Q){
	*Q->Tt = vect[0];
	Q->Check_param(vect, Q);
} /* End of Convert_vect_to_Q_matrix_E_SNP_F81(). */




/* These functions will convert parameters to vectors. */
void Convert_Q_matrix_to_vect_JC69(Q_matrix *Q, double *vect){
	vect[0] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_JC69(). */

void Convert_Q_matrix_to_vect_K80(Q_matrix *Q, double *vect){
	vect[0] = *Q->kappa;
	vect[1] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_K80(). */

void Convert_Q_matrix_to_vect_F81(Q_matrix *Q, double *vect){
	vect[0] = Q->pi[0];
	vect[1] = Q->pi[1];
	vect[2] = Q->pi[2];
	vect[3] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_F81(). */

void Convert_Q_matrix_to_vect_HKY85(Q_matrix *Q, double *vect){
	vect[0] = Q->pi[0];
	vect[1] = Q->pi[1];
	vect[2] = Q->pi[2];
	vect[3] = *Q->kappa;
	vect[4] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_HKY85(). */

void Convert_Q_matrix_to_vect_SNP_JC69(Q_matrix *Q, double *vect){
	vect[0] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_SNP_JC69(). */

void Convert_Q_matrix_to_vect_SNP_F81(Q_matrix *Q, double *vect){
	vect[0] = Q->pi[0];
	vect[1] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_SNP_F81(). */

void Convert_Q_matrix_to_vect_E_F81(Q_matrix *Q, double *vect){
	vect[0] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_E_F81(). */

void Convert_Q_matrix_to_vect_E_HKY85(Q_matrix *Q, double *vect){
	vect[0] = *Q->kappa;
	vect[1] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_E_HKY85(). */

void Convert_Q_matrix_to_vect_E_SNP_F81(Q_matrix *Q, double *vect){
	vect[0] = *Q->Tt;
} /* End of Convert_Q_matrix_to_vect_E_SNP_F81(). */




/* ----- For printing. ----- */
void Print_Q_matrix_JC69(Q_matrix *Q){
	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  Tt: %.8f\n", *Q->Tt);
} /* End of Print_Q_matrix_JC69(). */

void Print_Q_matrix_K80(Q_matrix *Q){
	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  kappa: %.8f, Tt: %.8f\n", *Q->kappa, *Q->Tt);
} /* End of Print_Q_matrix_K80(). */

void Print_Q_matrix_F81(Q_matrix *Q){
	int s;

	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  pi:");
	for(s = 0; s < *Q->ncode; s++){
		printf(" %.8f", Q->pi[s]);
	}
	printf("\n");
	printf("  Tt: %.8f\n", *Q->Tt);
} /* End of Print_Q_matrix_F81(). */

void Print_Q_matrix_HKY85(Q_matrix *Q){
	int s;

	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  pi:");
	for(s = 0; s < *Q->ncode; s++){
		printf(" %.8f", Q->pi[s]);
	}
	printf("\n");
	printf("  kappa: %.8f, Tt: %.8f\n", *Q->kappa, *Q->Tt);
} /* End of Print_Q_matrix_HKY85(). */

void Print_Q_matrix_SNP_JC69(Q_matrix *Q){
	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  Tt: %.8f\n", *Q->Tt);
} /* End of Print_Q_matrix_SNP_JC69(). */

void Print_Q_matrix_SNP_F81(Q_matrix *Q){
	int s;

	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  pi:");
	for(s = 0; s < NSNP; s++){
		printf(" %.8f", Q->pi[s]);
	}
	printf("\n");
	printf("  Tt: %.8f\n", *Q->Tt);
} /* End of Print_Q_matrix_SNP_F81(). */




/* ----- For copy. ----- */
void copy_Q_matrix(Q_matrix *Q_from, Q_matrix *Q_to){
	copy_double_RT(NCODE[*Q_from->code_type], NCODE[*Q_from->code_type], Q_from->Pt, Q_to->Pt);
	copy_double_RT(NCODE[*Q_from->code_type], NCODE[*Q_from->code_type], Q_from->log_Pt, Q_to->log_Pt);
	copy_double_1D(NCODE[*Q_from->code_type], Q_from->H, Q_to->H);
	copy_double_1D(NCODE[*Q_from->code_type], Q_from->pi, Q_to->pi);
	*Q_to->kappa = *Q_from->kappa;
	*Q_to->Tt = *Q_from->Tt;
	*Q_to->check_param = *Q_from->check_param;
} /* End of copy_Q_matrix(). */

void reset_Q_matrix(Q_matrix *Q){
	Q_matrix *org_Q = initialize_Q_matrix(*Q->code_type, *Q->substitution_model);
	copy_Q_matrix(org_Q, Q);
	free_Q_matrix(org_Q);
} /* End of reset_Q_matrix(). */




/* ----- For debug. ----- */
void print_log_Pt(Q_matrix *Q){
	int s_from, s_to;

	printf("log_Pt:\n");
	for(s_from = 0; s_from < NCODE[*Q->code_type]; s_from++){
		printf("    ");
		for(s_to = 0; s_to < NCODE[*Q->code_type]; s_to++){
			printf(" %f", Q->log_Pt[s_from][s_to]);
		}
		printf("\n");
	}
} /* End of print_log_Pt(). */

void print_Pt(Q_matrix *Q){
	int s_from, s_to;

	printf("Pt:\n");
	for(s_from = 0; s_from < NCODE[*Q->code_type]; s_from++){
		printf("    ");
		for(s_to = 0; s_to < NCODE[*Q->code_type]; s_to++){
			printf(" %f", Q->Pt[s_from][s_to]);
		}
		printf("\n");
	}
} /* End of print_Pt(). */

void print_H(Q_matrix *Q){
	int s_from;

	printf("H:\n");
	printf("    ");
	for(s_from = 0; s_from < NCODE[*Q->code_type]; s_from++){
		printf(" %f", Q->H[s_from]);
	}
	printf("\n");
} /* End of print_Pt(). */

void print_Q(Q_matrix *Q){
	int s;

	printf("Q_matrix: %s, n_param: %d\n", SUBSTITUTION_MODEL[*Q->substitution_model], *Q->n_param);
	printf("  pi:");
	for(s = 0; s < NCODE[*Q->code_type]; s++){
		printf(" %.8f", Q->pi[s]);
	}
	printf("\n");
	printf("  kappa: %.8f, Tt: %.8f\n", *Q->kappa, *Q->Tt);
} /* End of print_Q(). */


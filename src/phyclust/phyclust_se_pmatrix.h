/* For sequencing error models. */

/* This file contains declarations for phyclust.
 * Copy from "phyclust_qmatrix.h" since the similarity of structure.
 * Probability f_err is defined rather than rates.
 * Currently, this model is only for nucleotide data and unsupervised.
 */


#ifndef __PHYCLUST_SE_P_MATRIX_
#define __PHYCLUST_SE_P_MATRIX_

typedef struct _SE_P_matrix		SE_P_matrix;


/* Parameters of probability matrix in sequencing error model. */
struct _SE_P_matrix{
/* Fixed variables, used in initial and duplicate functions. */
	/* Define code type. */
	int	code_type;						/* NUCLEOTIDE/SNP. */
	int	ncode;							/* = NN. */
	int	ncode_wigap;					/* = NNG. */
	int	gap_index;						/* = NNG or NSNPG, indicates the gap index. */
	int	gap_flag;						/* = 0 or 1 for data without gap or with gap. */
	/* Configuration of P. */	
	int	se_model;						/* Sequencing error model. */
	int	n_param;						/* Number of parameters in f_err. */
	/* FP for f_err. */
	void	(*Check_param)(SE_P_matrix*);				/* Check paramaters. */
	void	(*Print_f_err)(SE_P_matrix*);				/* Print f_err matrix. */
	void	(*Convert_vect_to_f_err)(double*, SE_P_matrix*);	/* Convert vector to f_err. */
	void	(*Convert_f_err_to_vect)(SE_P_matrix*, double*);	/* Convert f_err to vector. */
	void	(*Copy_f_err)(SE_P_matrix*, SE_P_matrix*);		/* Copy f_err matrix. */
	/* For optermizing subroutine, restriction of f_err. */
	double	se_constant;						/* A constrain constant, dependent on models, default 1e-2. */
	double  lower_bound;						/* lower bound of parameters, default 1e-16. */
	double  upper_bound;						/* uper bound of parameters, default se_constant - 1e-16. */
	double  lower_bound_diag;					/* lower bound of diagonal parameters, default 1e-16. */
	double  upper_bound_diag;					/* uper bound of diagonal parameters, default 1 - 1e-16. */

/* Dynamical variables, used in copy functions only. */
	double	**f_err;						/* Sequencing error probabilities, dim = ncode * ncode_wigap. */
	int	check_param;						/* 0, fail. 1, reasonable parameter. */

	int	K;							/* Temporary storage for K. */
	double  ***log_conv;						/* Temporary storage for log of convolution. */
};

SE_P_matrix* initialize_SE_P_matrix(int code_type, int se_model, double se_constant, int gap_flag, int K);
void free_SE_P_matrix(SE_P_matrix *SE_P);
SE_P_matrix* duplicate_SE_P_matrix(SE_P_matrix *org_SE_P);
void assign_FP_to_SE_P_matrix(SE_P_matrix *SE_P);
void initialize_f_err(SE_P_matrix *SE_P);


/* These functions check vectors. */
void Check_param_f_err_se_convolution(SE_P_matrix *SE_P);
void Check_param_f_err_se_convolution_gap(SE_P_matrix *SE_P);

/* These functions convert a vector to parameters. */
void Convert_vect_to_f_err_se_convolution(double *vect, SE_P_matrix *SE_P);
void Convert_vect_to_f_err_se_convolution_gap(double *vect, SE_P_matrix *SE_P);

/* These functions convert parameters to a vector. */
void Convert_f_err_to_vect_se_convolution(SE_P_matrix *SE_P, double *vect);
void Convert_f_err_to_vect_se_convolution_gap(SE_P_matrix *SE_P, double *vect);

/* These functions print f_err matrix. */
void Print_f_err_common(SE_P_matrix *SE_P);
void Print_f_err_common_gap(SE_P_matrix *SE_P);


/* ----- For copy and dynamic update. ----- */
void Copy_f_err_common(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to);
void Copy_f_err_common_gap(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to);

void copy_SE_P_matrix(SE_P_matrix *SE_P_from, SE_P_matrix *SE_P_to);
void reset_SE_P_matrix(SE_P_matrix *SE_P);


/* ----- For debug. ----- */
void print_SE_P(SE_P_matrix *SE_P);

#endif	/* End of __PHYCLUST_SE_P_MATRIX_. */


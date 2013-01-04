/* This file contains declarations for phyclust. */


#ifndef __PHYCLUST_Q_MATRIX_
#define __PHYCLUST_Q_MATRIX_

typedef struct _Q_matrix		Q_matrix;


/* Parameters of transition rate matrix. */
struct _Q_matrix{
/* Fixed variables, used in initial and duplicate functions. */
	/* Define code type. */
	int	*code_type;					/* NUCLEOTIDE/SNP. */
	int	*ncode;						/* = NN or NSNP. */
	/* Configuration of Q. */	
	int	*substitution_model;				/* Substitution model. JC69, K80, or HKY85. */
	int	*n_param;					/* Number of parameters including Q and Tt. */
	/* FP for Q. */
	void	(*Update_log_Pt)(Q_matrix*);			/* A function point to log_Pt_<substitution_model>(). */
	void	(*Check_param)(double*, Q_matrix*);		/* Check paramaters. */
	void	(*Convert_vect_to_Q_matrix)(double*, Q_matrix*);	/* Convert vector to Q. */
	void	(*Convert_Q_matrix_to_vect)(Q_matrix*, double*);	/* Convert Q to vector. */
	void	(*Print_Q_matrix)(Q_matrix*);			/* Print Q and Tt. */
	/* For optermizing subroutine, restriction of pi, kappa, Tt. */
	double  *lower_bound;					/* lower bound of parameters, default 1e-10. */
	double  *upper_bound;					/* uper bound of parameters, default 1 - lower_bound. */
/* Dynamical variables, used in copy functions only. */
	double	**Pt;						/* Transition probabilities, dim = ncode * ncode. */
	double	**log_Pt;					/* Log of transition probabilities, dim = ncode * ncode. */
	double  *H;						/* Negative of entropy, rowSums(Pt * log_Pt). */
	double	*pi;						/* pi for A, G, C, and T, dim = ncode. */
	double	*kappa;						/* kappa. */
	double	*Tt;						/* Total evolution time. */
	int	*check_param;					/* 0, fail. 1, reasonable parameter. */
};

Q_matrix* initialize_Q_matrix(int code_type, int substitution_model);
void free_Q_matrix(Q_matrix *Q);
Q_matrix* duplicate_Q_matrix(Q_matrix *org_Q);
void assign_FP_to_Q_matrix(int substitution_model, Q_matrix *Q);
Q_matrix* repoint_Q_matrix(Q_matrix *Q_from);


/* These functions update log(P(t)). */
void Update_log_Pt_JC69(Q_matrix *Q);
void Update_log_Pt_K80(Q_matrix *Q);
void Update_log_Pt_F81(Q_matrix *Q);
void Update_log_Pt_HKY85(Q_matrix *Q);
void Update_log_Pt_SNP_JC69(Q_matrix *Q);
void Update_log_Pt_SNP_F81(Q_matrix *Q);

/* This function update H. */
void Update_H(Q_matrix *Q);

/* These functions check vectors. */
void Check_param_JC69(double *vect, Q_matrix *Q);
void Check_param_K80(double *vect, Q_matrix *Q);
void Check_param_F81(double *vect, Q_matrix *Q);
void Check_param_HKY85(double *vect, Q_matrix *Q);
void Check_param_SNP_JC69(double *vect, Q_matrix *Q);
void Check_param_SNP_F81(double *vect, Q_matrix *Q);
void Check_param_E_F81(double *vect, Q_matrix *Q);
void Check_param_E_HKY85(double *vect, Q_matrix *Q);
void Check_param_E_SNP_F81(double *vect, Q_matrix *Q);

/* These functions convert a vector to parameters. */
void Convert_vect_to_Q_matrix_JC69(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_K80(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_F81(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_HKY85(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_SNP_JC69(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_SNP_F81(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_E_F81(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_E_HKY85(double *vect, Q_matrix *Q);
void Convert_vect_to_Q_matrix_E_SNP_F81(double *vect, Q_matrix *Q);

/* These functions convert parameters to a vector. */
void Convert_Q_matrix_to_vect_JC69(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_K80(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_F81(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_HKY85(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_SNP_JC69(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_SNP_F81(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_E_F81(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_E_HKY85(Q_matrix *Q, double *vect);
void Convert_Q_matrix_to_vect_E_SNP_F81(Q_matrix *Q, double *vect);

/* These functions print Q matrix. */
void Print_Q_matrix_JC69(Q_matrix *Q);
void Print_Q_matrix_K80(Q_matrix *Q);
void Print_Q_matrix_F81(Q_matrix *Q);
void Print_Q_matrix_HKY85(Q_matrix *Q);
void Print_Q_matrix_SNP_JC69(Q_matrix *Q);
void Print_Q_matrix_SNP_F81(Q_matrix *Q);

/* ----- For copy. ----- */
void copy_Q_matrix(Q_matrix *Q_from, Q_matrix *Q_to);
void reset_Q_matrix(Q_matrix *Q);


/* ----- For debug. ----- */
void print_log_Pt(Q_matrix *Q);
void print_Pt(Q_matrix *Q);
void print_H(Q_matrix *Q);
void print_Q(Q_matrix *Q);

#endif	/* End of __PHYCLUST_Q_MATRIX_. */


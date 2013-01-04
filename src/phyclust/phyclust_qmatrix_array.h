/* This file contains declarations for phyclust. */


#ifndef __PHYCLUST_Q_MATRIX_ARRAY_
#define __PHYCLUST_Q_MATRIX_ARRAY_

#include "phyclust_qmatrix.h"

typedef struct _Q_matrix_array		Q_matrix_array;


/* Array of parameters of transition rate matrix. */
struct _Q_matrix_array{
/* Fixed variables in all Q_matrix. */
	/* Define code type. */
	int	code_type;		/* NUCLEOTIDE/SNP. */
	int	ncode;			/* = NN or NSNP. */
	/* Configuration of QA. */	
	int	K;			/* Number of clusters. */
	int	identifier;		/* EE, EV, VE, VV. */
	int	total_n_param;		/* Totoal number of parameters. */
	int	substitution_model;	/* Substitution model of Q_matrix. */
	int	n_param;		/* Number of parameters of Q_matrix. */
	int	check_param;		/* Check parameters of Q_matrix_array. */
	double	lower_bound;		/* Lower bound of pi. */
	double	upper_bound;		/* Upper bound of pi. */

	/* For computation used in em_step(). */
	void	(*Update_log_Pt)(Q_matrix_array*);		/* A function point to log_Pt_<substitution_model>(). */
	void	(*Check_param)(Q_matrix_array*);		/* Check paramaters. */
	void	(*Convert_vect_to_Q_matrix_array)(double*, Q_matrix_array*);	/* Convert vector to QA. */
	void	(*Convert_Q_matrix_array_to_vect)(Q_matrix_array*, double*);	/* Convert QA to vector. */
	void	(*Copy_Q_matrix_array)(Q_matrix_array*, Q_matrix_array*);

/* Memory storage. */
	Q_matrix	**Q;			/* Q_matrix array. */
	double		*tmp_vect;		/* For convert parameters between Q and vect. */
};

Q_matrix_array* initialize_Q_matrix_array(int code_type, int K, int substitution_model, int identifier);
void free_Q_matrix_array(Q_matrix_array *QA);
Q_matrix_array* duplicate_Q_matrix_array(Q_matrix_array *org_QA);

/* These functions update log(P(t)). */
void Update_log_Pt_common(Q_matrix_array *QA);
void Update_log_Pt_split(Q_matrix_array *QA);

/* These functions check vectors. */
void Check_param_common(Q_matrix_array *QA);
void Check_param_split(Q_matrix_array *QA);

/* These functions convert a vector to parameters. */
void Convert_vect_to_Q_matrix_array_EE(double *vect, Q_matrix_array *QA);
void Convert_vect_to_Q_matrix_array_EV(double *vect, Q_matrix_array *QA);
void Convert_vect_to_Q_matrix_array_VE(double *vect, Q_matrix_array *QA);
void Convert_vect_to_Q_matrix_array_VV(double *vect, Q_matrix_array *QA);

/* These functions convert parameters to a vector. */
void Convert_Q_matrix_array_to_vect_EE(Q_matrix_array *QA, double *vect);
void Convert_Q_matrix_array_to_vect_EV(Q_matrix_array *QA, double *vect);
void Convert_Q_matrix_array_to_vect_VE(Q_matrix_array *QA, double *vect);
void Convert_Q_matrix_array_to_vect_VV(Q_matrix_array *QA, double *vect);

/* ----- For copy. ----- */
void Copy_Q_matrix_array_EE(Q_matrix_array *QA_from, Q_matrix_array *QA_to);
void Copy_Q_matrix_array_EV(Q_matrix_array *QA_from, Q_matrix_array *QA_to);
void Copy_Q_matrix_array_VE(Q_matrix_array *QA_from, Q_matrix_array *QA_to);
void Copy_Q_matrix_array_VV(Q_matrix_array *QA_from, Q_matrix_array *QA_to);

void reset_Q_matrix_array(Q_matrix_array *QA);

/* ----- For debug. ----- */
void print_QA(Q_matrix_array *QA);

#endif	/* End of __PHYCLUST_Q_MATRIX_ARRAY_. */

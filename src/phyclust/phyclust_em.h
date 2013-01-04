/* This file contains declarations for em control, em phyclust struct, and
 * em function points. */


#ifndef __PHYCLUST_EM_
#define __PHYCLUST_EM_

#include "phyclust_struct.h"
#include "phyclust_qmatrix_array.h"

typedef struct _em_control		em_control;
typedef struct _em_phyclust_struct	em_phyclust_struct;
typedef struct _em_fp			em_fp;


/* Control em iterations and convergence in "phyclust_em_control.c". */
struct _em_control{
/* Fixed variables, used in initial and duplicate functions. */
	/* For em step. */
	int	exhaust_iter;			/* Exhaust em iterations. */
	int	fixed_iter;			/* Fixed em iterations for em+.EM. */
	int	short_iter;			/* Short em iterations. */
	int	EM_iter;			/* Long EM iterations. */
	double	short_eps;			/* Short em tolerance. */
	double	EM_eps;				/* Long EM tolerance. */

	/* For cm optim. */
	double	cm_reltol;			/* Relative tolerance. */
	int	cm_maxit;			/* Maximum iterations. */

	/* For NM optim. */
	double	nm_abstol_Mu_given_QA;		/* Absolute tolerance, update_flag = 0. */
	double	nm_abstol_QA_given_Mu;		/* Absolute tolerance, update_flag = 1. */
	double	nm_reltol_Mu_given_QA;		/* Relative tolerance. */
	double	nm_reltol_QA_given_Mu;		/* Relative tolerance. */
	int	nm_maxit_Mu_given_QA;		/* Maximum iterations. */
	int	nm_maxit_QA_given_Mu;		/* Maximum iterations. */
	/* For simplified optim. */
	int	est_non_seg_site;		/* Estimate non-segregating sites. */

	/* For initialization and EM steps. */
	int	max_init_iter;			/* Max initial steps. */
	int	min_n_class;			/* Min n for each class. */
	int	init_procedure;			/* Initialization procedure. */
	int	init_method;			/* Initialization method. */
	int	substitution_model;		/* Substition model. */
	int	edist_model;			/* Distance model for initialization. */
	int	identifier;			/* EE, EV, VE, VV. */
	int	code_type;			/* NUCLEOTIDE/SNP. */
	int	em_method;			/* EM method. */
	int	boundary_method;		/* Method to deal with boundary solutions. */
	double	Eta_lower_bound;		/* Lower bound of Eta, 1/N_X_org. */
	double	Eta_upper_bound;		/* Upper bound of Eta, 1 - Eta_lower_bound. */

/* Dynamical variables, used in copy functions only. */
	/* For report, reset by reset_em_control(). */
	double	converge_eps;			/* Convergent tolerance. */
	double	converge_error;			/* Convergent error as likelihood decreasing. */
	int	converge_flag;			/* Convergent flag. */
	int	converge_iter;			/* Convergent iteration. */
	int	converge_inner_iter;		/* Convergent inner iteration. */
	int	converge_cm_iter;		/* Convergent CM iteration. */
	int	update_flag;			/* 0 for update Mu, 1 for update Q. */

/* Extension applications. */
	/* For sequencing error models. */
	int	se_type;			/* For sequencing error. */
	int	se_model;			/* Sequencing error model. */
	double	se_constant;			/* constanf for SE_CONSTANT probability. */
}; /* End of struct _em_control. */

em_control* initialize_em_control();
void free_em_control(em_control *EMC);
em_control* duplicate_em_control(em_control *org_EMC);
void update_em_control(em_control *EMC);	/* update EMC when initial settings are changed. */
void reset_em_control(em_control *EMC);		/* For long EM steps. */


/* EM phyloclustering structure for EM steps only. */
struct _em_phyclust_struct{
/* Fixed variables, used in initial and duplicate functions, copy or point to pcs. */
	/* Define code type. */
	int		code_type;		/* NUCLEOTIDE/SNP. */
	int		ncode;			/* = NN or NSNP, indicates the dimension. */
	int		gap_index;		/* = NNG or NSNPG, indicates the gap index. */
	int		gap_flag;		/* = 0 or 1 for data without gap or with gap. */
	/* Constant points to phylcust_struct *pcs. */
	int		N_X_org;		/* Number of original sequences. */
	int		N_X;			/* Number of sequences. */
	int		N_seg_site;		/* Total segregating sites. */
	int		L;			/* Number of loci. */
	int		K;			/* Number of clusters. */
	int		**X_org;		/* Original data pointer, dim = N_X_org * L. */
	int		**X;			/* Data pointer, dim = N_X * L if compress = 1, dim = N_X_org * L otherwise. */
	int		*map_X_org_to_X;	/* Map indexes from X_org to X, dim = N_X_org. */
	int		*map_X_to_X_org;	/* Map indexes from X to X_org, dim = N_X. */
	int		*replication_X;		/* Count replications of each unique sequence, dim = N_X. */
	/* Summarized information for X. */
	int		*seg_site_id;		/* Segregating site id, dim = N_seg_site. */
	int		*class_id;		/* For manually initialization, dim = N_X_org. */

/* Dynamical variables, used in copy functions only, owned memory and need to be freed. */
	/* For EM. */
	int		*n_class;		/* Number of sequences for each class, dim = K. */
	int		**Mu;			/* Centers, dim = K * L. */
	double		**Z_modified;		/* Unnormalized log Z, dim = N_X * K. */
	double		**Z_normalized;		/* Normalized Z, dim = N_X * K. */
	double		*Eta;			/* Proportion, dim = K. */
	double		*log_Eta;		/* Log of proportion, dim = K. */
	double		logL_observed;		/* Observed logL. */

	/* Computing storage. */
	int		****count_Mu_X;		/* Used in update_Z_modified(), logL_observed(), initialize_count_Mu_X_and_gap(). */
	int		***count_Mu_X_gap;	/* Used in update_Z_modified(), logL_observed(), initialize_count_Mu_X_and_gap(). */

/* Extension applications. */
	/* For labels, only owned pointers memory and need to be freed.
	 * This can be an independent object, em_phyclust_label. */
	int		K_labeled;			/* Total of clusters with labeled sequences. */
	int		N_X_labeled;			/* Total of unique labeled sequences. */
	int		N_X_unlabeled;			/* Total of unique unlabeled sequences. */
	int		**X_labeled;			/* (Pointer) Point to empcs->X, length = N_X_labeled. */
	int		**X_unlabeled;			/* (Pointer) Point to empcs->X, length = N_X_unlabeled. */
	int		*label_semi;			/* (Pointer) Point to pcl->label->smple, length = N_X. */
	int		*label_index;			/* (Pointer) Point to pcl->label->index, length = N_X_labeled. */
	double		**Z_modified_labeled;		/* (Pointer) Point to empcs->Z_modified, dim = N_X_labeled. */
	double		**Z_modified_unlabeled;		/* (Pointer) Point to empcs->Z_modified, dim = N_X_unlabeled. */
	double		**Z_normalized_labeled;		/* (Pointer) Point to empcs->Z_normalized, dim = N_X_labeled. */
	double		**Z_normalized_unlabeled;	/* (Pointer) Point to empcs->Z_normalized, dim = N_X_unlabeled. */

	/* For sequencing error model. */
	int		se_type;			/* For sequencing error types. */
	double		***Z_err_modified;		/* Unnormalized log Z_err, dim = K * L * NCODE_WIGAP[code_type]. */
	double		***Z_err_normalized;		/* Normalized Z_err, dim = K * L * NCODE_WIGAP[code_type]. */
	double		*Eta_err;			/* Proportion, dim = 2. */
	double		*log_Eta_err;			/* Proportion, dim = 2. */
	SE_P_matrix	*SE_P;				/* Parameters of sequencing error models. */
}; /* End of struct _em_phyclust_struct. */

em_phyclust_struct* initialize_em_phyclust_struct(phyclust_struct *pcs);
void free_em_phyclust_struct(em_phyclust_struct *empcs);
em_phyclust_struct* duplicate_em_phyclust_struct(em_phyclust_struct *org_empcs);

void initialize_em_phyclust_label(em_phyclust_struct *empcs, phyclust_struct *pcs);
void free_em_phyclust_label(em_phyclust_struct *empcs);
void duplicate_em_phyclust_label(em_phyclust_struct *org_empcs, em_phyclust_struct *new_empcs);


/* E_step's. */
void e_step_with_stable_exp(int *K, double *a_Z_normalized, double *total_sum, double *scale_exp, int *flag_out_range);
void E_step_simple(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP);
void E_step_logL_observed(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP);
void E_step_logL_observed_label_semi(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP);
void E_step_logL_observed_label_general(em_phyclust_struct *empcs, Q_matrix_array *QA, em_fp *EMFP);

/* M_step's. */
int M_step_simple(em_phyclust_struct *empcs, Q_matrix_array *new_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);		/* tmp_QA unused */
/* All CM steps require to set EMC->update_flag = 1. */
int M_step_CM(em_phyclust_struct *empcs, Q_matrix_array *new_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);
int M_step_ACM(em_phyclust_struct *empcs, Q_matrix_array *new_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP,
		em_phyclust_struct *tmp_empcs, Q_matrix_array *tmp_QA);		/* tmp_QA unused */

/* Method to update Eta.
 * IGNORE = ignore the results and stop if Eta's are out of boundary.
 * ADJUST = adjust the Eta's to be fixed on the boundary. */
int Update_Eta_given_Z_ADJUST(em_phyclust_struct *empcs, em_control *EMC);
int Update_Eta_given_Z_IGNORE(em_phyclust_struct *empcs, em_control *EMC);

/* This special function updates empcs->Z_modified (log and unnormalized). */
void Update_Z_modified(em_phyclust_struct *empcs, Q_matrix_array *QA);

/* Check convergence. */
int Check_convergence_em(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, Q_matrix_array *new_QA,
		Q_matrix_array *org_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP);
int Check_convergence_org(em_phyclust_struct *new_empcs, em_phyclust_struct *org_empcs, Q_matrix_array *new_QA,
		Q_matrix_array *org_QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP);

/* Different EM methods. */
void Em_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
void Short_em_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);


/* Control function pointers in "phyclust_em_fp.c". */
struct _em_fp{
/* Deal with different initialization methods. */
	/* In "phyclust_init_method.c". */
	int	(*Update_init)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);	/* Initialization method. */

/* Deal with different EM algorithms. */
	/* In "phyclust_em_step.c". */
	int	(*M_step)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*, em_control*,
		       	em_fp*, em_phyclust_struct*, Q_matrix_array*);				/* M_step. */
	int	(*Check_convergence)(em_phyclust_struct*, em_phyclust_struct*, Q_matrix_array*,
		       	Q_matrix_array*, Q_matrix_array*, em_control*, em_fp*);			/* Check convergence. */
	void	(*Em_step)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);		/* Em_step. */
	void	(*Short_em_step)(em_phyclust_struct*, Q_matrix_array*, em_control*, em_fp*);	/* Short_em_step. */
	void	(*E_step_logL_observed)(em_phyclust_struct*, Q_matrix_array*, em_fp*);		/* LogL by E_step. */
	int	(*Update_Eta_given_Z)(em_phyclust_struct*, em_control*);			/* Update Eta given Z. */
	void	(*Update_Z_modified)(em_phyclust_struct*, Q_matrix_array*);			/* Updae Z_modified. */

	/* In "phyclust_logpL.h". */
	int	(*Maximize_logpL)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*,
		       	em_control*, em_fp*);							/* Maximize profile logL. */

	/* In "phyclust_em_tool.h". */
	double	(*LogL_observed)(em_phyclust_struct*, Q_matrix_array*);				/* Observed logL. */
	double	(*LogL_complete)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*);	/* Complete logL. */
	double	(*LogL_profile)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*);		/* Profile logL. */
	void	(*Copy_empcs)(em_phyclust_struct*, em_phyclust_struct*);			/* copy empcs to empcs. */
	void	(*Copy_empcs_to_pcs)(em_phyclust_struct*, phyclust_struct*);			/* copy empcs to pcs. */
	void	(*Copy_pcs_to_empcs)(phyclust_struct*, em_phyclust_struct*);			/* copy empcs to pcs. */

/* Deal with gap or non-gap sequences. */
	/* In "phyclust_logpL.c". */
	void	(*Update_Mu_given_QA)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*);	/* Update Mu given QA in logpL. */
	double	(*Compute_R)(em_phyclust_struct*, Q_matrix_array*, Q_matrix_array*);		/* Update function R(). */
}; /* End of struct _em_fp. */

em_fp* initialize_em_fp(em_control *EMC, phyclust_struct *pcs);
void free_em_fp(em_fp *EMFP);

#endif	/* End of __PHYCLUST_EM_. */

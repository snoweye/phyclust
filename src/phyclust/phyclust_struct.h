/* This file contains declarations for phyclust. */

#ifndef __PHYCLUST_STRUCT_
#define __PHYCLUST_STRUCT_

#include "phyclust_qmatrix_array.h"
#include "phyclust_label.h"
#include "phyclust_se_pmatrix.h"

typedef struct _phyclust_struct		phyclust_struct;


/* Phyloclustering structure for storing data and results. */
struct _phyclust_struct{
/* Fixed variables, used in initial and duplicate functions. */
	/* Define code type. */
	int		code_type;		/* NUCLEOTIDE/SNP. */
	int		ncode;			/* = NN or NSNP, indicates the dimension. */
	int		gap_index;		/* = NNG or NSNPG, indicates the gapgap index. */
	int		gap_flag;		/* = 0 or 1 for data without gap or with gap. */
	/* Storage. */
	int		n_param;		/* Number of parameters, including Mu's, and eta's. */
	int		N_X_org;		/* Number of original sequences. */
	int		N_X;			/* Number of sequences. */
	int		N_seg_site;		/* Total segregating sites, including all are gap. */
	int		L;			/* Number of loci. */
	int		K;			/* Number of clusters. */
	int		**X_org;		/* Original data pointer, dim = N_X_org * L. */
	int		**X;			/* Data pointer, dim = N_X * L. */
	int		*map_X_org_to_X;	/* Map indexes from X_org to X, dim = N_X_org. */
	int		*map_X_to_X_org;	/* Map indexes from X to X_org, dim = N_X. */
	int		*replication_X;		/* Count replications of each unique sequence, dim = N_X. */
	/* Summarizd information for X. */
	int		*seg_site_id;		/* Segregating site id, dim = N_seg_site. */

/* Dynamical variables, used in copy functions only. */
	/* Parameters. */
	int		**Mu;			/* Centers, dim = K * L. */
	double		*Eta;			/* Proportion, dim = K. */
	double		**Z_normalized;		/* Normalized Z, dim = N_X * K. */
	/* For summary. */
	double		logL_observed;		/* Observed logL. */
	double		logL_entropy;		/* Observed logL + entropy. */
	double		bic;			/* bic. */
	double		aic;			/* aic. */
	double		icl;			/* icl. */
	int		*class_id;		/* Class id of each sequence, dim = N_X. */
	int		*n_class;		/* Number of sequences for each class, dim = K. */

/* Extension applications. */
	/* Label. */
	phyclust_label	*label;			/* For labels. */

	/* For sequencing error if assumed. */
	int		se_type;		/* For sequencing error. */
	SE_P_matrix	*SE_P;			/* Parameters of sequencing error models. */
};

phyclust_struct* initialize_phyclust_struct(int code_type, int N_X_org, int L, int K);
void free_phyclust_struct(phyclust_struct *pcs);
void update_phyclust_struct(phyclust_struct *pcs);


/* ----- Summary tool. ----- */
void assign_class(phyclust_struct *pcs);
void update_ic(phyclust_struct *pcs, Q_matrix_array *QA);


/* ----- For debug. ----- */
void print_X(phyclust_struct *pcs);
void print_Mu(phyclust_struct *pcs);
void print_class_id(phyclust_struct *pcs);

#endif	/* End of __PHYCLUST_STRUCT_. */

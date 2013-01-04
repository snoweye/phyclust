/* This file contains declarations for phyclust. */


#ifndef __PHYCLUST_INIT_METHOD_
#define __PHYCLUST_INIT_METHOD_

#include "phyclust_em.h"
#include "phyclust_ape_nj.h"


/* Tools for initialization. */
int rdunif(int n);
void srswor(int n, int k, int *x);
int init_m_step(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
int check_all_min_n_class(int k, int *n_class, int min_n_class);
void assign_Mu_by_class(int N_X_org, int K, int L, int ncode, int gap_index, int *class_id, int **X_org, int **Mu);

/* consensus_Mu[L]. */
void find_consensus_Mu(int N_X_org, int L, int ncode, int gap_index, int **X_org, int *consensus_Mu);
void find_consensus_Mu_gap(int N_X_org, int L, int ncode, int gap_index, int **X_org, int *consensus_Mu);


/* Initialization functions.
 * All functions will update empcs, QA.
 * Basically, these functions only provides Mu and update QA.
 * Further updates for em are implemented "phyclust_init_procedure.c".
 * All functions will return 1 for fail, 0 for success. */

/* These functions with "_unique" pick centers from unique sequences
 * unlike other methods, then maps id back to usual "X_org". */

/* Randomly pick Mu's and assign by edist. */
int Update_init_random_Mu_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
int Update_init_random_Mu_unique_label(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By neighbor-joining. */
int Update_init_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By random neighbor-joining. */
int Update_init_random_nj_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By pam. */
int Update_init_pam(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By k-medoids. */
int Update_init_k_medoids(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By hand-coding or a prespecified input file. */
int Update_init_manually(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);

/* By sampling each locus. */
int Update_init_sampleL_unique(em_phyclust_struct *empcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);


/* File: "phyclust_init_method_nj.c". */
void search_largest_branch(nj_struct *njs, int *largest_branch_id);
void random_branch(nj_struct *njs, int *random_branch_id);
int assign_class_by_njs_branch(int K, nj_struct *njs, int *branch_id, int *class_id);

/* File: "phyclust_init_method_kmed.c". */
void assign_class_by_k_medoids(int N_X, int K, double **EDM, int *center_id, int *class_id);
void assign_class_unique_by_k_medoids(int N_X_org, int K, double **EDM, int N_X, int *map_X_to_X_org,
		int *center_id, int *class_id);

/* File: "phyclust_init_method_pam.c". */
void assign_class_by_pam(int N_X, int K, double **EDM_LT_pam, int *center_id, int *class_id);


/* ----- For debug. ----- */
void print_consensus_Mu(em_phyclust_struct *empcs, int *consensus_Mu);

#endif	/* End of __PHYCLUST_INIT_METHOD_. */

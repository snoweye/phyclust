/* This file contains declarations for em steps. */

#ifndef __PHYCLUST_EM_TOOL_
#define __PHYCLUST_EM_TOOL_

#include "phyclust_struct.h"
#include "phyclust_em.h"
#include "phyclust_qmatrix_array.h"


/* ----- Independent summary tool. ----- */
double LogL_observed(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_observed_label_semi(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_observed_label_general(em_phyclust_struct *empcs, Q_matrix_array *QA);
/* For debug only. */
double LogL_complete(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);		/* QA_H unused */
double LogL_complete_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);	/* QA_H != QA */
/* For M-step. */
double LogL_profile(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);		/* QA_H unused */
double LogL_profile_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);	/* QA_H != QA */

/* ----- Initialization tool. ----- */
void initialize_count_Mu_X_and_gap(em_phyclust_struct *empcs);
void reset_Mu_non_seg_site(em_phyclust_struct *empcs);

/* ----- Checking tool. ----- */
int is_finite(double x);

/* ----- For copy. ----- */
void copy_EMC(em_control *EMC_from, em_control *EMC_to);
void reassign_label_pointer(em_phyclust_struct *empcs);
void Copy_empcs(em_phyclust_struct *empcs_from, em_phyclust_struct *empcs_to);
void Copy_empcs_to_pcs(em_phyclust_struct *empcs, phyclust_struct *pcs);
/* For M-step lonely, or semi-supervised. */
void Copy_pcs_to_empcs(phyclust_struct *pcs, em_phyclust_struct *empcs);
void Copy_pcs_to_empcs_label(phyclust_struct *pcs, em_phyclust_struct *empcs);


/* ----- For debug. ----- */
void print_empcs(em_phyclust_struct *empcs);
void print_EMC(em_control *EMC);
void print_rich_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC);
void print_result(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC);
void print_Z_modified(em_phyclust_struct *empcs);
void print_Z_normalized(em_phyclust_struct *empcs);
void print_Eta(em_phyclust_struct *empcs);
void print_empcs_Mu(em_phyclust_struct *empcs);
void print_empcs_Mu_seg_site(em_phyclust_struct *empcs);
void print_count_Mu_X(em_phyclust_struct *empcs, int n_X, int k);
void print_count_Mu_X_gap(em_phyclust_struct *empcs, int n_X, int k);

#endif	/* End of __PHYCLUST_EM_TOOL_. */

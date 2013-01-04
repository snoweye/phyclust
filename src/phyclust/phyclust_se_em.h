/* For sequencing error models. */

/* This file contains declarations for EM with sequencing error models. */

#ifndef __PHYCLUST_SE_EM_
#define __PHYCLUST_SE_EM_

#include "phyclust_struct.h"
#include "phyclust_em.h"
#include "phyclust_qmatrix_array.h"


/* ----- Initial functions in "phyclust_se_em_phyclust_struct.c". ----- */
void initialize_em_phyclust_struct_se(em_phyclust_struct *empcs, phyclust_struct *pcs);
void free_em_phyclust_struct_se(em_phyclust_struct *empcs);
void duplicate_em_phyclust_struct_se(em_phyclust_struct *org_empcs, em_phyclust_struct *new_empcs);

/* ----- Em functions in "phyclust_se_em_step.c". ----- */
/* This special function updates empcs->Z_modified (log and unnormalized). */
void Update_Z_modified_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA);
void Update_Z_modified_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA);

/* ----- Control function pointers in "phyclust_se_em_fp.c". ----- */
void update_em_fp_se(em_fp *EMFP, em_control *EMC, phyclust_struct *pcs);


/* ----- Tool functions in "phyclust_se_convolution_em_tool.c". ----- */
/* Independent summary tool. */
double LogL_observed_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA);
double LogL_observed_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA);

/* For debug only. */
double LogL_complete_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);		/* QA_H == QA */
double LogL_complete_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);	/* QA_H == QA */

/* For M-step. */
double LogL_profile_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);		/* QA_H == QA */
double LogL_profile_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);	/* QA_H == QA */

/* For copy.*/
void Copy_empcs_se_convolution(em_phyclust_struct *empcs_from, em_phyclust_struct *empcs_to);

void Copy_empcs_to_pcs_se(em_phyclust_struct *empcs, phyclust_struct *pcs);
void Copy_pcs_to_empcs_se(phyclust_struct *pcs, em_phyclust_struct *empcs);

/* Utility function. */
void update_convolution_Pt_f_err(Q_matrix_array *QA, SE_P_matrix *SE_P);
void update_convolution_Pt_f_err_gap(Q_matrix_array *QA, SE_P_matrix *SE_P);
void print_convolution_Pt_f_err(double ***log_conv, int K, int nrow, int ncol);

#endif	/* End of __PHYCLUST_SE_EM_. */


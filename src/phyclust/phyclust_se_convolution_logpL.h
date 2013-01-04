/* This file contains declarations for em steps. */

#ifndef __PHYCLUST_SE_CONVOLUTION_LOGPL_
#define __PHYCLUST_SE_CONVOLUTION_LOGPL_

#include "phyclust_logpL.h"

double negative_logpL_Mu_given_QA_se_convolution(int m, double *vect, void *ex);
int Maximize_logpL_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP);

/* For update Mu given QA. */
/* QA_H unused */
void Update_Mu_given_QA_full_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
void Update_Mu_given_QA_skip_non_seg_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
/* QA_H unused */
void Update_Mu_given_QA_full_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
void Update_Mu_given_QA_skip_non_seg_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);

#endif	/* End of __PHYCLUST_SE_CONVOLUTION_LOGPL_. */

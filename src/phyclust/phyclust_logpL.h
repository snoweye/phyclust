/* This file contains declarations for em steps. */

#ifndef __PHYCLUST_LOGPL_
#define __PHYCLUST_LOGPL_

#include "phyclust_em.h"
#include "phyclust_struct.h"
#include "phyclust_qmatrix_array.h"

/* Internal used only, to convert between struct and vector for optim functions. */
typedef struct _ex_struct		ex_struct;


/* For maximizing subroutine. */
struct _ex_struct{
	em_phyclust_struct	*empcs;
	em_fp			*EMFP;
	Q_matrix_array		*QA;
	Q_matrix_array		*QA_H;
	double			*org_vect;
};

double negative_logpL_Mu_given_QA(int m, double *vect, void *ex);
double negative_logpL_QA_given_Mu(int m, double *vect, void *ex);
int Maximize_logpL(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP);

/* For update Mu given QA. */
/* QA_H unused */
void Update_Mu_given_QA_full(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
void Update_Mu_given_QA_skip_non_seg(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
/* QA_H != QA */
void Update_Mu_given_QA_full_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
void Update_Mu_given_QA_skip_non_seg_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);


/* For update the function R(). */
/* QA_H unused */
double Compute_R(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);
/* QA_H != QA */
double Compute_R_gap(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H);


/* ----- For debug. ----- */
void print_vect(int m, double *vect);

#endif	/* End of __PHYCLUST_LOGPL_. */

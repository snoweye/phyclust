#include <R.h>
#include <Rinternals.h>
#include "phyclust/phyclust.h"
#include "R_phyclust.h"

/* In file "R_phyclust_se.c". */
void copy_R_EMC_to_EMC_se(SEXP R_EMC, em_control *EMC);

typedef struct _emptr_se	*EMPTR_SE;

struct _emptr_se{
	/* int C_protect_length; */
	int *C_N_X_org, *C_N_X, *C_L, *C_K, *C_p,
	    *C_converge_flag, *C_converge_iter, *C_converge_inner_iter, *C_converge_cm_iter, *C_check_param,
	    *C_N_seg_site, *C_label_method;
	int *C_se_type, *C_se_model;
	double *C_logL, *C_bic, *C_aic, *C_icl, *C_pi, *C_kappa, *C_Tt,
	       *C_converge_eps, *C_converge_error;
	double *C_se_constant, *C_se_f_err;
};

EMPTR_SE allocate_emptr_se(void);
SEXP initialize_emptr_se(EMPTR_SE emptr_se, phyclust_struct *pcs);
void update_emptr_se(EMPTR_SE emptr_se, phyclust_struct *pcs, SEXP emobj);
void copy_all_to_emptr_se(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, EMPTR_SE emptr_se);


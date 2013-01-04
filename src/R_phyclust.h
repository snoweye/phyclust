#include <R.h>
#include <Rinternals.h>
#include "phyclust/phyclust.h"

/* In file "R_phyclust.c". */
SEXP getListElement(SEXP list, const char *str);
void copy_R_EMC_to_EMC(SEXP R_EMC, em_control *EMC);

typedef struct _emptr	*EMPTR;

struct _emptr{
	int C_protect_length;
	int *C_N_X_org, *C_N_X, *C_L, *C_K, *C_p,
	    *C_converge_flag, *C_converge_iter, *C_converge_inner_iter, *C_converge_cm_iter, *C_check_param,
	    *C_N_seg_site, *C_label_method;
	double *C_logL, *C_bic, *C_aic, *C_icl, *C_pi, *C_kappa, *C_Tt,
	       *C_converge_eps, *C_converge_error;
};

EMPTR allocate_emptr(void);
SEXP initialize_emptr(EMPTR emptr, phyclust_struct *pcs);
void copy_all_to_emptr(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, EMPTR emptr);


/* In file "./R_phyclust_struct.c". */
phyclust_struct* R_initialize_phyclust_struct(int code_type, int N_X_org, int L, int K);
void R_free_phyclust_struct(phyclust_struct *pcs);

/* In file "R_phyclust_label.c". */
void R_update_phyclust_label(phyclust_struct *pcs, SEXP R_label);


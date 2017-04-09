#include <R.h>
#include <R_ext/Rdynload.h>

#include "zzz.h"

static const R_CallMethodDef callMethods[] = {
	{"R_phyclust", (DL_FUNC) &R_phyclust, 7},
	{"R_phyclust_edist", (DL_FUNC) &R_phyclust_edist, 4},
	{"R_phyclust_em_step", (DL_FUNC) &R_phyclust_em_step, 11},
	{"R_phyclust_e_step", (DL_FUNC) &R_phyclust_e_step, 12},
	{"R_phyclust_m_step", (DL_FUNC) &R_phyclust_m_step, 10},
	{"R_phyclust_find_consensus", (DL_FUNC) &R_phyclust_find_consensus, 5},
	{"R_phyclust_logL", (DL_FUNC) &R_phyclust_logL, 11},
	{"R_phyclust_logPt", (DL_FUNC) &R_phyclust_logPt, 5},
	{"R_phyclust_se", (DL_FUNC) &R_phyclust_se, 7},
	{"R_phyclust_se_update", (DL_FUNC) &R_phyclust_se_update, 9},
	{"R_phyclust_update", (DL_FUNC) &R_phyclust_update, 9},
	{"R_RRand", (DL_FUNC) &R_RRand, 5},

	{"R_ms_main", (DL_FUNC) &R_ms_main, 2},
	{"R_paml_baseml_main", (DL_FUNC) &R_paml_baseml_main, 2},
	{"R_seq_gen_main", (DL_FUNC) &R_seq_gen_main, 2},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_phyclust(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_phyclust(). */

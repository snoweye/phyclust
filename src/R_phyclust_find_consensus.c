/* This file contains a function R_phyclust_find_consensus() called by R wraps
 * find_consensus() in "R/f_find_consensus.r", and this function calls the
 * relative functions
 * find_consensus_Mu() in "src/phyclust_init_method.c".
 *
 * Writen: Wei-Chen Chen on 2010/02/26. */

#include <R.h>
#include <Rinternals.h>
#include "phyclust/phyclust.h"

/* This function calls find_consensus_Mu() in "src/phyclust_init_method.c" and
 * is called by find.consensus() using .Call() in "R/f_find_consensus.r".
 * Input:
 *   R_N_X_org: SEXP[1], number of sequences.
 *   R_L: SEXP[1], length of sequences.
 *   R_code_type: SEXP[1], code type.
 *   R_WIGAP: SEXP[1], with gap.
 *   R_X_org: SEXP[1], sequences.
 * Output:
 *   ret: SEXP[R_L], an array contains Mu. */
SEXP R_phyclust_find_consensus(SEXP R_N_X_org, SEXP R_L, SEXP R_code_type,
		SEXP R_WIGAP, SEXP R_X_org){
	/* Declare variables for calling C. */
	int *C_N_X_org, *C_L, *C_code_type, *C_WIGAP, **C_X_org;

	/* Declare variables for R's returning. */
	SEXP ret;
	int *ret_ptr;

	/* Declare variables for processing. */
	int i, *tmp_ptr;
  
	/* Set initial values. */
	C_N_X_org = INTEGER(R_N_X_org);
	C_L = INTEGER(R_L);
	C_code_type = INTEGER(R_code_type);
	C_WIGAP = INTEGER(R_WIGAP);

	/* Assign data. */
	C_X_org = allocate_int_2D_AP(*C_N_X_org);
	tmp_ptr = INTEGER(R_X_org);
	for(i = 0; i < *C_N_X_org; i++){
		C_X_org[i] = tmp_ptr;
		tmp_ptr += *C_L;
	}

	/* For return. */
	PROTECT(ret = allocVector(INTSXP, *C_L));
	ret_ptr = INTEGER(ret);

	/* Compute. */
	if(*C_WIGAP == 0){
        	find_consensus_Mu(*C_N_X_org, *C_L,
			NCODE[*C_code_type],
			GAP_INDEX[*C_code_type], C_X_org, ret_ptr);
	} else{
        	find_consensus_Mu_gap(*C_N_X_org, *C_L,
			NCODE_WIGAP[*C_code_type],
			GAP_INDEX[*C_code_type], C_X_org, ret_ptr);
	}

	/* Free memory and release protectation. */
	UNPROTECT(1);

	return(ret);
} /* End of SEXP R_phyclust_find_consensus(). */


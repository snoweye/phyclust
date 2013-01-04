/* For sequencing error models. */

/* This file contains all functions required in em steps.
 *
 * WARNING: This is not a well optimized version and not quite extensible.
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_struct.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"
#include "phyclust_se_em.h"
#include "phyclust_se_pmatrix.h"


/* For sequencing error models. */
void initialize_em_phyclust_struct_se(em_phyclust_struct *empcs, phyclust_struct *pcs){
	if(empcs->code_type == NUCLEOTIDE && empcs->se_type == SE_YES){
		empcs->SE_P = duplicate_SE_P_matrix(pcs->SE_P);
	}
} /* End of initialize_em_phyclust_struct_se(). */

void free_em_phyclust_struct_se(em_phyclust_struct *empcs){
	if(empcs->code_type == NUCLEOTIDE && empcs->se_type == SE_YES){
		free_SE_P_matrix(empcs->SE_P);
	}
} /* End of free_em_phyclust_struct_se(). */

void duplicate_em_phyclust_struct_se(em_phyclust_struct *org_empcs, em_phyclust_struct *new_empcs){
	if(org_empcs->code_type == NUCLEOTIDE && org_empcs->se_type == SE_YES){
		new_empcs->SE_P = duplicate_SE_P_matrix(org_empcs->SE_P);
	}
} /* End of duplicate_em_phyclust_struct_se(). */


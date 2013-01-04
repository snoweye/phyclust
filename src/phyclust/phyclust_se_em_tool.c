/* For sequencing error models. */

/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"
#include "phyclust_se_em.h"
#include "phyclust_se_pmatrix.h"


/* ----- For copy. ----- */
void Copy_empcs_se_convolution(em_phyclust_struct *empcs_from, em_phyclust_struct *empcs_to){
	Copy_empcs(empcs_from, empcs_to);
	copy_SE_P_matrix(empcs_from->SE_P, empcs_to->SE_P);
} /* End of Copy_empcs_se_convolution(). */


void Copy_empcs_to_pcs_se(em_phyclust_struct *empcs, phyclust_struct *pcs){
	Copy_empcs_to_pcs(empcs, pcs);
	copy_SE_P_matrix(empcs->SE_P, pcs->SE_P);
} /* End of Copy_empcs_to_pcs_se(). */


void Copy_pcs_to_empcs_se(phyclust_struct *pcs, em_phyclust_struct *empcs){
	Copy_pcs_to_empcs(pcs, empcs);
	copy_SE_P_matrix(pcs->SE_P, empcs->SE_P);
} /* End of Copy_pcs_to_empcs_se(). */


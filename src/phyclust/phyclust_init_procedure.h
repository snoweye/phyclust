/* This file contains declarations for phyclust. */


#ifndef __PHYCLUST_INIT_PROCEDURE_
#define __PHYCLUST_INIT_PROCEDURE_

#include "phyclust_em.h"

void init_em_step(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);
/* For models using empirical pi's in Q. */
void update_Q_matrix_array(Q_matrix_array *QA, phyclust_struct *pcs);


/* Initialization procedures with the prespecified initialization method
 * indicated in init_em_step(). */

/* Exhaust all the initializations with EM. */
void exhaust_EM(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);	/* exhaustEM. */

/* Run EM on the best random initialization. */
void Rnd_EM(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);	/* RndEM. */

/* Run EM on the best random initialization with fiexed iterations of EM updates. */
void Rndp_EM(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);	/* RndpEM. */

/* Run EM on the best short em initialization. */
void em_EM(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP);		/* emEM. */

#endif	/* End of __PHYCLUST_INIT_PROCEDURE_. */

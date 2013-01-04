/* For sequencing error models. */

/* This file contains declarations for phyclust. */

#ifndef __PHYCLUST_SE_STRUCT_
#define __PHYCLUST_SE_STRUCT_

#include "phyclust_struct.h"
#include "phyclust_em.h"

void free_phyclust_se_struct(phyclust_struct *pcs);
void update_phyclust_se_struct(phyclust_struct *pcs, em_control *EMC);

#endif	/* End of __PHYCLUST_SE_STRUCT_. */

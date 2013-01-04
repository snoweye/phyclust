/* This file contains all header files. */

#ifndef __PHYCLUST_
#define __PHYCLUST_

/* Define constants. */
#include "phyclust_constant.h"

/* Define phyclust structures. */
#include "phyclust_struct.h"
#include "phyclust_label.h"
#include "phyclust_tool.h"

/* Define Q matrix and distance. */
#include "phyclust_qmatrix.h"
#include "phyclust_edist.h"
#include "phyclust_qmatrix_array.h"

/* Define EM related structures and functions. */
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"

/* Define initialization procedures and methods. */
#include "phyclust_init_procedure.h"
#include "phyclust_init_method.h"
#include "phyclust_ape_nj.h"

/* Define NM and BFGS. */
#include "phyclust_optim_nmmin.h"

/* Define input file format: phylip, fasta. */
#include "phyclust_file_input.h"

/* Define sequencing error models. */
#include "phyclust_se_struct.h"
#include "phyclust_se_pmatrix.h"
#include "phyclust_se_em.h"
#include "phyclust_se_convolution_logpL.h"

#endif	/* End of __PHYCLUST_. */

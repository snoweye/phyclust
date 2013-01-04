/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"


/* Initial a em controler. */
em_control* initialize_em_control(){
	em_control *EMC;

	EMC = (em_control*) malloc(sizeof(em_control));
	EMC->exhaust_iter = 1;
	EMC->fixed_iter = 5;
	EMC->short_iter = 100;
	EMC->EM_iter = 1000;
	EMC->short_eps = 1e-2;
	EMC->EM_eps = 1e-6;

	EMC->cm_reltol = 1e-8;
	EMC->cm_maxit = 5000;

	EMC->nm_abstol_Mu_given_QA = 1e-8;	/* for update_flag = 0. */
	EMC->nm_abstol_QA_given_Mu = 1e-8;	/* for update_flag = 1. */
	EMC->nm_reltol_Mu_given_QA = 1e-8;
	EMC->nm_reltol_QA_given_Mu = 1e-8;
	EMC->nm_maxit_Mu_given_QA = 500;
	EMC->nm_maxit_QA_given_Mu = 5000;
	EMC->est_non_seg_site = 0;		/* 1 for original case, 0 for skip. */

	EMC->max_init_iter = 50;
	EMC->min_n_class = 1;
	EMC->init_procedure = exhaustEM;
	EMC->init_method = randomMu;
	EMC->substitution_model = HKY85;
	EMC->edist_model = D_HAMMING;
	EMC->identifier = EE;
	EMC->code_type = NUCLEOTIDE;
	EMC->em_method = EM;
	EMC->boundary_method = ADJUST;

	EMC->Eta_lower_bound = 1e-8;		/* The default value should be 1 / N_X_org. */
	EMC->Eta_upper_bound = 1.0;		/* The default value should be 1 - 1 / N_X_org. */

	EMC->converge_eps = 0.0;
	EMC->converge_error = 0.0;
	EMC->converge_flag = 0;
	EMC->converge_iter = 0;
	EMC->converge_inner_iter = 0;
	EMC->converge_cm_iter = 0;
	EMC->update_flag = 0;		/* 0 to apply optimization on QA, Tt and update Mu given QA and Tt;
					   1 to apply optimization on QA and Tt given Mu. */

	EMC->se_type = SE_NO;
	EMC->se_model = SE_CONVOLUTION;
	EMC->se_constant = SE_CONSTANT;

	update_em_control(EMC);
	return(EMC);
} /* End of initialize_em_control(). */

void free_em_control(em_control *EMC){
	free(EMC);
} /* End of free_em_control(). */

em_control* duplicate_em_control(em_control *org_EMC){
	em_control *new_EMC;

	new_EMC = initialize_em_control();
	new_EMC->exhaust_iter = org_EMC->exhaust_iter;
	new_EMC->fixed_iter = org_EMC->fixed_iter;
	new_EMC->short_iter = org_EMC->short_iter;
	new_EMC->EM_iter = org_EMC->EM_iter;
	new_EMC->short_eps = org_EMC->short_eps;
	new_EMC->EM_eps = org_EMC->EM_eps;

	new_EMC->cm_reltol = org_EMC->cm_reltol;
	new_EMC->cm_maxit = org_EMC->cm_maxit;

	new_EMC->nm_abstol_Mu_given_QA = org_EMC->nm_abstol_Mu_given_QA;
	new_EMC->nm_abstol_QA_given_Mu = org_EMC->nm_abstol_QA_given_Mu;
	new_EMC->nm_reltol_Mu_given_QA = org_EMC->nm_reltol_Mu_given_QA;
	new_EMC->nm_reltol_QA_given_Mu = org_EMC->nm_reltol_QA_given_Mu;
	new_EMC->nm_maxit_Mu_given_QA = org_EMC->nm_maxit_Mu_given_QA;
	new_EMC->nm_maxit_QA_given_Mu = org_EMC->nm_maxit_QA_given_Mu;
	new_EMC->est_non_seg_site = org_EMC->est_non_seg_site;

	new_EMC->max_init_iter = org_EMC->max_init_iter;
	new_EMC->min_n_class = org_EMC->min_n_class;
	new_EMC->init_procedure = org_EMC->init_procedure;
	new_EMC->init_method = org_EMC->init_method;
	new_EMC->substitution_model = org_EMC->substitution_model;
	new_EMC->edist_model = org_EMC->edist_model;
	new_EMC->identifier = org_EMC->identifier;
	new_EMC->code_type = org_EMC->code_type;
	new_EMC->em_method = org_EMC->em_method;
	new_EMC->boundary_method = org_EMC->boundary_method;

	new_EMC->Eta_lower_bound = org_EMC->Eta_lower_bound;
	new_EMC->Eta_upper_bound = org_EMC->Eta_upper_bound;

	new_EMC->se_type = org_EMC->se_type;
	new_EMC->se_model = org_EMC->se_model;
	new_EMC->se_constant = org_EMC->se_constant;

	copy_EMC(org_EMC, new_EMC);
	return(new_EMC);
} /* End of duplicate_em_control(). */

void update_em_control(em_control *EMC){
	/* Check code_type. */
	switch(EMC->code_type){
		case SNP:
			EMC->edist_model = D_HAMMING;
			if(EMC->substitution_model != SNP_JC69 ||
			   EMC->substitution_model != SNP_F81 ||
			   EMC->substitution_model != E_SNP_F81){
				EMC->substitution_model = SNP_JC69;
			}
			break;
		case NUCLEOTIDE:
			if(EMC->substitution_model == SNP_JC69 ||
			   EMC->substitution_model == SNP_F81 ||
			   EMC->substitution_model == E_SNP_F81){
				EMC->substitution_model = JC69;
			}
			break;
		default:
			fprintf_stderr("PE: The code type is not found.\n");
			exit(1);
	}

	/* Change settings. */
	if(EMC->init_method == NJ || EMC->init_method == PAM ||
	   EMC->init_method == manualMu){
		EMC->exhaust_iter = 1;
		EMC->init_procedure = exhaustEM;
	}
	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
} /* End of update_em_control(). */

void reset_em_control(em_control *EMC){
	EMC->converge_eps = 0.0;
	EMC->converge_error = 0.0;
	EMC->converge_flag = 0;
	EMC->converge_iter = 0;
	EMC->converge_inner_iter = 0;
	EMC->converge_cm_iter = 0;
	EMC->update_flag = (EMC->em_method == EM) ? 0 : 1;
} /* End of reset_em_control(). */


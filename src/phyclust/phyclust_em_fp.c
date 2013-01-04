/* This file contains all function pointers required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_init_method.h"
#include "phyclust_logpL.h"
#include "phyclust_se_em.h"

/* Initial a em_fp. */
em_fp* initialize_em_fp(em_control *EMC, phyclust_struct *pcs){
	em_fp *EMFP;
	
	EMFP = (em_fp*) malloc(sizeof(em_fp));
	
	/* Assign init_method FP. */
	switch(EMC->init_method){
		case randomMu:
			EMFP->Update_init = &Update_init_random_Mu_unique;
			break;
		case NJ:
			EMFP->Update_init = &Update_init_nj_unique;
			break;
		case randomNJ:
			EMFP->Update_init = &Update_init_random_nj_unique;
			break;
		case PAM:
			EMFP->Update_init = &Update_init_pam;
			break;
		case kMedoids:
			EMFP->Update_init = &Update_init_k_medoids;
			break;
		case manualMu:
			EMFP->Update_init = &Update_init_manually;
			break;
		default:
			fprintf_stderr("PE: The initial method is not found.\n");
			exit(1);
	}

	/* Assign EM method. */
	switch(EMC->em_method){
		case EM:
			EMFP->M_step = &M_step_simple;
			EMFP->Check_convergence = &Check_convergence_em;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			EMFP->Update_Z_modified = &Update_Z_modified;
			EMFP->Maximize_logpL = &Maximize_logpL;
			break;
		case ECM:
			EMFP->M_step = &M_step_CM;
			EMFP->Check_convergence = &Check_convergence_org;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			EMFP->Update_Z_modified = &Update_Z_modified;
			EMFP->Maximize_logpL = &Maximize_logpL;
			break;
		case AECM:
			EMFP->M_step = &M_step_ACM;
			EMFP->Check_convergence = &Check_convergence_org;
			EMFP->Em_step = &Em_step;
			EMFP->Short_em_step = &Short_em_step;
			EMFP->Update_Z_modified = &Update_Z_modified;
			EMFP->Maximize_logpL = &Maximize_logpL;
			break;
		default:
			fprintf_stderr("PE: The EM method is not found.\n");
			exit(1);
	}

	/* Update functions. */
	switch(pcs->label->label_method){
		case NONE:
			EMFP->E_step_logL_observed = &E_step_logL_observed;
			EMFP->LogL_observed = &LogL_observed;
			EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs;
			break;
		case SEMI:
			EMFP->Update_init = &Update_init_random_Mu_unique_label;
			EMFP->E_step_logL_observed = &E_step_logL_observed_label_semi;
			EMFP->LogL_observed = &LogL_observed_label_semi;
			EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs_label;
			break;
		case GENERAL:
			EMFP->Update_init = &Update_init_random_Mu_unique_label;
			EMFP->E_step_logL_observed = &E_step_logL_observed_label_general;
			EMFP->LogL_observed = &LogL_observed_label_general;
			EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs_label;
			break;
		default:
			fprintf_stderr("PE: The label method is not found.\n");
			exit(1);
	}

	/* For boundary methods. */
	switch(EMC->boundary_method){
		case IGNORE:
			EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_IGNORE;
			break;
		case ADJUST:
			EMFP->Update_Eta_given_Z = &Update_Eta_given_Z_ADJUST;
			break;
		default:
			fprintf_stderr("PE: The boundary method is not found.\n");
			exit(1);
	}

	/* For GAPs. */
	if(pcs->gap_flag){
		EMFP->LogL_complete = &LogL_complete_gap;
		EMFP->LogL_profile = &LogL_profile_gap;
		EMFP->Compute_R = &Compute_R_gap;
		if(EMC->est_non_seg_site != 0){
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full_gap;
		} else{
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg_gap;
		}
	} else{
		EMFP->LogL_complete = &LogL_complete;
		EMFP->LogL_profile = &LogL_profile;
		EMFP->Compute_R = &Compute_R;
		if(EMC->est_non_seg_site != 0){
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full;
		} else{
			EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg;
		}
	}

	EMFP->Copy_empcs = &Copy_empcs;
	EMFP->Copy_empcs_to_pcs = &Copy_empcs_to_pcs;

	/* For sequencing error model. */
	update_em_fp_se(EMFP, EMC, pcs);
	
	return(EMFP);
} /* End of initialize_em_fp(). */


void free_em_fp(em_fp *EMFP){
	free(EMFP);
} /* End of free_em_fp(). */


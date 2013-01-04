/* For sequencing error models. */

/* This file contains all function pointers required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"
#include "phyclust_se_em.h"
#include "phyclust_se_convolution_logpL.h"

/* Initial a em_fp. */
void update_em_fp_se(em_fp *EMFP, em_control *EMC, phyclust_struct *pcs){
	if(pcs->se_type == SE_YES){
		if(EMC->em_method != EM){
			fprintf_stderr("PE: The em_method is not implemented.\n");
			exit(1);
		}
		if(pcs->label->label_method != NONE){
			fprintf_stderr("PE: The semi-supervised method with sqeuencing error is not implemented.\n");
			exit(1);
		}

		/* Functions need to be implemented for sequencing error model. */
		switch(EMC->se_model){
			case SE_CONVOLUTION:
				/* Overwrite any necessary function points in "phyclust_em_fp.c" since
				 * the sequencing error model may have different implementations.
				 * In the order of "struct _em_fp" in "phyclust_em.h".
				 * These may be assigned by the function "initialize_em_fp()", but
				 * I redundantly overwrite again in case of miss specify the function pointer. */

				/* The same as EM, no need for reimplementation if inner function pointers
				 * are implemented and specified appropriately. */
				EMFP->M_step = &M_step_simple;
				EMFP->Check_convergence = &Check_convergence_em;
				EMFP->Em_step = &Em_step;
				EMFP->Short_em_step = &Short_em_step;
				EMFP->E_step_logL_observed = &E_step_logL_observed;
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

				/* Different implementations for EM.
				 * The gap implementations are totally different to the original
				 * gap mechanism. Original, we can ignore GAPs in observed logL, but
				 * not complete logL. Here, we have to take into account the error probability
				 * int observed logL, and adjust for complete logL, and maximize profile logL.
				 * The memory copy mechanism is also needed to take care. */
				if(pcs->gap_flag){
					EMFP->Update_Z_modified = &Update_Z_modified_gap_se_convolution;
				} else{
					EMFP->Update_Z_modified = &Update_Z_modified_se_convolution;
				}

				EMFP->Maximize_logpL = &Maximize_logpL_se_convolution;

				if(pcs->gap_flag){
					EMFP->LogL_observed = &LogL_observed_gap_se_convolution;
					EMFP->LogL_complete = &LogL_complete_gap_se_convolution;
					EMFP->LogL_profile = &LogL_profile_gap_se_convolution;
				} else{
					EMFP->LogL_observed = &LogL_observed_se_convolution;
					EMFP->LogL_complete = &LogL_complete_se_convolution;
					EMFP->LogL_profile = &LogL_profile_se_convolution;
				}
				EMFP->Copy_empcs = &Copy_empcs_se_convolution;
				EMFP->Copy_pcs_to_empcs = &Copy_pcs_to_empcs_se;
				EMFP->Copy_empcs_to_pcs = &Copy_empcs_to_pcs_se;

				if(pcs->gap_flag){
					if(EMC->est_non_seg_site != 0){
						EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full_gap_se_convolution;
					} else{
						EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg_gap_se_convolution;
					}
				} else{
					if(EMC->est_non_seg_site != 0){
						EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_full_se_convolution;
					} else{
						EMFP->Update_Mu_given_QA = &Update_Mu_given_QA_skip_non_seg_se_convolution;
					}
				}

				/* The same as EM without gap. */
				EMFP->Compute_R = &Compute_R;
				break;
			default:
				fprintf_stderr("PE: The SE_P model is not found.\n");
				exit(1);
				break;
		}
	}
} /* End of update_em_fp_se(). */


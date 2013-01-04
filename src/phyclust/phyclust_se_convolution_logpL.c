/* This file contains functions about log profile likelihood. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_em.h"
#include "phyclust_em_tool.h"
#include "phyclust_logpL.h"
#include "phyclust_optim_nmmin.h"
#include "phyclust_tool.h"
#include "phyclust_se_em.h"
#include "phyclust_se_convolution_logpL.h"


/* This function will update Mu given QA for full length sequences. */
void Update_Mu_given_QA_full_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){
	int org_s_mu, new_s_mu, n_X, k, l, flag = 0;
	int K = empcs->K, L = empcs->L, Ncode = empcs->ncode, N_X = empcs->N_X;
	double tmp_psi_k, psi_k[empcs->ncode], z_k[empcs->ncode];

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				z_k[org_s_mu] = 0.0;
			}

			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->X[n_X][l] == MISSING_ALLELE){
					continue;
				}
				if(empcs->replication_X[n_X] == 1){
					z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
				} else{
					z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
				}
			}

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				psi_k[org_s_mu] = 0.0;
				for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
					psi_k[org_s_mu] += z_k[new_s_mu] * empcs->SE_P->log_conv[k][org_s_mu][new_s_mu];
				}
			}

			org_s_mu = empcs->Mu[k][l];
			tmp_psi_k = psi_k[org_s_mu];
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				if(psi_k[new_s_mu] > tmp_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					tmp_psi_k = psi_k[new_s_mu];
					flag |= 1;
				}
			}

			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]-- and
			 *      count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++
			 * for all n_X = 1, ..., N_X.
 			 * "ALWAYS and ONLY" unique sequences (decided by empcs->N_X) are counted.
			 * For non-unique data, empcs->replication_X will involve or multiply to the targets,
			 * such as computing likelihood, finding Mu. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					if(empcs->X[n_X][l] == MISSING_ALLELE){
						continue;
					}
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_se_convolution(). */

/* This function will update Mu given QA, but skip non-segregating sites. */
void Update_Mu_given_QA_skip_non_seg_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){
	int org_s_mu, new_s_mu, n_X, k, l, i, flag = 0;
	int K = empcs->K, N_seg_site = empcs->N_seg_site, Ncode = empcs->ncode, N_X = empcs->N_X;
	double tmp_psi_k, psi_k[empcs->ncode], z_k[empcs->ncode];

	update_convolution_Pt_f_err(QA, empcs->SE_P);

	for(k = 0; k < K; k++){
		for(i = 0; i < N_seg_site; i++){
			l = empcs->seg_site_id[i];

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				z_k[org_s_mu] = 0.0;
			}

			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->X[n_X][l] == MISSING_ALLELE){
					continue;
				}
				if(empcs->replication_X[n_X] == 1){
					z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
				} else{
					z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
				}
			}

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				psi_k[org_s_mu] = 0.0;
				for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
					psi_k[org_s_mu] += z_k[new_s_mu] * empcs->SE_P->log_conv[k][org_s_mu][new_s_mu];
				}
			}

			org_s_mu = empcs->Mu[k][l];
			tmp_psi_k = psi_k[org_s_mu];
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				if(psi_k[new_s_mu] > tmp_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					tmp_psi_k = psi_k[new_s_mu];
					flag |= 1;
				}
			}

			/* Dynamically update empcs->count_Mu_X:
			 * If empcs->Mu[k][l] were changed from org_s_mu to empcs->Mu[k][l]=new_s_mu,
			 * then count_Mu_X[][][][] should be updated by
			 *      count_Mu_X[n][k][org_s_mu][empcs->X[n][l]]-- and
			 *      count_Mu_X[n][k][empcs->Mu[k][l]][empcs->X[n][l]]++
			 * for all n = 1, ..., N_X and where l's should belong to segregating sites.
 			 * "ALWAYS and ONLY" unique sequences (decided by empcs->N_X) are counted.
			 * For non-unique data, empcs->replication_X will involve or multiply to the targets,
			 * such as computing likelihood, finding Mu. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					if(empcs->X[n_X][l] == MISSING_ALLELE){
						continue;
					}
					empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
					empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_skip_non_seg_se_convolution(). */


/* Missing version: This function will update Mu given QA for full length sequences. */
void Update_Mu_given_QA_full_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){
	int org_s_mu, new_s_mu, n_X, k, l, flag = 0;
	int K = empcs->K, L = empcs->L, Ncode = empcs->ncode, N_X = empcs->N_X;
	int fix_s_mu = empcs->SE_P->ncode;	/* Observed gap. */
	double tmp_psi_k, psi_k[empcs->ncode], z_k[empcs->ncode], z_k_gap;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				z_k[org_s_mu] = 0.0;
			}
			z_k_gap = 0.0;

			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->X[n_X][l] == MISSING_ALLELE){
					continue;
				}
				if(empcs->X[n_X][l] != empcs->gap_index){
					if(empcs->replication_X[n_X] == 1){
						z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
					} else{
						z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
					}
				} else{	/* For gap. */
					if(empcs->replication_X[n_X] == 1){
						z_k_gap += empcs->Z_normalized[n_X][k];
					} else{
						z_k_gap += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
					}
				}
			}

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				psi_k[org_s_mu] = 0.0;
				for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
					psi_k[org_s_mu] += z_k[new_s_mu] * empcs->SE_P->log_conv[k][org_s_mu][new_s_mu];
				}
				psi_k[org_s_mu] += z_k_gap * empcs->SE_P->log_conv[k][org_s_mu][fix_s_mu];	/* For gap. */
			}

			org_s_mu = empcs->Mu[k][l];
			tmp_psi_k = psi_k[org_s_mu];
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				if(psi_k[new_s_mu] > tmp_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					tmp_psi_k = psi_k[new_s_mu];
					flag |= 1;
				}
			}

			/* Dynamically update empcs->count_Mu_X: See non-gap version for details. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					if(empcs->X[n_X][l] == MISSING_ALLELE){
						continue;
					}
					if(empcs->X[n_X][l] != empcs->gap_index){
						empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
						empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
					} else{	/* For gap. */
						empcs->count_Mu_X_gap[n_X][k][org_s_mu]--;
						empcs->count_Mu_X_gap[n_X][k][empcs->Mu[k][l]]++;
					}
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_gap_se_convolution(). */


/* Missing version: This function will update Mu given QA, but skip non-segregating sites. */
void Update_Mu_given_QA_skip_non_seg_gap_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H){
	int org_s_mu, new_s_mu, n_X, k, l, i, flag = 0;
	int K = empcs->K, N_seg_site = empcs->N_seg_site, Ncode = empcs->ncode, N_X = empcs->N_X;
	int fix_s_mu = empcs->SE_P->ncode;	/* Observed gap. */
	double tmp_psi_k, psi_k[empcs->ncode], z_k[empcs->ncode], z_k_gap;

	update_convolution_Pt_f_err_gap(QA, empcs->SE_P);

	for(k = 0; k < K; k++){
		for(i = 0; i < N_seg_site; i++){
			l = empcs->seg_site_id[i];

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				z_k[org_s_mu] = 0.0;
			}
			z_k_gap = 0.0;

			for(n_X = 0; n_X < N_X; n_X++){
				if(empcs->X[n_X][l] == MISSING_ALLELE){
					continue;
				}
				if(empcs->X[n_X][l] != empcs->gap_index){
					if(empcs->replication_X[n_X] == 1){
						z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k];
					} else{
						z_k[empcs->X[n_X][l]] += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
					}
				} else{	/* For gap. */
					if(empcs->replication_X[n_X] == 1){
						z_k_gap += empcs->Z_normalized[n_X][k];
					} else{
						z_k_gap += empcs->Z_normalized[n_X][k] * empcs->replication_X[n_X];
					}
				}
			}

			for(org_s_mu = 0; org_s_mu < Ncode; org_s_mu++){
				psi_k[org_s_mu] = 0.0;
				for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
					psi_k[org_s_mu] += z_k[new_s_mu] * empcs->SE_P->log_conv[k][org_s_mu][new_s_mu];
				}
				psi_k[org_s_mu] += z_k_gap * empcs->SE_P->log_conv[k][org_s_mu][fix_s_mu];	/* For gap. */
			}

			org_s_mu = empcs->Mu[k][l];
			tmp_psi_k = psi_k[org_s_mu];
			flag = 0;
			for(new_s_mu = 0; new_s_mu < Ncode; new_s_mu++){
				if(new_s_mu == org_s_mu){
					continue;
				}
				if(psi_k[new_s_mu] > tmp_psi_k){
					empcs->Mu[k][l] = new_s_mu;
					tmp_psi_k = psi_k[new_s_mu];
					flag |= 1;
				}
			}

			/* Dynamically update empcs->count_Mu_X: See non-gap version for details. */
			if(flag){
				for(n_X = 0; n_X < N_X; n_X++){
					if(empcs->X[n_X][l] == MISSING_ALLELE){
						continue;
					}
					if(empcs->X[n_X][l] != empcs->gap_index){
						empcs->count_Mu_X[n_X][k][org_s_mu][empcs->X[n_X][l]]--;
						empcs->count_Mu_X[n_X][k][empcs->Mu[k][l]][empcs->X[n_X][l]]++;
					} else{	/* For gap. */
						empcs->count_Mu_X_gap[n_X][k][org_s_mu]--;
						empcs->count_Mu_X_gap[n_X][k][empcs->Mu[k][l]]++;
					}
				}
			}
		}
	}
} /* End of Update_Mu_given_QA_skip_non_seg_gap_se_convolution(). */




/* Negative of profile log-likelihood for minimizing, -R(Mu, QA, Tt, SE_P).
 * parameters: eta_1, ..., eta_(K-1), QA(model)
 * model JC69: Tt.
 *       K80: kappa, Tt.
 *       HKY85: pi_A, pi_G, pi_C, kappa, Tt.
 * vect stores parameters: QA(model).
 * vect stores parameters: SE_P. */
double negative_logpL_Mu_given_QA_se_convolution(int m, double *vect, void *ex){
	double ret = 0.0;
	ex_struct *in = (ex_struct*) ex;

	double *vect_pt = vect + in->QA->total_n_param;
	double *org_vect_pt = in->org_vect + in->QA->total_n_param;

	in->QA->Convert_Q_matrix_array_to_vect(in->QA, in->org_vect);		/* Backup vect in the previous step. */
	in->QA->Convert_vect_to_Q_matrix_array(vect, in->QA);			/* Update QA, Tt and check_param. */
	in->empcs->SE_P->Convert_f_err_to_vect(in->empcs->SE_P, org_vect_pt);	/* Backup vect in the previous step. */
	in->empcs->SE_P->Convert_vect_to_f_err(vect_pt, in->empcs->SE_P);	/* Update SE_P->f_err. */

	/* empcs->Z_normalized is fixed and should not be updated in M-step.
	 * QA and Tt are updated by NM. Mu are updated given QA, Tt and empcs->Z_normalized. */
	if(in->QA->check_param && in->empcs->SE_P->check_param){
		in->QA->Update_log_Pt(in->QA);					/* Update log(P(t)). */
		in->EMFP->Update_Mu_given_QA(in->empcs, in->QA, in->QA_H);	/* Update Mu. */
		in->EMFP->Update_Z_modified(in->empcs, in->QA);			/* Update empcs->Z_modified (log and unnormalized). */
		ret = -in->EMFP->Compute_R(in->empcs, in->QA, in->QA_H);	/* Compute R(Mu, QA, Tt). */
	} else{
		/* NM failed, restore to the original vect and return Inf to stop NM. */
		in->QA->Convert_vect_to_Q_matrix_array(in->org_vect, in->QA);		/* Restore to the original vect. */
		in->empcs->SE_P->Convert_vect_to_f_err(org_vect_pt, in->empcs->SE_P);	/* Restore to the original vect. */
		ret = Inf;
	}

	#if (EMDEBUG & 8) == 8
		printf("    Update Mu given QA\n");
		printf("      vect:");
		print_vect(m, vect);
		if(is_finite(ret)){
			printf("    -logpL: %.4f\n", ret);
		} else{
			printf("    -logpL: %.4e\n", ret);
		}
	#endif

	return(ret);
} /* End of negative_logpL_Mu_given_QA_se_convolution(). */


/* Maximize profile complete log-likelihood R(Mu, QA, Tt). */
int Maximize_logpL_se_convolution(em_phyclust_struct *empcs, Q_matrix_array *QA, Q_matrix_array *QA_H, em_control *EMC, em_fp *EMFP){
	ex_struct exs;
	nm_struct *nms;
	int ret_stop;

	/* (1, ..., QA->total_n_param) are for evolution models QA,
	 * (QA->total_n_param + 1, ..., mempcs->SE_P->n_param) are for sequencing error models. */
	int total_param = QA->total_n_param + empcs->SE_P->n_param;
	double *vect = allocate_double_1D(total_param);
	double *org_vect = allocate_double_1D(total_param);
	double *vect_pt = vect + QA->total_n_param;

	QA->Convert_Q_matrix_array_to_vect(QA, vect);
	empcs->SE_P->Convert_f_err_to_vect(empcs->SE_P, vect_pt);

	exs.empcs = empcs;
	exs.EMFP = EMFP;
	exs.QA = QA;
	exs.QA_H = QA_H;
	exs.org_vect = org_vect;

	nms = initialize_nm_struct(total_param);
	nms->Bvec = vect;
	nms->ex = &exs;

	nms->fminfn = &negative_logpL_Mu_given_QA_se_convolution;
	nms->abstol = EMC->nm_abstol_Mu_given_QA;
	nms->reltol = EMC->nm_reltol_Mu_given_QA;
	nms->maxit = EMC->nm_maxit_Mu_given_QA;

	#if (EMDEBUG & 4) == 4
		int i, j;
		double tmp_i, tmp_logpL, tmp_R, tmp_Q, tmp_obs;
		tmp_logpL = -nms->fminfn(nms->n_param, nms->Bvec, nms->ex);
		tmp_R = EMFP->LogL_profile(empcs, QA, QA_H);
		tmp_Q = EMFP->LogL_complete(empcs, QA, QA_H);
		tmp_obs = EMFP->LogL_observed(empcs, QA);
		printf("    maximize_logpL:\n");
		printf("    init. coefs:");
		print_vect(nms->n_param, vect);
		if(is_finite(tmp_logpL) && is_finite(tmp_R) && is_finite(tmp_Q) && is_finite(tmp_obs)){
			printf("        logpL: %.6f\n", tmp_logpL);
			printf("        R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
			for(i = 0 ; i < empcs->N_X; i++){
				tmp_i = 0.0;
				for(j = 0; j < empcs->K; j++){
					tmp_i += empcs->Z_normalized[i][j] * empcs->log_Eta[j];
				}
				if(empcs->replication_X[i] == 1){
					tmp_logpL += tmp_i;
				} else{
					tmp_logpL += tmp_i * empcs->replication_X[i];
				}
			}
			printf("        R+head terms = %.6f. [indep fcn]\n", tmp_logpL);
		} else{
			printf("        logpL: %.6e\n", tmp_logpL);
			printf("        R = %.6e, Q = %.6e, Obs = %.6e. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		}
	#endif

	#if (EMDEBUG & 32) == 32
		printf("    maximize_logpL:\n");
		printf("    init. f_err:");
		print_SE_P(empcs->SE_P);
	#endif

	ret_stop = phyclust_optim_nmmin(nms);
	if(ret_stop > 0){
		return(ret_stop);
	}

	#if (EMDEBUG & 4) == 4
		tmp_logpL = -nms->fminfn(nms->n_param, nms->Bvec, nms->ex);
		tmp_R = EMFP->LogL_profile(empcs, QA, QA_H);
		tmp_Q = EMFP->LogL_complete(empcs, QA, QA_H);
		tmp_obs = EMFP->LogL_observed(empcs, QA);
		printf("    conv. coefs:");
		print_vect(nms->n_param, vect);
		if(is_finite(tmp_logpL) && is_finite(tmp_R) && is_finite(tmp_Q) && is_finite(tmp_obs)){
			printf("        logpL: %.6f\n", tmp_logpL);
			printf("        R = %.6f, Q = %.6f, Obs = %.6f. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
			for(i = 0 ; i < empcs->N_X; i++){
				tmp_i = 0.0;
				for(j = 0; j < empcs->K; j++){
					tmp_i += empcs->Z_normalized[i][j] * empcs->log_Eta[j];
				}
				if(empcs->replication_X[i] == 1){
					tmp_logpL += tmp_i;
				} else{
					tmp_logpL += tmp_i * empcs->replication_X[i];
				}
			}
			printf("        R+head terms = %.6f.\n", tmp_logpL);
		} else{
			printf("        logpL: %.6e.\n", tmp_logpL);
			printf("        R = %.6e, Q = %.6e, Obs = %.6e. [indep fcn]\n", tmp_R, tmp_Q, tmp_obs);
		}
		printf("    convergence: %d\n", *nms->fail);
	#endif

	#if (EMDEBUG & 32) == 32
		printf("    conv. f_err:");
		print_SE_P(empcs->SE_P);
	#endif
	
	EMC->converge_inner_iter += *nms->fncount;

	free(vect);
	free(org_vect);
	free_nm_struct(nms);
	return(ret_stop);
} /* End of Maximize_logpL_se_convolution(). */


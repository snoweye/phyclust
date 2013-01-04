/* This file contains all fucntions for initialization. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_init_procedure.h"
#include "phyclust_init_method.h"
#include "phyclust_em_tool.h"

/*
 * Possible combinations checked in init_em_step() only.
 *
 * CODE_TYPE	SUBSTITUTION		EDIST
 * ---------	------------		-----
 * NUCLEOTIDE	JC69,K80,F81,HKY85	D_JC69,D_K80,D_HAMMING
 * SNP		SNP_JC69,SNP_F81	D_HAMMING
 *
 * INIT_METHOD	INIT_PROCEDURE		IDENTIFIER	LABEL_METHOD
 * -----------	--------------		----------	------------
 * randomMu	all			all		all
 * NJ		exhaustEM(iter=1)	all		NONE
 * randomNJ	all			all		NONE
 * PAM		exhaustEM(iter=1)	all		NONE
 * kMedoids	all			all		NONE
 * manualMu	exhaustEM(iter=1)	all		NONE
 * sampleL	all			all		NONE
 */

/* Collection of all initialization procedures. */
void init_em_step(phyclust_struct *pcs, Q_matrix_array *QA, em_control *EMC, em_fp *EMFP){
	double lower_bound_org, lower_bound;

	if(pcs->K * EMC->min_n_class >= pcs->N_X){
		fprintf_stderr("PE: K is too huge.\n");
		exit(1);
	}
	lower_bound_org = (double) EMC->min_n_class / (double) pcs->N_X_org;
	lower_bound = 1.0 / (double) pcs->N_X;
	EMC->Eta_lower_bound = (lower_bound_org > lower_bound) ? lower_bound_org : lower_bound;
	if(pcs->K > 1){
		EMC->Eta_upper_bound = 1.0 - EMC->Eta_lower_bound;
	} else{
		EMC->Eta_upper_bound = 1.0;
	}

	/* Check exceptions. */
	if(pcs->label->label_method != NONE){
		EMC->init_method = randomMu;
	}
	if(EMC->init_method == randomNJ && pcs->K == 1){
		EMC->exhaust_iter = 1;
		EMC->init_procedure = exhaustEM;
	}
	update_em_control(EMC);
	update_Q_matrix_array(QA, pcs);

	#if (EMDEBUG & 1) == 1
		printf("init proc: %s, method: %s, sub model: %s, dist: %s.\n",
		INIT_PROCEDURE[EMC->init_procedure], INIT_METHOD[EMC->init_method],
		SUBSTITUTION_MODEL[EMC->substitution_model], EDISTANCE_MODEL[EMC->edist_model]);
	#endif
	switch(EMC->init_procedure){
		case exhaustEM:
			exhaust_EM(pcs, QA, EMC, EMFP);
			break;
		case emEM:
			em_EM(pcs, QA, EMC, EMFP);
			break;
		case RndEM:
			Rnd_EM(pcs, QA, EMC, EMFP);
			break;
		case RndpEM:
			Rnd_EM(pcs, QA, EMC, EMFP);
			break;
		default:
			fprintf_stderr("PE: The initial procedure is not found.\n");
			exit(1);
	}
} /* End of init_em_step(). */


/* For models using empirical pi's in Q. */
void update_Q_matrix_array(Q_matrix_array *QA, phyclust_struct *pcs){
	int s, n_X_org, l, k, flag = 0;
	double pi[QA->ncode], total = 0.0, sum = 0.0;

	for(s = 0; s < QA->ncode; s++){
		pi[s] = 0;
	}
	for(n_X_org = 0; n_X_org < pcs->N_X_org; n_X_org++){
		for(l = 0; l < pcs->L; l++){
			if(pcs->X_org[n_X_org][l] == pcs->gap_index){	/* For gaps. */
				continue;
			}
			pi[pcs->X_org[n_X_org][l]]++;
			total++;
		}
	}
	for(s = 0; s < (QA->ncode - 1); s++){
		pi[s] = pi[s] / total;
		sum += pi[s];
	}
	pi[s] = 1.0 - sum;
	for(s = 0; s < QA->ncode; s++){
		flag |= ((pi[s] <= QA->lower_bound) || (pi[s] >= QA->upper_bound));
	}
	if(flag){
		fprintf_stderr("PE: Empirical pi's:");
		for(s = 0; s < QA->ncode; s++){
			printf(" %e", pi[s]);
		}
		printf("\n");
		exit(1);
	} else{
		for(k = 0; k < QA->K; k++){
			for(s = 0; s < QA->ncode; s++){
				QA->Q[k]->pi[s] = pi[s];
			}
		}
		QA->Update_log_Pt(QA);
	}
} /* End of update_Q_matrix_array(). */




/* Initialization procedure. */
void exhaust_EM(phyclust_struct *pcs, Q_matrix_array *org_QA, em_control *org_EMC, em_fp *EMFP){
	int init_flag = 0, init_iter = 0;
	int converge_iter = 0, converge_inner_iter = 0, converge_cm_iter = 0;
	int iter = 1, exhaust_iter = org_EMC->exhaust_iter;
	Q_matrix_array *new_QA;
	em_control *new_EMC;
	em_phyclust_struct *org_empcs, *new_empcs;

	new_QA = duplicate_Q_matrix_array(org_QA);
	new_EMC = duplicate_em_control(org_EMC);
	org_empcs = initialize_em_phyclust_struct(pcs);
	new_empcs = initialize_em_phyclust_struct(pcs);

	#if verbosity_exhaust_EM > 0
		printf("Iteration %d/%d", iter, org_EMC->exhaust_iter);
	#endif
	#if verbosity_exhaust_EM > 1
		printf(":\n\t");
	#elif verbosity_exhaust_EM > 0
		printf("\n");
	#endif
	init_flag = EMFP->Update_init(new_empcs, new_QA, new_EMC, EMFP);

	if(exhaust_iter == 1 && init_flag > 0){
	 	free_Q_matrix_array(new_QA);
		free_em_control(new_EMC);
		free_em_phyclust_struct(org_empcs);
		free_em_phyclust_struct(new_empcs);
		fprintf_stderr("PE: Initialization error. (%s, %s)\n", INIT_PROCEDURE[org_EMC->init_procedure],
				INIT_METHOD[org_EMC->init_method]);
		exit(1);
	}

	EMFP->Em_step(new_empcs, new_QA, new_EMC, EMFP);
	EMFP->Copy_empcs(new_empcs, org_empcs);
	org_QA->Copy_Q_matrix_array(new_QA, org_QA);

	copy_EMC(new_EMC, org_EMC);
	converge_iter += new_EMC->converge_iter;
	converge_inner_iter += new_EMC->converge_inner_iter;
	converge_cm_iter += new_EMC->converge_cm_iter;

	for(iter = 1; iter < exhaust_iter; iter++){
		#if verbosity_exhaust_EM > 0
			printf("\nIteration %d/%d", iter+1, org_EMC->exhaust_iter);
		#endif
		#if verbosity_exhaust_EM > 1
			printf(":\n\t");
		#elif verbosity_exhaust_EM > 0
			printf("\n");
		#endif
		init_flag = 1;
		init_iter = 0;
		while(init_flag > 0 && init_iter < org_EMC->max_init_iter){
			init_flag = EMFP->Update_init(new_empcs, new_QA, new_EMC, EMFP);
			init_iter++;
		}
		if(init_flag > 0){
			iter++;
			continue;
		}

		EMFP->Em_step(new_empcs, new_QA, new_EMC, EMFP);
		converge_iter += new_EMC->converge_iter;
		converge_inner_iter += new_EMC->converge_inner_iter;
		converge_cm_iter += new_EMC->converge_cm_iter;

		if(new_empcs->logL_observed > org_empcs->logL_observed &&
		   new_EMC->converge_flag < 2){
			EMFP->Copy_empcs(new_empcs, org_empcs);
			org_QA->Copy_Q_matrix_array(new_QA, org_QA);
			copy_EMC(new_EMC, org_EMC);
		}
	}

	if(org_empcs->logL_observed == -Inf){
		free_Q_matrix_array(new_QA);
		free_em_control(new_EMC);
		free_em_phyclust_struct(org_empcs);
		free_em_phyclust_struct(new_empcs);
		fprintf_stderr("PE: Initialization error. (%s, %s)\n", INIT_PROCEDURE[org_EMC->init_procedure],
				INIT_METHOD[org_EMC->init_method]);
		exit(1);
	}

	org_EMC->converge_iter = converge_iter;
	org_EMC->converge_inner_iter = converge_inner_iter;
	org_EMC->converge_cm_iter = converge_cm_iter;
	#if verbosity_exhaust_EM > 0
		printf("\n");
	#endif

	/* Copy results from empcs to pcs. */
	EMFP->Copy_empcs_to_pcs(org_empcs, pcs);

	free_Q_matrix_array(new_QA);
	free_em_control(new_EMC);
	free_em_phyclust_struct(org_empcs);
	free_em_phyclust_struct(new_empcs);
} /* End of exhaust_EM(). */


void Rnd_EM(phyclust_struct *pcs, Q_matrix_array *org_QA, em_control *org_EMC, em_fp *EMFP){
	int init_flag = 0, init_iter = 0;
	int converge_iter = 0, converge_inner_iter = 0, converge_cm_iter = 0;
	int iter = 0, short_iter = org_EMC->short_iter, EM_iter = org_EMC->EM_iter;
	double EM_eps = org_EMC->EM_eps;
	Q_matrix_array *new_QA;
	em_control *new_EMC;
	em_phyclust_struct *org_empcs, *new_empcs;

	new_QA = duplicate_Q_matrix_array(org_QA);

	org_EMC->EM_iter = 1;
	org_EMC->EM_eps = Inf;
	new_EMC = duplicate_em_control(org_EMC);
	org_empcs = initialize_em_phyclust_struct(pcs);
	new_empcs = initialize_em_phyclust_struct(pcs);

	org_empcs->logL_observed = -Inf;
	for(iter = 0; iter < short_iter; iter++){
		init_flag = 1;
		init_iter = 0;
		while(init_flag > 0 && init_iter < org_EMC->max_init_iter){
			init_flag = EMFP->Update_init(new_empcs, new_QA, new_EMC, EMFP);
			init_iter++;
		}
		if(init_flag > 0){
			iter++;
			continue;
		}

		EMFP->Em_step(new_empcs, new_QA, new_EMC, EMFP);
		converge_iter += new_EMC->converge_iter;
		converge_inner_iter += new_EMC->converge_inner_iter;
		converge_cm_iter += new_EMC->converge_cm_iter;

		if(new_empcs->logL_observed > org_empcs->logL_observed){
			EMFP->Copy_empcs(new_empcs, org_empcs);
			org_QA->Copy_Q_matrix_array(new_QA, org_QA);
			copy_EMC(new_EMC, org_EMC);
		}
	}

	if(org_empcs->logL_observed == -Inf){
		free_Q_matrix_array(new_QA);
		free_em_control(new_EMC);
		free_em_phyclust_struct(org_empcs);
		free_em_phyclust_struct(new_empcs);
		fprintf_stderr("PE: Initialization error. (%s, %s)\n", INIT_PROCEDURE[org_EMC->init_procedure],
				INIT_METHOD[org_EMC->init_method]);
		exit(1);
	}

	org_EMC->EM_iter = EM_iter;
	org_EMC->EM_eps = EM_eps;
	EMFP->Em_step(org_empcs, org_QA, org_EMC, EMFP);
	org_EMC->converge_iter += converge_iter;
	org_EMC->converge_inner_iter += converge_inner_iter;
	org_EMC->converge_cm_iter += converge_cm_iter;

	/* Copy results from empcs to pcs. */
	EMFP->Copy_empcs_to_pcs(org_empcs, pcs);

	free_Q_matrix_array(new_QA);
	free_em_control(new_EMC);
	free_em_phyclust_struct(org_empcs);
	free_em_phyclust_struct(new_empcs);
} /* End of Rnd_EM(). */


void Rndp_EM(phyclust_struct *pcs, Q_matrix_array *org_QA, em_control *org_EMC, em_fp *EMFP){
	int init_flag = 0, init_iter = 0;
	int converge_iter = 0, converge_inner_iter = 0, converge_cm_iter = 0;
	int iter = 0, short_iter = org_EMC->short_iter, EM_iter = org_EMC->EM_iter, fixed_iter = org_EMC->fixed_iter;
	double EM_eps = org_EMC->EM_eps;
	Q_matrix_array *new_QA;
	em_control *new_EMC;
	em_phyclust_struct *org_empcs, *new_empcs;

	new_QA = duplicate_Q_matrix_array(org_QA);

	org_EMC->EM_iter = fixed_iter;
	org_EMC->EM_eps = Inf;
	new_EMC = duplicate_em_control(org_EMC);
	org_empcs = initialize_em_phyclust_struct(pcs);
	new_empcs = initialize_em_phyclust_struct(pcs);

	org_empcs->logL_observed = -Inf;
	while(iter < short_iter){
		init_flag = 1;
		init_iter = 0;
		while(init_flag > 0 && init_iter < org_EMC->max_init_iter){
			init_flag = EMFP->Update_init(new_empcs, new_QA, new_EMC, EMFP);
			init_iter++;
		}
		if(init_flag > 0){
			iter++;
			continue;
		}

		EMFP->Em_step(new_empcs, new_QA, new_EMC, EMFP);
		converge_iter += new_EMC->converge_iter;
		converge_inner_iter += new_EMC->converge_inner_iter;
		converge_cm_iter += new_EMC->converge_cm_iter;

		if(new_empcs->logL_observed > org_empcs->logL_observed){
			EMFP->Copy_empcs(new_empcs, org_empcs);
			org_QA->Copy_Q_matrix_array(new_QA, org_QA);
			copy_EMC(new_EMC, org_EMC);
		}
		iter += fixed_iter;
	}

	if(org_empcs->logL_observed == -Inf){
		free_Q_matrix_array(new_QA);
		free_em_control(new_EMC);
		free_em_phyclust_struct(org_empcs);
		free_em_phyclust_struct(new_empcs);
		fprintf_stderr("PE: Initialization error. (%s, %s)\n", INIT_PROCEDURE[org_EMC->init_procedure],
				INIT_METHOD[org_EMC->init_method]);
		exit(1);
	}

	org_EMC->EM_iter = EM_iter;
	org_EMC->EM_eps = EM_eps;
	EMFP->Em_step(org_empcs, org_QA, org_EMC, EMFP);
	org_EMC->converge_iter += converge_iter;
	org_EMC->converge_inner_iter += converge_inner_iter;
	org_EMC->converge_cm_iter += converge_cm_iter;

	/* Copy results from empcs to pcs. */
	EMFP->Copy_empcs_to_pcs(org_empcs, pcs);

	free_Q_matrix_array(new_QA);
	free_em_control(new_EMC);
	free_em_phyclust_struct(org_empcs);
	free_em_phyclust_struct(new_empcs);
} /* End of Rndp_EM(). */


void em_EM(phyclust_struct *pcs, Q_matrix_array *org_QA, em_control *org_EMC, em_fp *EMFP){
	int init_flag = 0, init_iter = 0;
	int converge_iter = 0, converge_inner_iter = 0, converge_cm_iter = 0;
	int short_iter = org_EMC->short_iter;
	double short_eps = org_EMC->short_eps;
	Q_matrix_array *new_QA;
	em_control *new_EMC;
	em_phyclust_struct *org_empcs, *new_empcs;

	new_QA = duplicate_Q_matrix_array(org_QA);
	new_EMC = duplicate_em_control(org_EMC);
	org_empcs = initialize_em_phyclust_struct(pcs);
	new_empcs = initialize_em_phyclust_struct(pcs);

	org_empcs->logL_observed = -Inf;
	while(new_EMC->short_iter > 0){
		init_flag = 1;
		init_iter = 0;
		while(init_flag > 0 && init_iter < org_EMC->max_init_iter){
			init_flag = EMFP->Update_init(new_empcs, new_QA, new_EMC, EMFP);
			init_iter++;
		}
		if(init_flag > 0){
			new_EMC->short_iter--;
			continue;
		}

		EMFP->Short_em_step(new_empcs, new_QA, new_EMC, EMFP);
		converge_iter += new_EMC->converge_iter;
		converge_inner_iter += new_EMC->converge_inner_iter;
		converge_cm_iter += new_EMC->converge_cm_iter;

		if(new_empcs->logL_observed > org_empcs->logL_observed){
			EMFP->Copy_empcs(new_empcs, org_empcs);
			org_QA->Copy_Q_matrix_array(new_QA, org_QA);
			copy_EMC(new_EMC, org_EMC);
		}
		new_EMC->short_iter -= new_EMC->converge_iter;
	}

	if(org_empcs->logL_observed == -Inf){
		free_Q_matrix_array(new_QA);
		free_em_control(new_EMC);
		free_em_phyclust_struct(org_empcs);
		free_em_phyclust_struct(new_empcs);
		fprintf_stderr("PE: Initialization error. (%s, %s)\n", INIT_PROCEDURE[org_EMC->init_procedure],
				INIT_METHOD[org_EMC->init_method]);
		exit(1);
	}

	org_EMC->short_iter = short_iter;
	org_EMC->short_eps = short_eps;
	EMFP->Em_step(org_empcs, org_QA, org_EMC, EMFP);
	org_EMC->converge_iter += converge_iter;
	org_EMC->converge_inner_iter += converge_inner_iter;
	org_EMC->converge_cm_iter += converge_cm_iter;

	/* Copy results from empcs to pcs. */
	EMFP->Copy_empcs_to_pcs(org_empcs, pcs);

	free_Q_matrix_array(new_QA);
	free_em_control(new_EMC);
	free_em_phyclust_struct(org_empcs);
	free_em_phyclust_struct(new_empcs);
} /* End of em_EM(). */


#include <R.h>
#include <Rinternals.h>

SEXP R_phyclust(SEXP R_N_X_org, SEXP R_L, SEXP R_K, SEXP R_X, SEXP R_EMC, SEXP R_manual_id,
		SEXP R_label);
SEXP R_phyclust_edist(SEXP R_edist_model, SEXP R_N_X, SEXP R_L, SEXP R_X);
SEXP R_phyclust_em_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label);
SEXP R_phyclust_e_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_Z_state, SEXP R_label);
SEXP R_phyclust_m_step(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_vect, SEXP R_Z_normalized,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label);
SEXP R_phyclust_find_consensus(SEXP R_N_X_org, SEXP R_L, SEXP R_code_type,
		SEXP R_WIGAP, SEXP R_X_org);
SEXP R_phyclust_logL(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_K,
		SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_substitution_model, SEXP R_identifier, SEXP R_code_type,
		SEXP R_label);
SEXP R_phyclust_logPt(SEXP R_pi, SEXP R_kappa, SEXP R_Tt,
		SEXP R_code_type, SEXP R_substitution_model);
SEXP R_phyclust_se(SEXP R_N_X_org, SEXP R_L, SEXP R_K, SEXP R_X, SEXP R_EMC,
		SEXP R_manual_id, SEXP R_label);
SEXP R_phyclust_se_update(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_EMC,
		SEXP R_K, SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_label);
SEXP R_phyclust_update(SEXP R_N_X_org, SEXP R_L, SEXP R_X, SEXP R_EMC,
		SEXP R_K, SEXP R_Eta, SEXP R_Mu, SEXP R_vect,
		SEXP R_label);
SEXP R_RRand(SEXP R_N, SEXP R_TRUK, SEXP R_PREDK, SEXP R_trcl, SEXP R_prcl);

SEXP R_ms_main(SEXP R_argv, SEXP R_ms_file);
SEXP R_paml_baseml_main(SEXP R_argv, SEXP R_file_name);
SEXP R_seq_gen_main(SEXP R_argv, SEXP R_seq_gen_file_name);


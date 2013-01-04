/* This file contains constants for phyclust. */

#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


#ifndef __PHYCLUST_CONSTANT_
#define __PHYCLUST_CONSTANT_

#include <float.h>
#define Inf DBL_MAX

#ifdef __GNUC__
#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
#else
#define VARIABLE_IS_NOT_USED
#endif


/* Nucleotides. */
#define NN 4								/* Number of nucleotides. */
#define NNG 5								/* Index of nucleotides with gap. */
enum {A, G, C, T, GAP};							/* Nucleotides. */
static const char NUCLEOTIDE_CODE[NNG]  = {'A', 'G', 'C', 'T', '-'};
static const char NUCLEOTIDE_lower[NNG]  = {'a', 'g', 'c', 't', '-'};
static const char NUCLEOTIDE_ID[NNG]    = {'0', '1', '2', '3', '4'};

/* SNPs. */
#define NSNP 2								/* Number of SNPs. */
#define NSNPG 3								/* Index of SNPs with unknown. */
enum {NL, SU, UN};							/* Normal/Substitution/Unknown. */
static const char SNP_CODE[NSNPG] = {'0', '1', '-'};
static const char SNP_ID[NSNPG] = {'1', '2', '-'};

/* Codons. */
#define N_CODON 64
#define N_CODON_G 65
// static const char CODON_CODE[N_CODON] = {};
// static const char CODON_lower[N_CODON] = {};
// static const char CoDON_ID[N_AMINO_ACID_G] = {};

/* Amino Acids. */
#define N_AMINO_ACID 21
#define N_AMINO_ACID_G 22
// enum {NL, SU, UN};
// static const char AMINO_ACID_CODE[N_AMINO_ACID] = {};
// static const char AMINO_ACID_ID[N_AMINO_ACID_G] = {};

/* Common missing Allele for all code types. */
#define MISSING_ALLELE -1


/* CODEs. */
#define N_CODE_TYPE 4							/* Total number of code type. */
enum {NUCLEOTIDE, SNP, CODON, AMINO_ACID};				/* Number of codes. */
static char VARIABLE_IS_NOT_USED *CODE_TYPE[N_CODE_TYPE] =
	{"NUCLEOTIDE", "SNP", "CODON", "AMINO_ACID"};
static const int NCODE[N_CODE_TYPE] =
	{NN, NSNP, N_CODON, N_AMINO_ACID};				/* Indicate the dimension. */
static const int NCODE_WIGAP[N_CODE_TYPE] =
	{NNG, NSNPG, N_CODON_G, N_AMINO_ACID_G};			/* Indicate the dimension. */
#define PRINT_CODE_TYPE 0						/* 0 for CODE, 1 for CODE_ID. */
static const int GAP_INDEX[N_CODE_TYPE] =
	{NN, NSNP, N_CODON, N_AMINO_ACID};				/* Indicate the gap index. */


/* EM procedures. */
#define N_INIT_PROCEDURE 4
enum {exhaustEM, emEM, RndEM, RndpEM};					/* Initialization procedures. */
static char VARIABLE_IS_NOT_USED *INIT_PROCEDURE[N_INIT_PROCEDURE] =
	{"exhaustEM", "emEM", "RndEM", "RndpEM"};

/* Initialization methods. */
#define N_INIT_METHOD 7
enum {randomMu, NJ, randomNJ, PAM, kMedoids, manualMu, sampleL};	/* Initialization methods. */
static char VARIABLE_IS_NOT_USED *INIT_METHOD[N_INIT_METHOD] =
	{"randomMu", "NJ", "randomNJ", "PAM", "K-Medoids", "manualMu", "sampleL"};


/* Substitution models. */
#define N_SUB_MODEL 9							/* Number of evolution models. */
enum {JC69, K80, F81, HKY85, SNP_JC69, SNP_F81,
	E_F81, E_HKY85, E_SNP_F81};
static char VARIABLE_IS_NOT_USED *SUBSTITUTION_MODEL[N_SUB_MODEL] =
	{"JC69", "K80", "F81", "HKY85", "SNP_JC69", "SNP_F81",
		"E_F81", "E_HKY85", "E_SNP_F81"};

/* Distance models. */
#define N_EDIST 4
enum {D_JC69, D_K80, D_HAMMING, D_HAMMING_WOGAP};			/* Distance models. */
static char VARIABLE_IS_NOT_USED *EDISTANCE_MODEL[N_EDIST] =
	{"D_JC69", "D_K80", "D_HAMMING", "D_HAMMING_WOGAP"};

/* Sequencing error models for nucleotide only. */
#define N_SE_TYPE 2
enum {SE_NO, SE_YES};							/* Sequencing error types. */
static char VARIABLE_IS_NOT_USED *SE_TYPE[N_SE_TYPE] =
	{"SE_NO", "SE_YES"};
#define N_SE_MODEL 1
enum {SE_CONVOLUTION};							/* Seaquencing error models. */
static char VARIABLE_IS_NOT_USED *SE_MODEL[N_SE_MODEL] =
	{"SE_CONVOLUTION"};
#define SE_CONSTANT 1e-2

//#define N_SE_P_MODEL 2
//enum {SE_P_NONE, SE_CONVOLUTION, SE_P_UNSTRUCT, SE_P_SYMMETRIC};			/* Seaquencing error models. */
//static char VARIABLE_IS_NOT_USED *SE_P_MODEL[N_SE_P_MODEL] =
//	{"SE_P_NONE", "SE_P_CONVOLUTION", "SE_P_UNSTRUCT", "SE_P_SYMMETRIC"};


/* Identifier for constraints. */
#define N_IDENTIFIER 4
enum {EE, EV, VE, VV};							/* First letter for Q, second for Tt. */
static char VARIABLE_IS_NOT_USED *IDENTIFIER[N_IDENTIFIER] =
	{"EE", "EV", "VE", "VV"};


/* EM methods. */
#define N_EM_METHOD 4
enum {EM, ECM, AECM, APECM};
static char VARIABLE_IS_NOT_USED *EM_METHOD[N_EM_METHOD] =
	{"EM", "ECM", "AECM", "APECM"};


/* Boundary methods. */
#define N_BOUNDARY_METHOD 2
enum {ADJUST, IGNORE};
static char VARIABLE_IS_NOT_USED *BOUNDARY_METHOD[N_BOUNDARY_METHOD] =
	{"ADJUST", "IGNORE"};


/* Label methods. */
#define N_LABEL_METHOD 3
enum {NONE, SEMI, GENERAL};
static char VARIABLE_IS_NOT_USED *LABEL_METHOD[N_LABEL_METHOD] =
	{"NONE", "SEMI", "GENERAL"};


/* EM Debugging only. */
/* Ex: 1 for em steps, 2 for m steps, 3 for em and m steps. */
#define EMDEBUG 0		/* 0  = 000000 for no output,
				 * 1  = 000001 for em steps,
				 * 2  = 000010 for m step,
				 * 4  = 000100 for nm step,
				 * 8  = 001000 for -logpL.
				 * 16 = 010000 for log_conv in se_convolution.
				 * 32 = 100000 for f_err in se_convolution. */
#define INITDEBUG 0		/* 1 */		/* 0 for no output, >0 for tracing initialization. */
#define verbosity_em_step 0	/* 3 */		/* 0 for no output, >0 for print information. */
#define verbosity_exhaust_EM 0	/* 1 */		/* 0 for no output, >0 for print information. */
#define PRINT_ERROR 0		/* 1 */		/* 0 for no output, >0 for error messages. */


/* For error messages. */
#define fprintf_stderr(...) fprintf(stderr, __VA_ARGS__)


/* Deal with R. */
#ifdef __HAVE_R_

#undef printf
#define printf Rprintf

#undef exit
#define exit(a) error("%d\n", a)

#undef fprintf_stderr
#define fprintf_stderr REprintf

#endif


#endif	/* End of __PHYCLUST_CONSTANT_. */

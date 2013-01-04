/* This file contains declarations for phylip and fasta. */

#ifndef __PHYCLUST_FILE_INPUT_
#define __PHYCLUST_FILE_INPUT_

static const int NAME_LENGTH = 10;		/* length of sequence name. */


typedef struct _input_struct		input_struct;


/* Input structure for storing data and results. */
struct _input_struct{
	/* Define code type. */
	int	code_type;			/* NUCLEOTIDE/SNP. */
	int	ncode;				/* = NN or NSNP, indicates the dimension. */
	/* Constant points to input_struct *ins. */
	int	N_X_org;			/* Number of sequences. */
	int	L;				/* Number of loci. */
	int	**X_org;			/* Data, dim = N_X_org * L. */
	char	**X_name;			/* Name, dim = N_X_org * (NAME_LENGTH + 1). */
};

input_struct* initialize_input_struct(int code_type, int N_X, int L);
void free_input_struct(input_struct *ins);

/* File: phyclust_file_input.c */
int nucleotide_to_id(char x);
input_struct* read_input_phylip(char *file_name);
input_struct* read_input_fasta(char *file_name);
input_struct* read_input_snp(char *file_name);

#endif	/* End of __PHYCLUST_FILE_INPUT_. */

/* This file contains functions to call seq_gen_main() in "seq-gen.c". */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

FILE *R_seq_gen_file_pointer;
int seq_gen_main(int argc, char **argv);

SEXP R_seq_gen_main(SEXP R_argv, SEXP R_seq_gen_file_name){
	int argc, i;
	const char **argv;
	const char *seq_gen_file_name;

	/* Obtain arguments. */
	argc = length(R_argv);
	argv = (const char**) malloc(argc * sizeof(char *));
	if(argv == NULL){
		error("Memory allocation fails!\n");
	}

	for(i = 0; i < argc; i++){
		argv[i] = CHAR(STRING_ELT(R_argv, i)); 
	}

	/* Open a file pointer. */
	seq_gen_file_name = CHAR(STRING_ELT(R_seq_gen_file_name, 0));
	R_seq_gen_file_pointer = fopen(seq_gen_file_name, "w");

	/* Run. */
	GetRNGstate();		/* Get the seed from R. */
	seq_gen_main(argc, (char**) argv);
	PutRNGstate();		/* Update the seed of R. */

	/* Close the file pointer and free memory. */
	fclose(R_seq_gen_file_pointer);
	free(argv);
	return(R_NilValue);
} /* End of R_seq_gen_main(). */

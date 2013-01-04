/* This file contains a wrap for R to call ms_main() in "ms.c". */

#include <R.h>
#include <Rinternals.h>
#include "ms.h"

SEXP R_ms_main(SEXP R_argv, SEXP R_ms_file){
	int argc, i;
	const char **argv;

	/* Obtain arguments. */
	argc = length(R_argv);
	argv = (const char**) malloc(argc * sizeof(char *));
	if(argv == NULL){
		error("Memory allocation fails!\n");
	}

	for(i = 0; i < argc; i++){
		argv[i] = CHAR(STRING_ELT(R_argv, i)); 
	}

	/* Get a file name, but open file inside ms_main(). */
	R_ms_file_name = CHAR(STRING_ELT(R_ms_file, 0));

	/* Required by ms_main() and segtre_mig(). */
	maxsites = SITESINC;
	seglimit = SEGINC;

	/* Run. */
	GetRNGstate();		/* Get the seed from R. */
	ms_main(argc, (char**) argv);
	PutRNGstate();		/* Update the seed of R. */

	/* Free memory. */
	free(argv);
	return(R_NilValue);
} /* End of R_ms_main(). */

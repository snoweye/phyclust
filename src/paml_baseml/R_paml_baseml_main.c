/* This file contains functions to call paml_baseml_main(). in "baseml.c". */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
//WCC #include <unistd.h>

FILE *R_paml_baseml_file_pointer;
int paml_baseml_main(int argc, char **argv);

//WCC SEXP R_paml_baseml_main(SEXP R_argv, SEXP R_temp_dir_name){
SEXP R_paml_baseml_main(SEXP R_argv, SEXP R_file_name){
	int argc, i;
	const char **argv;
	const char *file_name;
//WCC	const char *temp_dir_name;

	argc = length(R_argv);
	argv = (const char**) malloc(argc * sizeof(char *));
	if(argv == NULL){
		error("Memory allocation fails!\n");
	}

	for(i = 0; i < argc; i++){
		argv[i] = CHAR(STRING_ELT(R_argv, i)); 
	}

	file_name = CHAR(STRING_ELT(R_file_name, 0));
	R_paml_baseml_file_pointer = fopen(file_name, "w");

/* WCC
	temp_dir_name = CHAR(STRING_ELT(R_temp_dir_name, 0));
	if(chdir(temp_dir_name) != 0){
		printf("The directory is not found.\n");
		exit(1);
	}
*/

	GetRNGstate();		/* Get the seed from R. */
	paml_baseml_main(argc, (char**) argv);
	PutRNGstate();		/* Update the seed of R. */

	fclose(R_paml_baseml_file_pointer);
	free(argv);
	return(R_NilValue);
} /* End of R_paml_baseml_main(). */

/* This file contains functions for memory allocations. */

#include <stdlib.h>
#include <stdio.h>
#include "phyclust_constant.h"
#include "phyclust_tool.h"


/* Allocate a pointer array with double precision. */
double** allocate_double_2D_AP(int n_X){
	int i;
	double **pointerarray = (double **) malloc(n_X * sizeof(double *));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		pointerarray[i] = NULL;
	}

	return(pointerarray);
} /* End of allocate_double_2D_AP(). */

/* Allocate a pointer array with integer precision. */
int** allocate_int_2D_AP(int n_X){
	int i;
	int **pointerarray = (int **) malloc(n_X * sizeof(int *));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		pointerarray[i] = NULL;
	}

	return(pointerarray);
} /* End of allocate_int_2D_AP(). */

/* Allocate a pointer array with char. */
char** allocate_char_2D_AP(int n_X){
	int i;
	char **pointerarray = (char **) malloc(n_X * sizeof(char *));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		pointerarray[i] = NULL;
	}

	return(pointerarray);
} /* End of allocate_char_2D_AP(). */




/* Allocate a pointer array with double precision. */
double* allocate_double_1D(int n_X){
	int i;
	double *array = (double *) malloc(n_X * sizeof(double));

	if(array == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		array[i] = 0.0;
	}

	return(array);
} /* End of allocate_double_1D(). */

/* Allocate a pointer array with integer precision. */
int* allocate_int_1D(int n_X){
	int i;
	int *array = (int *) malloc(n_X * sizeof(int));

	if(array == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		array[i] = 0;
	}

	return(array);
} /* End of allocate_int_1D(). */

/* Allocate a pointer array with char. */
char* allocate_char_1D(int n_X){
	int i;
	char *array = (char *) malloc((n_X + 1) * sizeof(char));

	if(array == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(i = 0; i < n_X; i++){
		array[i] = '0';
	}
	array[n_X] = '\0';

	return(array);
} /* End of allocate_char_1D(). */




/* Allocate a square array with double. */
double** allocate_double_SQ(int n_X){
	int i, j;
	double **pointerarray = allocate_double_2D_AP(n_X);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < n_X; i++){
		pointerarray[i] = allocate_double_1D(n_X);
		if(pointerarray[i] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(j = 0; j < n_X; j++){
			pointerarray[i][j] = 0.0;
		}
	}

	return(pointerarray);
} /* End of allocate_double_SQ(). */

/* Allocate a upper triangular array with double. */
double** allocate_double_UT(int n_X){
	int i, j;
	double **pointerarray = allocate_double_2D_AP(n_X);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < n_X; i++){
		pointerarray[i] = allocate_double_1D(n_X - i);
		if(pointerarray[i] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(j = 0; j < n_X - i; j++){
			pointerarray[i][j] = 0.0;
		}
	}

	return(pointerarray);
} /* End of allocate_double_UT(). */

/* Allocate a rectriangular array with double. */
double** allocate_double_RT(int nrow, int ncol){
	int i, j;
	double **pointerarray = allocate_double_2D_AP(nrow);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < nrow; i++){
		pointerarray[i] = allocate_double_1D(ncol);
		if(pointerarray[i] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(j = 0; j < ncol; j++){
			pointerarray[i][j] = 0.0;
		}
	}

	return(pointerarray);
} /* End of allocate_double_RT(). */

/* Allocate a rectriangular array with int. */
int** allocate_int_RT(int nrow, int ncol){
	int i, j;
	int **pointerarray = allocate_int_2D_AP(nrow);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < nrow; i++){
		pointerarray[i] = allocate_int_1D(ncol);
		if(pointerarray[i] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(j = 0; j < ncol; j++){
			pointerarray[i][j] = 0;
		}
	}

	return(pointerarray);
} /* End of allocate_int_RT(). */

/* Allocate a rectriangular array with int. */
char** allocate_char_RT(int nrow, int ncol){
	int i, j;
	char **pointerarray = allocate_char_2D_AP(nrow);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < nrow; i++){
		pointerarray[i] = allocate_char_1D(ncol + 1);
		if(pointerarray[i] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(j = 0; j < ncol; j++){
			pointerarray[i][j] = '0';
		}
		pointerarray[i][ncol] = '\0';
	}

	return(pointerarray);
} /* End of allocate_int_RT(). */


void free_double_RT(int nrow, double **RT){
	int i;

	for(i = 0; i < nrow; i++){
		free(RT[i]);
	}
	free(RT);
} /* End of free_double_RT(). */

void free_int_RT(int nrow, int **RT){
	int i;

	for(i = 0; i < nrow; i++){
		free(RT[i]);
	}
	free(RT);
} /* End of free_int_RT(). */

void free_char_RT(int nrow, char **RT){
	int i;

	for(i = 0; i < nrow; i++){
		free(RT[i]);
	}
	free(RT);
} /* End of free_char_RT(). */




void copy_int_2D_AP(int length, int **from, int **to){
	int i;

	for(i = 0; i < length; i++){
		to[i] = from[i];
	}
} /* End of copy_int_2D_AP(). */

void copy_double_RT(int nrow, int ncol, double **from, double **to){
	int i, j;

	for(i = 0; i < nrow; i++){
		for(j = 0; j < ncol; j++){
			to[i][j] = from[i][j];
		}
	}
} /* End of copy_double_RT(). */

void copy_int_RT(int nrow, int ncol, int **from, int **to){
	int i, j;

	for(i = 0; i < nrow; i++){
		for(j = 0; j < ncol; j++){
			to[i][j] = from[i][j];
		}
	}
} /* End of copy_int_RT(). */

void copy_double_1D(int length, double *from, double *to){
	int i;

	for(i = 0; i < length; i++){
		to[i] = from[i];
	}
} /* End of copy_double_1D(). */

void copy_int_1D(int length, int *from, int *to){
	int i;

	for(i = 0; i < length; i++){
		to[i] = from[i];
	}
} /* End of copy_int_1D(). */




/* Allocate a sequential upper triangular array with double. */
double** allocate_s_double_UT(int n_X){
	int i, total = n_X * (n_X + 1) / 2;
	double **pointerarray = allocate_double_2D_AP(n_X);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	pointerarray[0] = allocate_double_1D(total);
	if(pointerarray[0] == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < total; i++){
		pointerarray[0][i] = 0.0;
	}

	total = n_X + 1;
	for(i = 1; i < n_X; i++){
		pointerarray[i] = pointerarray[i - 1] + (total - i);
	}

	return(pointerarray);
} /* End of allocate_s_double_UT(). */

double** allocate_s_double_LT(int n_X){
	int i, total = n_X * (n_X + 1) / 2;
	double **pointerarray = allocate_double_2D_AP(n_X);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	pointerarray[0] = allocate_double_1D(total);
	if(pointerarray[0] == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < total; i++){
		pointerarray[0][i] = 0.0;
	}

	for(i = 1; i < n_X; i++){
		pointerarray[i] = pointerarray[i - 1] + i;
	}

	return(pointerarray);
} /* End of allocate_s_double_LT(). */

/* dys structure, Kaufman and Rousseeuw (1990), p105.
 * 0
 * 9 0		=> dys = 0 9 1 4 3 2 7 ....
 * 1 4 0		 ^
 * 3 2 7 0		 The first element is 0 and then by LT form
 * ...			 due to the f2c code.
 */
double** allocate_s_double_LT_pam(int n_X){
	int i, total = n_X * (n_X + 1) / 2 + 1;
	double **pointerarray = allocate_double_2D_AP(n_X);

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	pointerarray[0] = allocate_double_1D(total);
	if(pointerarray[0] == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}

	for(i = 0; i < total; i++){
		pointerarray[0][i] = 0.0;
	}
	pointerarray[0]++;

	for(i = 1; i < n_X; i++){
		pointerarray[i] = pointerarray[i - 1] + i;
	}

	return(pointerarray);
} /* End of allocate_s_double_LT_pam(). */




/* Special storage. */
int**** allocate_int_RT_4D(int N_X, int K, int nrow, int ncol){
	int n_X, k;
	int ****pointerarray = (int ****) malloc(N_X * sizeof(int ***));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(n_X = 0; n_X < N_X; n_X++){
		pointerarray[n_X] = (int ***) malloc(K * sizeof(int **));
		if(pointerarray[n_X] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(k = 0; k < K; k++){
			pointerarray[n_X][k] = allocate_int_RT(nrow, ncol);
		}
	}

	return(pointerarray);
} /* End of allocate_int_RT_4D(). */

void free_int_RT_4D(int N_X, int K, int nrow, int ****RT4D){
	int n_X, k;

	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			free_int_RT(nrow, RT4D[n_X][k]);
		}
		free(RT4D[n_X]);
	}
	free(RT4D);
} /* End of free_int_RT_4D(). */

void copy_int_RT_4D(int N_X, int K, int nrow, int ncol, int ****from, int ****to){
	int n_X, k;

	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			copy_int_RT(nrow, ncol, from[n_X][k], to[n_X][k]);
		}
	}
} /* End of copy_int_RT_4D(). */

int*** allocate_int_RT_3D(int N_X, int K, int ncode){
	int n_X, k;
	int ***pointerarray = (int ***) malloc(N_X * sizeof(int **));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(n_X = 0; n_X < N_X; n_X++){
		pointerarray[n_X] = (int **) malloc(K * sizeof(int *));
		if(pointerarray[n_X] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(k = 0; k < K; k++){
			pointerarray[n_X][k] = allocate_int_1D(ncode);
		}
	}

	return(pointerarray);
} /* End of allocate_int_RT_3D(). */

void free_int_RT_3D(int N_X, int K, int ***RT3D){
	int n_X, k;

	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			free(RT3D[n_X][k]);
		}
		free(RT3D[n_X]);
	}
	free(RT3D);
} /* End of free_int_RT_3D(). */

void copy_int_RT_3D(int N_X, int K, int ncode, int ***from, int ***to){
	int n_X, k;

	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			copy_int_1D(ncode, from[n_X][k], to[n_X][k]);
		}
	}
} /* End of copy_int_RT_3D(). */

double*** allocate_double_RT_3D(int K, int L, int ncode){
	int k, l;
	double ***pointerarray = (double ***) malloc(K * sizeof(double **));

	if(pointerarray == NULL){
		fprintf_stderr("PE: Memory allocation fails!\n");
		exit(1);
	}
	for(k = 0; k < K; k++){
		pointerarray[k] = (double **) malloc(L * sizeof(double *));
		if(pointerarray[k] == NULL){
			fprintf_stderr("PE: Memory allocation fails!\n");
			exit(1);
		}
		for(l = 0; l < L; l++){
			pointerarray[k][l] = allocate_double_1D(ncode);
		}
	}

	return(pointerarray);
} /* End of allocate_double_RT_3D(). */

void free_double_RT_3D(int K, int L, double ***RT3D){
	int k, l;

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			free(RT3D[k][l]);
		}
		free(RT3D[k]);
	}
	free(RT3D);
} /* End of free_double_RT_3D(). */

void copy_double_RT_3D(int K, int L, int ncode, double ***from, double ***to){
	int k, l;

	for(k = 0; k < K; k++){
		for(l = 0; l < L; l++){
			copy_double_1D(ncode, from[k][l], to[k][l]);
		}
	}
} /* End of copy_double_RT_3D(). */


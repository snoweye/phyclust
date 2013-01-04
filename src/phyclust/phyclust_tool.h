/* This file contains functions for memory allocations. */

#ifndef __PHYCLUST_TOOL_
#define __PHYCLUST_TOOL_

/* Array pointer. */
double** allocate_double_2D_AP(int n_X);	/* double 2D array pointer. */
int** allocate_int_2D_AP(int n_X);		/* int 2D array pointer. */
char** allocate_char_2D_AP(int n_X);		/* char 2D array pointer. */

/* Full array. */
double* allocate_double_1D(int n_X);		/* double 1D array. */
int* allocate_int_1D(int n_X);			/* int 1D array. */
char* allocate_char_1D(int n_X);		/* char 1D array. */

/* Full array. */
double** allocate_double_SQ(int n_X);				/* double square array. */
double** allocate_double_UT(int n_X);				/* double upper triangular array, n/(n-1) for w/o diagonal. */
double** allocate_double_RT(int nrow, int ncol);		/* double rectangle array. */
int** allocate_int_RT(int nrow, int ncol);			/* int rectangle array. */
char** allocate_char_RT(int nrow, int ncol);			/* char rectangle array. */
void free_double_RT(int nrow, double **RT);			/* free double RT. */
void free_int_RT(int nrow, int **RT);				/* free int RT. */
void free_char_RT(int nrow, char **RT);				/* free char RT. */

/* Copy. */
void copy_int_2D_AP(int length, int **from, int **to);					/* copy a 2D array pointer. */
void copy_double_RT(int nrow, int ncol, double **from, double **to);			/* copy a rectangle array. */
void copy_int_RT(int nrow, int ncol, int **from, int **to);				/* copy a rectangle array. */
void copy_double_1D(int length, double *from, double *to);				/* copy a 1D array. */
void copy_int_1D(int length, int *from, int *to);					/* copy a 1D array. */

/* Sequential array. */
double** allocate_s_double_UT(int n_X);		/* sequential double upper triangular array. */
double** allocate_s_double_LT(int n_X);		/* sequential double lower triangular array. */
double** allocate_s_double_LT_pam(int n_X);	/* sequential double lower triangular array for PAM. */

/* Special storage. */
int**** allocate_int_RT_4D(int N_X, int K, int nrow, int ncol);				/* int 4D RT. */
void free_int_RT_4D(int N_X, int K, int nrow, int ****RT4D);				/* free 4D RT. */
void copy_int_RT_4D(int N_X, int K, int nrow, int ncol, int ****from, int ****to);	/* copy 4D RT. */
int*** allocate_int_RT_3D(int N_X, int K, int ncode);					/* int 3D RT. */
void free_int_RT_3D(int N_X, int K, int ***RT3D);					/* free 3D RT. */
void copy_int_RT_3D(int N_X, int K, int ncode, int ***from, int ***to);			/* copy 3D RT. */
double*** allocate_double_RT_3D(int K, int L, int ncode);				/* double 3D RT. */
void free_double_RT_3D(int K, int L, double ***RT);					/* free double 3D RT. */
void copy_double_RT_3D(int K, int L, int ncode, double ***from, double ***to);		/* copy double 3D RT. */

#endif	/* End of __PHYCLUST_TOOL_. */

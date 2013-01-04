/* This file contains declarations for phyclust. */

#ifndef __PHYCLUST_EDIST_STRUCT_
#define __PHYCLUST_EDIST_STRUCT_

enum {UT, LT_pam};				/* upper/lower triangular form. */

typedef struct _edist_struct	edist_struct;

/* Parameters of transition rate matrix. */
struct _edist_struct{
	int	form;				/* UT or LT. */
	int	N_X;				/* number of observations. */
	double	**EDM;				/* edist matrix, length = N_X * (N_X - 1) / 2. */
	double  (*get_pair_edist)(edist_struct*, int, int);	/* get distance for pair sequence. */
};

edist_struct* initialize_edist_struct(int form, int N_X);
void free_edist_struct(edist_struct *eds);

/* These two functions should be objectized into edist_struct. */
edist_struct* initialize_edist_struct_UT(int edist_model, int N_X, int L, int **X);
edist_struct* initialize_edist_struct_LT_pam(int edist_model, int N_X, int L, int **X);

/* Evolution distance for X and Mu with length = L. */
double (*get_edist_D(int edist_model))(int, int*, int*);
double edist_D_JC69(int L, int *x, int *mu);
double edist_D_K80(int L, int *x, int *mu);
double edist_D_HAMMING(int L, int *x, int *mu);
double edist_D_HAMMING_WOGAP(int L, int *x, int *mu);

/* Tools. */
double get_pair_edist_UT(edist_struct *eds, int u, int v);
double get_pair_edist_LT_pam(edist_struct *eds, int u, int v);
/* TODO: */ double* get_edist_given_mu(int edist_model, int L, int **X, int *mu);

/* ----- For debug. ----- */
void print_edist_matrix(int n, int N_X, double **EDM, int type);

#endif	/* End of __PHYCLUST_EDIST_STRUCT_. */

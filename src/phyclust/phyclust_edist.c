/* This file contains updating functions for edist matrix/structure. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_edist.h"
#include "phyclust_em_tool.h"
#include "phyclust_tool.h"


/* Initial a edist structure by given form and N_X. */
edist_struct* initialize_edist_struct(int form, int N_X){
	edist_struct *eds;

	eds = (edist_struct*) malloc(sizeof(edist_struct));
	eds->form = form;
	eds->N_X = N_X;
	if(form == UT){
		eds->EDM = allocate_s_double_UT(N_X - 1);
		eds->get_pair_edist = &get_pair_edist_UT;
	} else if(form == LT_pam){
		eds->EDM = allocate_s_double_LT_pam(N_X - 1);
		eds->get_pair_edist = &get_pair_edist_LT_pam;
	} else{
		fprintf_stderr("PE: The form of edist structure is not found.\n");
		exit(1);
	}
	return(eds);
} /* End of initialize_edist_struct(). */

void free_edist_struct(edist_struct *eds){
	if(eds->form == LT_pam){
		eds->EDM[0]--;
	}
	free(eds->EDM[0]);
	free(eds->EDM);
	free(eds);
} /* End of free_edist_struct(). */




/* In a upper triangular without diagonal format. */
edist_struct* initialize_edist_struct_UT(int edist_model, int N_X, int L, int **X){
	int i, j, I = N_X - 1;
	edist_struct *eds = initialize_edist_struct(UT, N_X);
	double (*edist_D)(int, int*, int*) = get_edist_D(edist_model);

	for(i = 0; i < I; i++){
		for(j = 0; j < (I - i); j++){
			eds->EDM[i][j] = edist_D(L, X[i], X[i + j + 1]);
		}
	}

	return(eds);
} /* End of initialize_edist_struct_UT(). */

/* In a lower triangular without diagonal format. */
edist_struct* initialize_edist_struct_LT_pam(int edist_model, int N_X, int L, int **X){
	int i, j, I = N_X - 1;
	edist_struct *eds = initialize_edist_struct(LT_pam, N_X);
	double (*edist_D)(int, int*, int*) = get_edist_D(edist_model);

	for(i = 0; i < I; i++){
		for(j = 0; j < i + 1; j++){
			eds->EDM[i][j] = edist_D(L, X[j], X[i + 1]);
		}
	}

	return(eds);
} /* End of initialize_edist_struct_LT_pam(). */




/* Evolution distance.
 * Input: X and Mu with length = L. */
double (*get_edist_D(int edist_model))(int, int*, int*){
	switch(edist_model){
		case D_JC69:
			return(&edist_D_JC69);
			break;
		case D_K80:
			return(&edist_D_K80);
			break;
		case D_HAMMING:
			return(&edist_D_HAMMING);
			break;
		case D_HAMMING_WOGAP:
			return(&edist_D_HAMMING_WOGAP);
			break;
		default:
			fprintf_stderr("PE: Evolution distance model is not found.\n");
			exit(1);
			return(NULL);
			break;
	}
} /* End of (*get_edist_D(int edist_model))(). */

double edist_D_JC69(int L, int *x, int *mu){
	int i, diff = 0, L_NOGAP = L;
	double d;

	for(i = 0; i < L; i++){
		if(x[i] == GAP || mu[i] == GAP ||
			x[i] == MISSING_ALLELE || mu[i] == MISSING_ALLELE){
			L_NOGAP--;
			continue;
		}
		if(x[i] != mu[i]){
			diff++;
		}
	}

	d = 1.0 - 4.0/3.0 * (double) diff / (double) L_NOGAP;
	if(d > 0){
		/*
		d = -0.75 * log(d);
		return((d > 0) ? d : 0.0);
		*/
		return(-0.75 * log(d));
	} else{
		return(Inf);
	}
} /* End of edist_D_JC69(). */

double edist_D_K80(int L, int *x, int *mu){
	int i, L_NOGAP = L;
	double PQ = 0.0, P = 0, Q = 0, d_1, d_2;

	for(i = 0; i < L; i++){
		if(x[i] == GAP || mu[i] == GAP ||
			x[i] == MISSING_ALLELE || mu[i] == MISSING_ALLELE){
			L_NOGAP--;
			continue;
		}
		if(x[i] != mu[i]){
			PQ++;
		}
		if((x[i] == A && mu[i] == G) || (x[i] == G && mu[i] == A) ||
			(x[i] == C && mu[i] == T) || (x[i] == T && mu[i] == C)){
			P++;
		}
	}
	Q = (PQ - P) / (double) L_NOGAP;
	P = P / (double) L_NOGAP;

	d_1 = 1 - 2 * P - Q;
        d_2 = 1 - 2 * Q;
	if(d_1 > 0 && d_2 > 0){
		/*
		d_1 = -0.5 * log(d_1) - 0.25 * log(d_2);
		return((d_1 > 0) ? d_1 : 0.0);
		*/
		return(-0.5 * log(d_1) - 0.25 * log(d_2));
	} else{
		return(Inf);
	}
} /* End of edist_D_K80(). */

double edist_D_HAMMING(int L, int *x, int *mu){
	int i, diff = 0;

	for(i = 0; i < L; i++){
		if(x[i] != mu[i]){
			diff++;
		}
	}

	return((double) diff);
} /* End of edist_D_HAMMING(). */

double edist_D_HAMMING_WOGAP(int L, int *x, int *mu){
	int i, diff = 0;

	for(i = 0; i < L; i++){
		if(x[i] == GAP || mu[i] == GAP ||
			x[i] == MISSING_ALLELE || mu[i] == MISSING_ALLELE){
			continue;
		}
		if(x[i] != mu[i]){
			diff++;
		}
	}

	return((double) diff);
} /* End of edist_D_HAMMING_WOGAP(). */




/* Tools. */
double get_pair_edist_UT(edist_struct *eds, int u, int v){
	if(u > v){
		return(eds->EDM[v][u - v - 1]);
	} else if(u < v){
		return(eds->EDM[u][v - u - 1]);
	} else{
		return(0.0);
	}
} /* get_pair_edist_UT(). */

double get_pair_edist_LT_pam(edist_struct *eds, int u, int v){
	if(u > v){
		return(eds->EDM[u][v]);
	} else if(u < v){
		return(eds->EDM[v][u]);
	} else{
		return(0.0);
	}
} /* get_pair_edist_LT_pam(). */




/* ----- For debug. ----- */
void print_edist_matrix_UT(int first_N_X, int N_X, double **EDM){
	int i, j, I = first_N_X - 1;

	printf("  ");
	for(i = 1; i < first_N_X; i++){
		printf(" %8d", i);
	}
	printf("\n");
	for(i = 0; i < I; i++){
		printf("%2d", i);
		for(j = 0; j < i; j++){
			printf("         ");
		}
		for(j = 0; j < (I - i); j++){
			if(is_finite(EDM[i][j])){
				printf(" %8.4f", EDM[i][j]);
			} else{
				printf(" %8.1e", EDM[i][j]);
			}
		}
		printf("\n");
	}
} /* End of print_edist_matrix_UT(). */

void print_edist_matrix_full(int first_N_X, int N_X, double **EDM){
	int i, j;

	printf("  ");
	for(i = 0; i < first_N_X; i++){
		printf(" %8d", i);
	}
	printf("\n");
	for(i = 0; i < first_N_X; i++){
		printf("%2d", i);
		for(j = 0; j < i; j++){
			if(is_finite(EDM[j][i - j - 1])){
				printf(" %8.4f", EDM[j][i - j - 1]);
			} else{
				printf(" %8.1e", EDM[j][i - j - 1]);
			}
		}
		printf("         ");
		for(j = 0; j < (first_N_X - i - 1); j++){
			if(is_finite(EDM[i][j])){
				printf(" %8.4f", EDM[i][j]);
			} else{
				printf(" %8.1e", EDM[i][j]);
			}
		}
		printf("\n");
	}
} /* End of print_edist_matrix_full(). */

void print_edist_matrix(int first_N_X, int N_X, double **EDM, int type){
	if(first_N_X > N_X){
		first_N_X = N_X;
	}
	switch(type){
		case 0:
			print_edist_matrix_UT(first_N_X, N_X, EDM);
			break;
		case 1:
			print_edist_matrix_full(first_N_X, N_X, EDM);
			break;
		default:
			printf("Printing method is not found.\n");
	};
} /* End of print_dist_matrix(). */


/* This file contains tool fucntions for PAM initialization.
 * Since this file will call pam(), the variables n and N are the same as
 * n_X and N_X in other files. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_init_method.h"
#include "phyclust_tool.h"

void phyclust_pam(int nn, int kk, double *dys, int *nsend, 
		 int/*logical*/ *nrepr, int *nelem, 
		 double *radus, double *damer, double *avsyl,
		 double *ttsyl, double *obj, int *med, int *ncluv, int *nisol);

/* EDM: array(N * (N - 1) / 2), dim = LT.
 * center_id: array(K).
 * class_id: array(N).
 */

void print_pam(int N, int K, int *center_id, int *class_id){
	int n, k;

	printf("  class_id: ");
	for(n = 0; n < N; n++){
		printf("%d ", class_id[n]);
	}
	printf("\n");
	printf(" center_id: ");
	for(k = 0; k < K; k++){
		printf("%d ", center_id[k]);
	}
	printf("\n");
} /* End of print_pam(). */

void assign_class_by_pam(int N_X, int K, double **EDM_LT_pam, int *center_id, int *class_id){
	int i;
	double *dys, *radus, *damer, *avsyl, *ttsyl, *obj;
	int *nsend, *nrepr, *nelem, /* *med, *ncluv, */ *nisol;

	radus = allocate_double_1D(N_X);
	damer = allocate_double_1D(N_X);
	avsyl = allocate_double_1D(N_X);
	ttsyl = allocate_double_1D(1);	/* set ttsyl = 0. */
	obj = allocate_double_1D(2);	/* set obj = c(0, 0). */
	nsend = allocate_int_1D(N_X);
	nrepr = allocate_int_1D(N_X);
	nelem = allocate_int_1D(N_X);
	/* med = center_id; */		/* set med[...] = 0 for non initial, and store center id. */
	/* ncluv = class_id; */		/* store class id. (1,...,kk) */
	nisol = allocate_int_1D(1);	/* set nisol = 1 for swap, 0 o/w. */

	/* dys structure, Kaufman and Rousseeuw (1990), p105.
	 * 0
	 * 9 0		=> dys = 0 9 1 4 3 2 7 ....
	 * 1 4 0		 ^
	 * 3 2 7 0		 The first element is 0 and then by LT form
	 * ...			 due to the f2c code.
	 */
	dys = EDM_LT_pam[0] - 1;
	for(i = 0; i < N_X; i++){
		nsend[i] = nrepr[i] = nelem[i] = class_id[i] /* ncluv[i] */ = 0;
		radus[i] = damer[i] = avsyl[i] = 0;
	}
	for(i = 0; i < K; i++){
		center_id[i] = 0; /* med[i] = 0; */
	}
	ttsyl[0] = 0;
	obj[0] = obj[1] = 0;
	nisol[0] = 1;

	phyclust_pam(N_X, K, dys, nsend, nrepr, nelem, radus, damer, avsyl,
			ttsyl, obj, center_id, class_id, /* med, ncluv, */ nisol);

	for(i = 0; i < N_X; i++){
		class_id[i]--;
	}
	for(i = 0; i < K; i++){
		center_id[i]--;
	}

	free(radus);
	free(damer);
	free(avsyl);
	free(ttsyl);
	free(obj);
	free(nsend);
	free(nrepr);
	free(nelem);
	free(nisol);
} /* End of assign_class_by_pam(). */


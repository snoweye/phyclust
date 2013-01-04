/* This file contains tool fucntions for K-Medoids initialization.
 * Pick centers from unique sequences and evaluate total cost on full data set
 * to avoid degenerated clusters. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_init_method.h"

/* EDM: array(N_X_org * (N_X_org - 1) / 2), dim = UT.
 * center_id: array(K).
 * new_center_id: array(K).
 * total_cost: scalar.
 * new_total_cost: scalar.
 * class_id: array(N_X_org).
 * new_class_id: array(N_X_org). */

void assign_class_id_compute_total_cost(int N_X, int K, double **EDM, int *center_id, int *new_class_id,
		double *new_total_cost){
	int n_X, k;
	double tmp_dist, dist;

	/* For each point, identify cluster and total cost. */
	*new_total_cost = 0.0;
	for(n_X = 0; n_X < N_X; n_X++){
		new_class_id[n_X] = center_id[0];
		if(n_X < center_id[0]){
			dist = EDM[n_X][center_id[0] - n_X - 1];
		} else if(n_X > center_id[0]){
			dist = EDM[center_id[0]][n_X - center_id[0] - 1];
		} else{
			continue;
		}

		for(k = 1; k < K; k++){
			if(n_X < center_id[k]){
				tmp_dist = EDM[n_X][center_id[k] - n_X - 1];
			} else if(n_X > center_id[k]){
				tmp_dist = EDM[center_id[k]][n_X - center_id[k] - 1];
			} else{
				new_class_id[n_X] = center_id[k];
				dist = 0.0;
				break;
			}

			if(tmp_dist < dist){
				new_class_id[n_X] = center_id[k];
				dist = tmp_dist;
			}
		}

		*new_total_cost += dist;
	}
} /* End of assign_class_id_compute_total_cost(). */

void print_kmed(int N_X, int K, int *center_id, int *class_id, double total_cost, int *new_center_id){
	int n_X, k;

	printf("total_cost: %8.4f\n", total_cost);
	printf("  class_id: ");
	for(n_X = 0; n_X < N_X; n_X++){
		printf("%d ", class_id[n_X]);
	}
	printf("\n");
	printf(" center_id: ");
	for(k = 0; k < K; k++){
		printf("%d ", center_id[k]);
	}
	printf("\n");
	printf("new_center: ");
	for(k = 0; k < K; k++){
		printf("%d ", new_center_id[k]);
	}
	printf("\n");
} /* End of print_kmed(). */




void classify_by_EDM(int N_X, int K, double **EDM, int *center_id, int *class_id, int *new_center_id){
	int i, n_X, k;
	double dist, min_cost, class_cost[N_X];

	/* Compute pairwised cost within cluster. */
	for(n_X = 0; n_X < N_X; n_X++){
		class_cost[n_X] = 0.0;
	}
	for(n_X = 0; n_X < (N_X - 1); n_X++){
		for(i = (n_X + 1); i < N_X; i++){
			if(class_id[n_X] == class_id[i]){
				dist = EDM[n_X][i - n_X - 1];
				class_cost[n_X] += dist;
				class_cost[i] += dist;
			}
		}
	}

	/* Assign new centers. */
	for(k = 0; k < K; k++){
		new_center_id[k] = center_id[k];
		min_cost = Inf;
		for(n_X = 0; n_X < N_X; n_X++){
			if(class_id[n_X] == center_id[k]){
				if(class_cost[n_X] < min_cost){
					min_cost = class_cost[n_X];
					new_center_id[k] = n_X;
				}
			}
		}
	}
} /* End of classify_by_EDM(). */

void assign_class_by_k_medoids(int N_X, int K, double **EDM, int *center_id, int *class_id){
	int max_kmed_iter = 100;
	int i, n_X, k;
	int new_center_id[K], new_class_id[N_X];
	double total_cost = 0.0, new_total_cost = 0.0;

	srswor(N_X, K, center_id);
	assign_class_id_compute_total_cost(N_X, K, EDM, center_id, class_id, &total_cost);

	for(i = 0; i < max_kmed_iter; i++){
		classify_by_EDM(N_X, K, EDM, center_id, class_id, new_center_id);
		#if INITDEBUG > 0
			printf("iter = %d\n", i);
			print_kmed(N_X, K, center_id, class_id, total_cost, new_center_id);
		#endif
		assign_class_id_compute_total_cost(N_X, K, EDM, new_center_id, new_class_id, &new_total_cost);

		if(new_total_cost < total_cost){
			total_cost = new_total_cost;
			for(k = 0; k < K; k++){
				center_id[k] = new_center_id[k];
			}
			for(n_X = 0; n_X < N_X; n_X++){
				class_id[n_X] = new_class_id[n_X];
			}
		} else{
			break;
		}
	}

	/* Recordeing class_id. */
	for(n_X = 0; n_X < N_X; n_X++){
		for(k = 0; k < K; k++){
			if(class_id[n_X] == center_id[k]){
				class_id[n_X] = k;
				break;
			}
		}
	}
} /* End of assign_class_by_k_medoids(). */




void classify_unique_by_EDM(int N_X_org, int K, double **EDM, int N_X, int *map_X_to_X_org,
		int *center_id, int *class_id, int *new_center_id){
	int i, n_X_org, n_X, k, tmp_X_org[N_X];
	double dist, min_cost, class_cost[N_X_org];

	/* Compute pairwised cost within cluster. */
	for(n_X = 0; n_X < N_X; n_X++){
		tmp_X_org[n_X] = map_X_to_X_org[n_X];
	}
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		class_cost[n_X_org] = 0.0;
	}
	for(n_X_org = 0; n_X_org < (N_X_org - 1); n_X_org++){
		for(i = (n_X_org + 1); i < N_X_org; i++){
			if(class_id[n_X_org] == class_id[i]){
				dist = EDM[n_X_org][i - n_X_org - 1];
				class_cost[n_X_org] += dist;
				class_cost[i] += dist;
			}
		}
	}

	/* Assign new centers. */
	for(k = 0; k < K; k++){
		new_center_id[k] = center_id[k];
		min_cost = Inf;
		for(n_X = 0; n_X < N_X; n_X++){
			if(class_id[tmp_X_org[n_X]] == center_id[k]){
				if(class_cost[n_X] < min_cost){
					min_cost = class_cost[n_X];
					new_center_id[k] = n_X;
				}
			}
		}
	}
} /* End of classify_unique_by_EDM(). */

void assign_class_unique_by_k_medoids(int N_X_org, int K, double **EDM, int N_X, int *map_X_to_X_org,
		int *center_id, int *class_id){
	int max_kmed_iter = 1000;
	int i, n_X_org, k;
	int new_center_id[K], new_class_id[N_X_org];
	double total_cost = 0.0, new_total_cost = 0.0;

	srswor(N_X, K, center_id);
	for(k = 0; k < K; k++){
		center_id[k] = map_X_to_X_org[center_id[k]];
	}
	assign_class_id_compute_total_cost(N_X_org, K, EDM, center_id, class_id, &total_cost);

	for(i = 0; i < max_kmed_iter; i++){
		classify_unique_by_EDM(N_X_org, K, EDM, N_X, map_X_to_X_org, center_id, class_id,
				new_center_id);
		#if INITDEBUG > 0
			printf("iter = %d\n", i);
			print_kmed(N_X_org, K, center_id, class_id, total_cost, new_center_id);
		#endif
		assign_class_id_compute_total_cost(N_X_org, K, EDM, new_center_id, new_class_id, &new_total_cost);

		if(new_total_cost < total_cost){
			total_cost = new_total_cost;
			for(k = 0; k < K; k++){
				center_id[k] = new_center_id[k];
			}
			for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
				class_id[n_X_org] = new_class_id[n_X_org];
			}
		} else{
			break;
		}
	}

	/* Recordeing class_id. */
	for(n_X_org = 0; n_X_org < N_X_org; n_X_org++){
		for(k = 0; k < K; k++){
			if(class_id[n_X_org] == center_id[k]){
				class_id[n_X_org] = k;
				break;
			}
		}
	}
} /* End of assign_class_unique_by_k_medoids(). */


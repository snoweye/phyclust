/* This file contains tool fucntions for NJ initialization.
 * Since this file will call ape_nj(), the variables n and N are the same as
 * n_X and N_X in other files. */

#include <stdio.h>
#include <stdlib.h>
#include "phyclust_constant.h"
#include "phyclust_init_method.h"
#include "phyclust_tool.h"

int search_njs_edge1(int n, int from, nj_struct *njs, int *anc_id1, int *anc_id2);
int search_njs_edge2(int n, int from, nj_struct *njs, int *anc_id1, int *anc_id2);

int search_njs_edge1(int n, int from, nj_struct *njs, int *anc_id1, int *anc_id2){
	int i, ret = -1;

	if(n == *anc_id1 || n == *anc_id2){
		ret = n;
	} else{
		for(i = 0; i < njs->n_edge; i++){
			if(ret == -1 && njs->edge1[i] == n && njs->edge2[i] != from && njs->edge2[i] > njs->N){
				ret = search_njs_edge1(njs->edge2[i], n, njs, anc_id1, anc_id2);
				if(ret == -1){
					ret = search_njs_edge2(njs->edge2[i], n, njs, anc_id1, anc_id2);
				}
			}
		}
	}

	return(ret);
} /* End of search_njs_edge1(). */

int search_njs_edge2(int n, int from, nj_struct *njs, int *anc_id1, int *anc_id2){
	int i, ret = -1;

	if(n == *anc_id1 || n == *anc_id2){
		ret = n;
	} else{
		for(i = 0; i < njs->n_edge; i++){
			if(ret == -1 && njs->edge2[i] == n && njs->edge1[i] != from){
				ret = search_njs_edge2(njs->edge1[i], n, njs, anc_id1, anc_id2);
				if(ret == -1){
					ret = search_njs_edge1(njs->edge1[i], n, njs, anc_id1, anc_id2);
				}
			}
		}
	}

	return(ret);
} /* End of search_njs_edge2(). */

int search_njs(int n, nj_struct *njs, int *anc_id1, int *anc_id2){
	int i, ret = -1;

	for(i = 0; i < njs->n_edge; i++){
		if(njs->edge2[i] == n){
			ret = search_njs_edge1(njs->edge1[i], n, njs, anc_id1, anc_id2);
			if(ret == -1){
				ret = search_njs_edge2(njs->edge1[i], n, njs, anc_id1, anc_id2);
			}
			break;
		}
	}

	return(ret);
} /* End of search_njs(). */

void search_largest_branch(nj_struct *njs, int *largest_branch_id){
	int i, j, swap_id, tmp_id;

	for(j = 0; j < njs->n_internal_edge; j++){
		largest_branch_id[j] = -1;
	}
	for(i = 0; i < njs->n_edge; i++){
		if(njs->edge2[i] <= njs->N){		/* Escape for leaf nodes. */
			continue;
		} else{
			tmp_id = i;
			for(j = 0; j < njs->n_internal_edge; j++){
				if(largest_branch_id[j] != -1){
					if(njs->edge_length[tmp_id] > njs->edge_length[largest_branch_id[j]]){
						swap_id = largest_branch_id[j];
						largest_branch_id[j] = tmp_id;
						tmp_id = swap_id;
					}
				} else{
					largest_branch_id[j] = tmp_id;
					break;
				}
			}
		}
	}
} /* End of search_largest_branch(). */

/* random_id has length njs->n_internal_edge. */
void random_branch(nj_struct *njs, int *random_branch_id){
	int i, j, tmp_id;

	srswor(njs->n_internal_edge, njs->n_internal_edge, random_branch_id);
	for(i = 0; i < njs->n_internal_edge; i++){
		tmp_id = random_branch_id[i];
		#if INITDEBUG > 1
			printf("i = %d, random_branch_id = %d ", i, tmp_id);
		#endif
		for(j = 0; j < njs->n_edge; j++){
			if(njs->edge2[j] > njs->N){
				tmp_id--;
			}
			if(tmp_id < 0){
				break;
			}
		}
		random_branch_id[i] = j;
		#if INITDEBUG > 1
			printf("edge_id = %d\n", j);
		#endif
	}
} /* End of random_branch(). */

void print_nj_id(int N, int *class_id){
	int n;

	printf("label: ");
	for(n = 0; n < N; n++){
		printf("%2d ", n + 1);
	}
	printf("\n");
	printf("class: ");
	for(n = 0; n < N; n++){
		printf("%2d ", class_id[n]);
	}
	printf("\n");
} /* End of print_nj_id(). */

void print_nj_id_new(int N, int *new_class_id){
	int n;

	printf("  new: ");
	for(n = 0; n < N; n++){
		printf("%2d ", new_class_id[n]);
	}
	printf("\n");
} /* End of print_nj_id_new(). */




/* branch_id: array(njs->n_internal_edge).
 * class_id: array(njs->N).
 * return 0 for successful assignments.
 *        1 for exhaust all internal_edge.
 *        2 for error.
 * */
int assign_class_by_njs_branch(int K, nj_struct *njs, int *branch_id, int *class_id){
	int N = njs->N;
	int i, n, k, ret, new_class_id[N], K_id[K], max_K, flag;

	if(K < 1 || K > njs->n_internal_edge){
		fprintf_stderr("PE: K is out of range (%d, %d)\n", 2, njs->n_internal_edge);
		exit(1);
	}

	for(n = 0; n < N; n++){
		class_id[n] = 0;
	}
	if(K == 1){
		return(0);
	}

	for(i = 0; i < njs->n_internal_edge; i++){
		#if INITDEBUG > 0
			printf("anc_id1 = %d, anc_id2 = %d\n", njs->edge1[branch_id[i]], njs->edge2[branch_id[i]]);
		#endif
		for(n = 0; n < N; n++){
			/* ape_nj() uses 1:N for sequence id. */
			ret = search_njs(n + 1, njs, &njs->edge1[branch_id[i]], &njs->edge2[branch_id[i]]);
			if(ret == njs->edge1[branch_id[i]]){
				new_class_id[n] = 0;
			} else if(ret == njs->edge2[branch_id[i]]){
				new_class_id[n] = 1;
			} else{
				fprintf_stderr("PE: Sequence (%d) is not assigned.\n", n + 1);
				new_class_id[n] = -1;
				return(2);
			}
		}
		#if INITDEBUG > 0
			print_nj_id(N, class_id);
			print_nj_id_new(N, new_class_id);
		#endif

		max_K = 0;
		for(k = 0; k < K; k++){
			K_id[k] = -1;
		}
		for(n = 0; n < N; n++){
			class_id[n] = (class_id[n] << 1) | new_class_id[n];
			flag = 0;
			for(k = 0; k < K; k++){
				if(K_id[k] == -1){
					K_id[k] = class_id[n];
					max_K++;
					flag = 1;
					break;
				} else if(K_id[k] == class_id[n]){
					flag = 1;
					break;
				}
			}
			if(flag == 0 && k == K){
				fprintf_stderr("PE: Classes are over assigned.\n");
				return(3);
			}
		}
		#if INITDEBUG > 0
			printf("Update:\n");
			print_nj_id(N, class_id);
		#endif

		if(max_K == K){
			for(n = 0; n < N; n++){
				for(k = 0; k < K; k++){
					if(class_id[n] == K_id[k]){
						class_id[n] = k;
						break;
					}
				}
			}
			#if INITDEBUG > 0
				printf("Final:\n");
				print_nj_id(N, class_id);
			#endif
			return(0);
		}
	}

	return(1);
} /* End of assign_class_by_njs_branch(). */


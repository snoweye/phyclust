/* This file contains all functions required in em steps .*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "phyclust_constant.h"
#include "phyclust_label.h"
#include "phyclust_tool.h"


/* Initial a phyclust label. */
phyclust_label* initialize_phyclust_label(){
	phyclust_label *pcl;

	pcl = (phyclust_label*) malloc(sizeof(phyclust_label));
	pcl->label_method = NONE; 

	pcl->N_index_org = 0;		/* Assigned by update_phyclust_label(). */
	pcl->N_index = 0;		/* Assigned by update_phyclust_label(). */

	pcl->semi_org = NULL;		/* Assigned by update_phyclust_label(). */
	pcl->semi = NULL;		/* Assigned by update_phyclust_label(). */

	pcl->index_org = NULL;		/* Assigned by update_phyclust_label(). */
	pcl->index = NULL;		/* Assigned by update_phyclust_label(). */

	pcl->prob_org = NULL;		/* Assigned by update_phyclust_label(). */
	pcl->prob = NULL;		/* Assigned by update_phyclust_label(). */

	return(pcl);
} /* End of initialize_phyclust_label(). */

void free_phyclust_label(phyclust_label *pcl){
	if(pcl->label_method == SEMI){
		free(pcl->semi);
		free(pcl->index);
		free(pcl->prob_org);
		free(pcl->prob);
	} else if(pcl->label_method == GENERAL){
		free(pcl->index);
		free(pcl->prob_org);
		free(pcl->prob);
	}
	free(pcl);
} /* End of free_phyclust_label(). */

/* tmp_prob should be a read-only object. */
void update_phyclust_label(phyclust_label *pcl, int label_method, int N_label, int *label_semi,
		int *label_index, double *tmp_prob, int *map_X_org_to_X, int K){
	int i, i2, n_index, flag, tmp_i;

	if(N_label > 0){
		pcl->label_method = label_method;
		pcl->N_index_org = N_label;
		pcl->N_index = 0;

		/* Only the information of the first labeled sequence
		 * among all other identical sequences will be used. */
		n_index = 0;
		for(i = 0; i < pcl->N_index_org; i++){
			flag = 0;
			tmp_i = map_X_org_to_X[label_index[i]];
			for(i2 = 0; i2 < i; i2++){
				if(tmp_i == map_X_org_to_X[label_index[i2]]){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				n_index++;
			}
		}
		pcl->N_index = n_index;

		if(label_method == SEMI){
			pcl->semi_org = label_semi;
			pcl->semi = allocate_int_1D(pcl->N_index);
		}

		pcl->index_org = label_index;
		pcl->index = allocate_int_1D(pcl->N_index);

		pcl->prob_org = allocate_double_2D_AP(pcl->N_index_org);
		pcl->prob = allocate_double_2D_AP(pcl->N_index);

		n_index = 0;
		for(i = 0; i < pcl->N_index_org; i++){
			/* Assign pointers to obtain data from R. */
			pcl->prob_org[i] = tmp_prob;
			tmp_prob += K;

			if(n_index < pcl->N_index){
				flag = 0;
				tmp_i = map_X_org_to_X[label_index[i]];
				for(i2 = 0; i2 < i; i2++){
					if(tmp_i == map_X_org_to_X[label_index[i2]]){
						flag = 1;
						break;
					}
				}
				if(flag == 0){
					if(label_method == SEMI){
						pcl->semi[n_index] = label_semi[i];
					}
					pcl->index[n_index] = tmp_i;
					pcl->prob[n_index] = pcl->prob_org[i];
					n_index++;
				}
			}
		}
	}
} /* End of update_phyclust_label(). */


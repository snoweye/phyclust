/* This file contains declarations for semi-supervised clustering. */


#ifndef __PHYCLUST_LABEL_
#define __PHYCLUST_LABEL_

typedef struct _phyclust_label		phyclust_label;


/* Control label and index.
 * This is a read-only object and a reference for labeling.
 *
 * For computing, use empcs to indicate labeled and unlabed
 * for N_X, **Z_modified, and **Z_normalized.
 *
 * pcl should be initialized nomatter for un- or semi-supervised
 * clustering, point to NULL if no labels (unsupervised).
 *
 * Now, we only support a semi labeling method. For more general labeling,
 * this will be implemented later for general semi-supervised clustering,
 * since this requires more complicated computation of likelihood by
 * solving systems of linear equations. */
struct _phyclust_label{
	int		label_method;		/* Either NONE, SEMI, or GENERAL. */

	/* Storage. */
	int		N_index_org;		/* Length of index_org. */
	int		N_index;		/* Length of index. */

	/* Simple index, take values from 0, ..., (K_label - 1)
	 * These arraies indicate which cluster the sequence belongs to.
	 * -1 should be removed from label$semi.
	 * This summarizes label$semi corresponding label$index. */
	int		*semi_org;		/* original label, length = N_index_org. */
	int		*semi;			/* Unique label, length = N_index. */

	/* Index to the sequences. These arraies contain labeled seqence index.
	 * index_org takes value from 0 to (N_X_org - 1).
	 * index takes value from 0 to (N_X - 1). */
	int		*index_org;		/* Original index, length = N_index_org. */
	int		*index;			/* Unique index, length = N_index. */

	/* Storage for probability.
	 * This will be copied to empcs->Z_normalized_labeled. */
	double		**prob_org;		/* Original labeled probability, dim = N_index_org * K. */
	double		**prob;			/* Unique labeled probability, repoint to prob_org, length = N_index. */
};

phyclust_label* initialize_phyclust_label(void);
void free_phyclust_label(phyclust_label *pcl);
void update_phyclust_label(phyclust_label *pcl, int label_method, int N_label, int *label_semi,
	int *label_index, double *tmp_prob, int *map_X_org_to_X, int K);

#endif	/* End of __PHYCLUST_LABEL. */


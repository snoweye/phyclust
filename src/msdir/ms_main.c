/* This file contains a wrap for R to call ms_main() in "ms.c". */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ms.h"

/* This function modified from main() in "ms.c". */
void ms_main(int argc, char *argv[]){
	int i, k, howmany, segsites ;
	char **list, **cmatrix(), **tbsparamstrs ;
	double probss, tmrca, ttot ;

	R_ms_file_pointer = fopen(R_ms_file_name, "w");

	ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
	tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;

	for( i=0; i<argc; i++) fprintf(R_ms_file_pointer, "%s ",argv[i]);
	for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
	for( i = 1; i<argc ; i++)
			if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;
	
	count=0;

	if( ntbs > 0 )  for( k=0; k<ntbs; k++)  i = scanf(" %s", tbsparamstrs[k] );
	getpars( argc, argv, &howmany) ;   // results are stored in global variable, pars
	
	if( pars.mp.segsitesin ==  0 ) {
	     list = cmatrix(pars.cp.nsam,maxsites+1);
	     posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
	}
	else {
	     list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
	     posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
	     if( pars.mp.theta > 0.0 ){
		    segfac = 1.0 ;
		    for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
		 }
	}

	while(howmany - count++) {
		if((ntbs > 0) && (count > 1)){
			for(k = 0; k < ntbs; k++){ 
				if(scanf(" %s", tbsparamstrs[k]) == EOF){
					exit(0);
				}
			}
			free_pars();
			getpars(argc, argv, &howmany);
		}
	   
		fprintf(R_ms_file_pointer, "\n//");
		if(ntbs > 0){
			for(k = 0; k < ntbs; k++){
			       fprintf(R_ms_file_pointer, "\t%s", tbsparamstrs[k]);
			}
		}
		fprintf(R_ms_file_pointer, "\n");
		segsites = gensam(list, &probss, &tmrca, &ttot); 
		if(pars.mp.timeflag){
		       fprintf(R_ms_file_pointer, "time:\t%15.10lf\t%15.10lf\n", tmrca, ttot);
		}
		if((segsites > 0) || (pars.mp.theta > 0.0)){
			if((pars.mp.segsitesin > 0) && (pars.mp.theta > 0.0)){
				fprintf(R_ms_file_pointer, "prob: %g\n", probss);
			}
			fprintf(R_ms_file_pointer, "segsites: %d\n", segsites);
			if(segsites > 0){
				fprintf(R_ms_file_pointer, "positions: ");
			}
			for(i = 0; i < segsites; i++){
				fprintf(R_ms_file_pointer, "%15.10lf ", posit[i] );
			}
			fprintf(R_ms_file_pointer, "\n");
			if(segsites > 0){
				for(i = 0; i < pars.cp.nsam; i++){
				       fprintf(R_ms_file_pointer, "%s\n", list[i]);
				}
			}
		}
	}

	free(posit);
	free_char_2D_AP(tbsparamstrs, argc);
	free_char_2D_AP(list, pars.cp.nsam);
	free_pars();
	fclose(R_ms_file_pointer);
} /* End of ms_main(). */


void free_char_2D_AP(char **char_2D_AP, int nrow){
	int i;

	for(i = 0; i < nrow; i++){
		free(char_2D_AP[i]);
	}
	free(char_2D_AP);
} /* End of free_char_2D_AP(). */

void free_pars(){
	int i;
	
	if(pars.cp.config != NULL){
		free(pars.cp.config);
	}
	if(pars.cp.mig_mat != NULL){
		for(i = 0; i < pars.cp.npop; i++){
			free(pars.cp.mig_mat[i]);
		}
		free(pars.cp.mig_mat);
	}
	if(pars.cp.size != NULL){
		free(pars.cp.size);
	}
	if(pars.cp.alphag != NULL){
		free(pars.cp.alphag);
	}
	if(pars.cp.deventlist != NULL){
		free(pars.cp.deventlist);
	}
} /* End of free_pars(). */


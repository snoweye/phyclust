/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.   
        usage: ms nsam howmany  -t  theta  [options]
		or
	       ms nsam howmany -s segsites  [options] 

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary 
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be 
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate 
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The 
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No. 
           howmany: howmany samples to generate.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation. 
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program. 
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
	Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
	Changed many command line options, allowing arbitrary migration matrix, and subpopulation
	   sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
	   split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
	!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
	Fixed bug which resulted in incorrect results for the case where
             rho = 0.0 and gene conversion rate > 0.0. This case was not handled
	    correctly in early versions of the program. 5 Apr 2004.  (Thanks to
	    Vincent Plagnol for pointing out this problem.) 
	Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
		 1 Jul 2004.
	Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
	Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
	Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
	Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
	    to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
	Added -seeds option 4 Nov 2006
	Added "tbs" arguments feature 4 Nov 2006
	Added -L option.  10 May 2007
	Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
	fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
          machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
          Oct. 17, 2007.
***************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"

//WCC #define SITESINC 10 

unsigned maxsites = SITESINC ;

/* WCC move to "ms.h".
struct node{
	int abv;
	int ndes;
	float time;
	};

struct segl {
	int beg;
	struct node *ptree;
	int next;
	};
*/

double *posit ;
double segfac ;
int count, ntbs, nseeds ;
struct params pars ;	

/* WCC
main(argc,argv)
        int argc;
        char *argv[];
{
	int i, k, howmany, segsites ;
	char **list, **cmatrix(), **tbsparamstrs ;
	FILE *pf, *fopen() ;
	double probss, tmrca, ttot ;
	void seedit( const char * ) ;
	void getpars( int argc, char *argv[], int *howmany )  ;
	int gensam( char **list, double *probss, double *ptmrca, double *pttot ) ;


	ntbs = 0 ;   // these next few lines are for reading in parameters from a file (for each sample)
	tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;

	for( i=0; i<argc; i++) printf("%s ",argv[i]);
	for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
	for( i = 1; i<argc ; i++)
			if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;
	
	count=0;

	if( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] );
	getpars( argc, argv, &howmany) ;   // results are stored in global variable, pars
	
	if( !pars.commandlineseedflag ) seedit( "s");
	pf = stdout ;

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

    while( howmany-count++ ) {
	   if( (ntbs > 0) && (count >1 ) ){
	         for( k=0; k<ntbs; k++){ 
			    if( scanf(" %s", tbsparamstrs[k]) == EOF ){
			       if( !pars.commandlineseedflag ) seedit( "end" );
				   exit(0);
				}
			 }
			 getpars( argc, argv, &howmany) ;
	   }
	   
       fprintf(pf,"\n//");
	   if( ntbs >0 ){
			for(k=0; k< ntbs; k++) printf("\t%s", tbsparamstrs[k] ) ;
		}
		printf("\n");
        segsites = gensam( list, &probss, &tmrca, &ttot ) ; 
  		if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
        if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
   	       if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) 
		       fprintf(pf,"prob: %g\n", probss ) ;
           fprintf(pf,"segsites: %d\n",segsites);
		   if( segsites > 0 )	fprintf(pf,"positions: ");
		   for( i=0; i<segsites; i++)
              fprintf(pf,"%6.4lf ",posit[i] );
           fprintf(pf,"\n");
	       if( segsites > 0 ) 
	          for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
	    }
    }
	if( !pars.commandlineseedflag ) seedit( "end" );

}
*/



	int 
gensam( char **list, double *pprobss, double *ptmrca, double *pttot ) 
{
//WCC	int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
	int nsegs, i, k, seg, ns = 0, start, end, len, segsit ;
//WCC	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ; /* used to be: [MAXSEG];  */
	struct segl *seglst; /* used to be: [MAXSEG];  */
//WCC	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
	double nsinv,  tseg, tt;
	double *pk;
	int *ss;
	int segsitesin,nsites;
	double theta, es ;
	int nsam, mfreq ;
//WCC	void prtree( struct node *ptree, int nsam);
//WCC	int make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
//WCC	void ndes_setup( struct node *, int nsam );


	nsites = pars.cp.nsites ;
//WCC	nsinv = 1./nsites;
	nsinv = 1./(double)nsites;
	seglst = segtre_mig(&(pars.cp),  &nsegs ) ;
	
	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin ;
	theta = pars.mp.theta ;
	mfreq = pars.mp.mfreq ;

	if( pars.mp.treeflag ) {
	  	ns = 0 ;
	    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
	      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
		     end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		     start = seglst[seg].beg ;
		     len = end - start + 1 ;
//WCC		     fprintf(stdout,"[%d]", len);
		     fprintf(R_ms_file_pointer, "[%d]", len);
	      }
	      prtree( seglst[seg].ptree, nsam ) ;
	      if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) ) 
	  	      free(seglst[seg].ptree) ;
	    }
	}

	if( pars.mp.timeflag ) {
      tt = 0.0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
		  *ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
		if( mfreq == 1 ) tt += ttime(seglst[seg].ptree,nsam)*tseg ;
		else tt += ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
		if( (segsitesin == 0) && ( theta == 0.0 )  ) 
	  	      free(seglst[seg].ptree) ;
	    }
		*pttot = tt ;
	 }	
	
    if( (segsitesin == 0) && ( theta > 0.0)   ) {
	  ns = 0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len*(theta/nsites) ;
		if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
                else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
		segsit = poisso( tseg*tt );
		if( (segsit + ns) >= maxsites ) {
			maxsites = segsit + ns + SITESINC ;
			posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
			  biggerlist(nsam, list) ; 
		}
		make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
		free(seglst[seg].ptree) ;
		locate(segsit,start*nsinv, len*nsinv,posit+ns);   
		ns += segsit;
	  }
    }
   else if( segsitesin > 0 ) {

        pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
        if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


	  tt = 0.0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
               if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
               else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
                 tt += pk[k] ;
	  }
	  if( theta > 0.0 ) { 
	    es = theta * tt ;
	    *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
	  }
	  if( tt > 0.0 ) {
          for (k=0;k<nsegs;k++) pk[k] /= tt ;
          mnmial(segsitesin,nsegs,pk,ss);
	  }
	  else
            for( k=0; k<nsegs; k++) ss[k] = 0 ;
	  ns = 0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
		 end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		 start = seglst[seg].beg ;
		 len = end - start + 1 ;
		 tseg = len/(double)nsites;
		 make_gametes(nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);

		 free(seglst[seg].ptree) ;
		 locate(ss[k],start*nsinv, len*nsinv,posit+ns);   
		 ns += ss[k] ;
	  }
	  free(pk);
	  free(ss);

    }

	/* Add by Wei-Chen Chen */
	free(seglst);

	for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
	return( ns ) ;
}

	void 
ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

//WCC	int
	void
biggerlist(int nsam,  char **list )
//WCC	int nsam ;
//WCC	char ** list ;
{
	int i;

/*  fprintf(stderr,"maxsites: %d\n",maxsites);  */	
	for( i=0; i<nsam; i++){
	   list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
	   if( list[i] == NULL ) perror( "realloc error. bigger");
	   }
}
	   


/* allocates space for gametes (character strings) */
	char **
cmatrix(int nsam, int len)
//WCC	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}



//WCC	int
	void
locate(int n, double beg, double len, double *ptr)
//WCC	int n;
//WCC	double beg, len, *ptr;
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}

int NSEEDS = 3 ;

  void
getpars(int argc, char *argv[], int *phowmany )
{
	int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
	double migr, mij, psize, palpha ;
//WCC	void addtoelist( struct devent *pt, struct devent *elist ); 
//WCC	void argcheck( int arg, int argc, char ** ) ;
//WCC	int commandlineseed( char ** ) ;
//WCC	void free_eventlist( struct devent *pt, int npop );
	struct devent *ptemp , *pt ;
	FILE *pf ;
	char ch3 ;
	

  if( count == 0 ) {
//WCC	if( argc < 4 ){ fprintf(stderr,"Too few command line arguments\n"); usage();}
	if( argc < 4 ){ REprintf("Too few command line arguments\n"); usage();}
	pars.cp.nsam = atoi( argv[1] );
//WCC	if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); usage();}
	if( pars.cp.nsam <= 0 ) { REprintf("First argument error. nsam <= 0. \n"); usage();}
	*phowmany = atoi( argv[2] );
//WCC	if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); usage();}
	if( *phowmany  <= 0 ) { REprintf("Second argument error. howmany <= 0. \n"); usage();}
	pars.commandlineseedflag = 0 ;
	pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
	pars.cp.track_len = 0. ;
	pars.cp.npop = npop = 1 ;
	pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
	pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
	pars.cp.mig_mat[0][0] =  0.0 ;
	pars.mp.segsitesin = 0 ;
	pars.mp.treeflag = 0 ;
 	pars.mp.timeflag = 0 ;
       pars.mp.mfreq = 1 ;
	pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
	(pars.cp.config)[0] = pars.cp.nsam ;
	pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
	(pars.cp.size)[0] = 1.0  ;
	pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
	(pars.cp.alphag)[0] = 0.0  ;
	pars.cp.nsites = 2 ;
  }
  else{
	npop = pars.cp.npop ;
	free_eventlist( pars.cp.deventlist, npop );
  }
  	pars.cp.deventlist = NULL ;

	arg = 3 ;

	while( arg < argc ){
//WCC		if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage();}
		if( argv[arg][0] != '-' ) { REprintf(" argument should be -%s ?\n", argv[arg]); usage();}
		switch ( argv[arg][1] ){
			case 'f' :
//WCC				if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
				if( ntbs > 0 ) { REprintf(" can't use tbs args and -f option.\n"); exit(1); }
				arg++;
				argcheck( arg, argc, argv);
				pf = fopen( argv[arg], "r" ) ;
//WCC				if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); exit(0);}
				if( pf == NULL ) {REprintf(" no parameter file %s\n", argv[arg] ); exit(0);}
				arg++;
				argc++ ;
				argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
				argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
				argstart = arg ;
				while( fscanf(pf," %s", argv[arg]) != EOF ) {
					arg++;
					argc++;
					argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
				        argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
					}
				fclose(pf);
				argc--;
				arg = argstart ;
				break;
			case 'r' : 
				arg++;
				argcheck( arg, argc, argv);
				pars.cp.r = atof(  argv[arg++] );
				argcheck( arg, argc, argv);
				pars.cp.nsites = atoi( argv[arg++]);
				if( pars.cp.nsites <2 ){
//WCC					fprintf(stderr,"with -r option must specify both rec_rate and nsites>1\n");
					REprintf("with -r option must specify both rec_rate and nsites>1\n");
					usage();
					}
				break;		
			case 'c' : 
				arg++;
				argcheck( arg, argc, argv);
				pars.cp.f = atof(  argv[arg++] );
				argcheck( arg, argc, argv);
				pars.cp.track_len = atof( argv[arg++]);
				if( pars.cp.track_len <1. ){
//WCC					fprintf(stderr,"with -c option must specify both f and track_len>0\n");
					REprintf("with -c option must specify both f and track_len>0\n");
					usage();
					}
				break;		
			case 't' : 
				arg++;
				argcheck( arg, argc, argv);
				pars.mp.theta = atof(  argv[arg++] );
				break;
			case 's' : 
				arg++;
				argcheck( arg, argc, argv);
				if( argv[arg-1][2] == 'e' ){  /* command line seeds */
//WCC					pars.commandlineseedflag = 1 ;
//WCC					if( count == 0 ) nseeds = commandlineseed(argv+arg );
//WCC					arg += nseeds ;
				}
				else {
				    pars.mp.segsitesin = atoi(  argv[arg++] );
				}
				break;
			case 'F' : 
				arg++;
				argcheck( arg, argc, argv);
				pars.mp.mfreq = atoi(  argv[arg++] );
                                if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
//WCC                                    fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
                                    REprintf(" mfreq must be >= 2 and <= nsam/2.\n");
                                    usage();
                                    }
				break;
			case 'T' : 
				pars.mp.treeflag = 1 ;
				arg++;
				break;
			case 'L' : 
				pars.mp.timeflag = 1 ;
				arg++;
				break;
			case 'I' : 
			    arg++;
			    if( count == 0 ) {
				argcheck( arg, argc, argv);
			       	pars.cp.npop = atoi( argv[arg]);
			        pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
				npop = pars.cp.npop ;
				}
			    arg++;
			    for( i=0; i< pars.cp.npop; i++) {
				argcheck( arg, argc, argv);
				pars.cp.config[i] = atoi( argv[arg++]);
				}
			    if( count == 0 ){
				pars.cp.mig_mat = 
                                        (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
				pars.cp.mig_mat[0] = 
                                         (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
				for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] = 
                                         (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
				pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
				pars.cp.alphag = 
                                          (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
			        for( i=1; i< pars.cp.npop ; i++) {
				   (pars.cp.size)[i] = (pars.cp.size)[0]  ;
				   (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
				   }
			        }
			     if( (arg <argc) && ( argv[arg][0] != '-' ) ) {
				argcheck( arg, argc, argv);
				migr = atof(  argv[arg++] );
				}
			     else migr = 0.0 ;
			     for( i=0; i<pars.cp.npop; i++) 
				    for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
			     for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
			     break;
			case 'm' :
//WCC			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			     if( npop < 2 ) { REprintf("Must use -I option first.\n"); usage();}
			     if( argv[arg][2] == 'a' ) {
				    arg++;
				    for( pop = 0; pop <npop; pop++)
				      for( pop2 = 0; pop2 <npop; pop2++){
					     argcheck( arg, argc, argv);
					     pars.cp.mig_mat[pop][pop2]= atof( argv[arg++] ) ;
					  }
				    for( pop = 0; pop < npop; pop++) {
					  pars.cp.mig_mat[pop][pop] = 0.0 ;
					  for( pop2 = 0; pop2 < npop; pop2++){
					    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
					  }
				    }	
				}
			    else {
		             arg++;
			         argcheck( arg, argc, argv);
		             i = atoi( argv[arg++] ) -1;
			         argcheck( arg, argc, argv);
		             j = atoi( argv[arg++] ) -1;
			         argcheck( arg, argc, argv);
		             mij = atof( argv[arg++] );
		             pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
		             pars.cp.mig_mat[i][j] = mij;
			    }
				break;
			case 'n' :
//WCC			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			     if( npop < 2 ) { REprintf("Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, argc, argv);
			    pop = atoi( argv[arg++] ) -1;
			    argcheck( arg, argc, argv);
			    psize = atof( argv[arg++] );
			    pars.cp.size[pop] = psize ;
			   break;
			case 'g' :
//WCC			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			     if( npop < 2 ) { REprintf("Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, argc, argv);
			    pop = atoi( argv[arg++] ) -1;
//WCC			    if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    if( arg >= argc ) { REprintf("Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( argv[arg++] );
			    pars.cp.alphag[pop] = palpha ;
			   break;
			case 'G' :
			    arg++;
//WCC			    if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    if( arg >= argc ) { REprintf("Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( argv[arg++] );
			    for( i=0; i<pars.cp.npop; i++) 
			       pars.cp.alphag[i] = palpha ;
			   break;
			case 'e' :
			    pt = (struct devent *)malloc( sizeof( struct devent) ) ;
			    pt->detype = argv[arg][2] ;
			    ch3 = argv[arg][3] ;
			    arg++;
			    argcheck( arg, argc, argv);
			    pt->time = atof( argv[arg++] ) ;
			    pt->nextde = NULL ;
			    if( pars.cp.deventlist == NULL ) 
				    pars.cp.deventlist = pt ;
			    else if ( pt->time < pars.cp.deventlist->time ) { 
				    ptemp = pars.cp.deventlist ;
				    pars.cp.deventlist = pt ;
				    pt->nextde = ptemp ;	
				}	
			    else
				   addtoelist( pt, pars.cp.deventlist ) ;
			    switch( pt->detype ) {
				case 'N' :
			          argcheck( arg, argc, argv);
				      pt->paramv = atof( argv[arg++] ) ;
				      break;
				case 'G' :
//WCC				  if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
				  if( arg >= argc ) { REprintf("Not enough arg's after -eG.\n"); usage(); }
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'M' :
				    argcheck( arg, argc, argv);
				    pt->paramv = atof( argv[arg++] ) ;
				    break;
				case 'n' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'g' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
//WCC				  if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
				  if( arg >= argc ) { REprintf("Not enough arg's after -eg.\n"); usage(); }
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 's' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->paramv = atof( argv[arg++] ) ;
				  break;
				case 'm' :
				  if( ch3 == 'a' ) {
				     pt->detype = 'a' ;
				     argcheck( arg, argc, argv);
				     npop2 = atoi( argv[arg++] ) ;
				     pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
				     for( pop =0; pop <npop2; pop++){
					   (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
					   for( i=0; i<npop2; i++){
					     if( i == pop ) arg++;
					     else {
				               argcheck( arg, argc, argv); 
					       (pt->mat)[pop][i] = atof( argv[arg++] ) ;
					     }
					   }
				     }
				     for( pop = 0; pop < npop2; pop++) {
					    (pt->mat)[pop][pop] = 0.0 ;
					    for( pop2 = 0; pop2 < npop2; pop2++){
					       if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
					    }
				     }	
				  }
				  else {
			            argcheck( arg, argc, argv);
				        pt->popi = atoi( argv[arg++] ) -1 ;
			            argcheck( arg, argc, argv);
				        pt->popj = atoi( argv[arg++] ) -1 ;
			            argcheck( arg, argc, argv);
				        pt->paramv = atof( argv[arg++] ) ;
				  }
				  break;
				case 'j' :
			          argcheck( arg, argc, argv);
				  pt->popi = atoi( argv[arg++] ) -1 ;
			          argcheck( arg, argc, argv);
				  pt->popj = atoi( argv[arg++] ) -1 ;
				  break;
//WCC				default: fprintf(stderr,"e event\n");  usage();
				default: REprintf("e event\n");  usage();
			    }
			 break;
//WCC			default: fprintf(stderr," option default\n");  usage() ;
			default: REprintf(" option default\n");  usage() ;
			}
		}
		if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) ) {
//WCC			fprintf(stderr," either -s or -t or -T option must be used. \n");
			REprintf(" either -s or -t or -T option must be used. \n");
			usage();
			exit(1);
			}
		sum = 0 ;
		for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
		if( sum != pars.cp.nsam ) {
//WCC			fprintf(stderr," sum sample sizes != nsam\n");
			REprintf(" sum sample sizes != nsam\n");
			usage();
			exit(1);
			}
}


	void
argcheck( int arg, int argc, char *argv[] )
{
	if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
//WCC	   fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
	   REprintf("not enough arguments after %s\n", argv[arg-1] ) ;
//WCC	   fprintf(stderr,"For usage type: ms<return>\n");
	   REprintf("For usage type: ms<return>\n");
	   exit(0);
	  }
}
	
/* This function is rewrote for R.
	int
usage()
{
fprintf(stderr,"usage: ms nsam howmany \n");
fprintf(stderr,"  Options: \n"); 
fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
fprintf(stderr,"\t -T          (Output gene tree.)\n");
fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

exit(1);
} */

//WCC	int
//WCC usage()
	void
usage(void)
{
printf("usage: ms nsam howmany \n");
printf("  Options: \n"); 
printf("   -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
printf("   -s segsites   ( fixed number of segregating sites)\n");
printf("   -T          (Output gene tree.)\n");
printf("   -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
printf("   -r rho nsites     (rho here is 4Nc)\n");
printf("     -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
printf("        if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
printf("   -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
printf("   -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
printf("     -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
printf("     -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
printf("     -n i size_i   (popi has size set to size_i*N0 \n");
printf("     -g i alpha_i  (If used must appear after -M option.)\n"); 
printf("   The following options modify parameters at the time 't' specified as the first argument:\n");
printf("   -eG t alpha  (Modify growth rate of all pop's.)\n");     
printf("   -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
printf("   -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
printf("   -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
printf("   -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
printf("   -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
printf("   -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
printf("   -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
printf("       proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
printf("       Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
printf("   -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
printf("       size, alpha and M are unchanged.\n");  
//WCC fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
//WCC fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");
exit(1);
}


	void
addtoelist( struct devent *pt, struct devent *elist ) 
{
//WCC	struct devent *plast, *pevent, *ptemp  ;
	struct devent *plast = NULL, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}

	void 
free_eventlist( struct devent *pt, int npop )
{
   struct devent *next ;
   int pop ;
   
   while( pt != NULL){
	  next = pt->nextde ;
	  if( pt->detype == 'a' ) {
	     for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
		 free( pt->mat );
	  }
	  free(pt);
	  pt = next ;
   }
}

	
/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/

#define STATE1 '1'
#define STATE2 '0'

//WCC	int
	void
make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
	int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt), 
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);
		for( tip=0; tip < nsam ; tip++) {
		   if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
		   else list[tip][j] = STATE2 ;
		   }
		}
}


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

	double
ttime(struct node *ptree, int nsam)
//WCC	struct node *ptree;
//WCC	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}


	double
ttimemf(struct node *ptree, int nsam, int mfreq)
//WCC	struct node *ptree;
//WCC	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}


	void
prtree(struct node *ptree, int nsam)
//WCC	struct node *ptree;
//WCC	int nsam;
{
//WCC	double t;
	int i, *descl, *descr ;
//WCC	void parens( struct node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

   if( descl[noden] == -1 ) {
//WCC	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	fprintf(R_ms_file_pointer, "s%d:%15.12lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
   else{
//WCC	printf("(");
	fprintf(R_ms_file_pointer, "(");
	parens( ptree, descl,descr, descl[noden] ) ;
//WCC	printf(",");
	fprintf(R_ms_file_pointer, ",");
	parens(ptree, descl, descr, descr[noden] ) ;
//WCC	if( (ptree+noden)->abv == 0 ) printf(");\n"); 
	if( (ptree+noden)->abv == 0 ) fprintf(R_ms_file_pointer, ");\n"); 
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
//WCC	  printf("):%5.3lf", time );
	  fprintf(R_ms_file_pointer, "):%15.12lf", time );
	  }
        }
}

/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/

	int
pickb(int nsam, struct node *ptree, double tt)
//WCC	int nsam;
//WCC	struct node *ptree;
//WCC	double tt;
{
//WCC	double x, y, ran1();
	double x, y, ran1(void);
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		if( y >= x ) return( i ) ;
		}
	return( i );
}

	int
pickbmf(int nsam, int mfreq, struct node *ptree, double tt )
//WCC	int nsam, mfreq;
//WCC	struct node *ptree;
//WCC	double tt;
{
//WCC	double x, y, ran1();
	double x, y, ran1(void);
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	  if( y >= x ) return( i ) ;
	}
	return( i );
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

	int
tdesn(struct node *ptree, int tip, int node )
//WCC	struct node *ptree;
//WCC	int tip, node;
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


/* pick2()  */

	int
pick2(int n, int *i, int *j)
//WCC	int n, *i, *j;
{
//WCC	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

//WCC	int
	void
ordran(int n, double pbuf[])
//WCC	int n;
//WCC	double pbuf[];
{
	ranvec(n,pbuf);
	order(n,pbuf);
//WCC	return;
}


	int
mnmial(int n, int nclass, double p[], int rv[])
//WCC	int n, nclass, rv[];
//WCC	double p[];
{
//WCC	double ran1();
	double ran1(void);
	double x, s;
//WCC	int i, j;
	int i, j = 0;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return(j);
}

//WCC        int
	void
order(int n, double pbuf[])
//WCC        int n;
//WCC        double pbuf[];
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
}


//WCC	int
	void
ranvec(int n, double pbuf[])
//WCC	int n;
//WCC	double pbuf[];
{
	int i;
//WCC	double ran1();
	double ran1(void);

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

//WCC	return;
}



	int
poisso(double u)
//WCC	double u;
{
//WCC	double  cump, ru, ran1(), p, gasdev() ;
	double  cump, ru, ran1(void), p, gasdev(double, double) ;
	int i=1;

	if( u > 30. ) return( (int)(0.5 + gasdev(u,u)) );
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */

double gasdev(double m, double v)
//WCC	double m, v;
{
//WCC	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
//WCC	double ran1();
	double ran1(void);

//WCC	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
//WCC		iset=1;
		return( m + sqrt(v)*v2*fac);
//WCC	} else {
//WCC		iset=0;
//WCC		return( m + sqrt(v)*gset ) ;
//WCC	}
}


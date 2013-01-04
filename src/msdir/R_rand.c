/* This file contains a function, ran1(), to replace the original random
 * number generation in "ms", instead of R function runif().
 * The seed should be preset and reset in R.
 * Wrote: Wei-Chen Chen 2009-10-03. */
#ifdef __HAVE_R_ 
	#include <R.h>
	#include <Rmath.h>
#else
	#include <stdlib.h>
#endif

double ran1(){
	#ifdef __HAVE_R_ 
		return(runif(0, 1));
	#else
		return((double) rand() / ((double) RAND_MAX + 1.0));
	#endif
} /* End of ran1(). */


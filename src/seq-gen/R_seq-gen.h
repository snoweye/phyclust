//WCC:add
#ifdef __HAVE_R_ 

#include <R.h>
#include <Rinternals.h>

#undef printf
#define printf Rprintf

#undef exit
#define exit(a) error("%d\n", a)

#endif


//WCC:add
#ifdef __HAVE_R_ 

#include <R.h>
#include <Rinternals.h>

#undef printf
#define printf Rprintf

#undef exit
#define exit(a) error("%d\n", a)

#undef puts
#define puts(a) Rprintf(a)

#undef putchar
#define putchar(a) Rprintf(a)

#endif


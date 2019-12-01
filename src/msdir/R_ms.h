//WCC:add
#ifdef __HAVE_R_ 

#include <R.h>
#include <Rinternals.h>

#undef printf
#define printf Rprintf

#undef exit
#define exit(a) error("%d\n", a)

extern const char *R_ms_file_name;
extern FILE *R_ms_file_pointer;

/* In "ms_main.c". */
void ms_main(int argc, char *argv[]);

/* In "R_rand.c". */
double ran1();

#endif


#ifndef __MS_
#define __MS_

struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct devent *nextde;
	} ;
struct c_params {
	int npop;
	int nsam;
	int *config;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
	struct devent *deventlist ;
	} ;
struct m_params {
	 double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
	 } ;
struct params { 
	struct c_params cp;
	struct m_params mp;
	int commandlineseedflag ;
	};

/* In "ms.c". */
#define SITESINC 10
extern unsigned int maxsites;

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

extern double *posit ;
extern double segfac ;
extern int count, ntbs, nseeds ;
extern struct params pars ;	

int gensam(char **list, double *pprobss, double *ptmrca, double *pttot);
void ndes_setup(struct node *ptree, int nsam);
void biggerlist(int nsam, char **list);
char** cmatrix(int nsam, int len);
void locate(int n, double beg, double len, double *ptr);
void getpars(int argc, char *argv[], int *phowmany);
void free_char_2D_AP(char **list, int nsam);
void free_pars();
void argcheck(int arg, int argc, char *argv[]);
void usage();
void addtoelist(struct devent *pt, struct devent *elist);
void free_eventlist(struct devent *pt, int npop);
void make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list);
double ttime(struct node *ptree, int nsam);
double ttimemf(struct node *ptree, int nsam, int mfreq);
void prtree(struct node *ptree, int nsam);
void parens(struct node *ptree, int *descl, int *descr, int noden);
int pickb(int nsam, struct node *ptree, double tt);
int pickbmf(int nsam, int mfreq, struct node *ptree, double tt);
int tdesn(struct node *ptree, int tip, int node);
int pick2(int n, int *i, int *j);
void ordran(int n, double pbuf[]);
int mnmial(int n, int nclass, double p[], int rv[]);
void order(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
int poisso(double u);
double gasdev(double m, double v);

/* In "streec.c". */
#define SEGINC 80;
extern unsigned int seglimit;

struct segl* segtre_mig(struct c_params *cp, int *pnsegs);
int re(int nsam);
int cleftr(int nsam);
int cinr(int nsam, int nsites);
int xover(int nsam,int ic, int is);
int ca(int nsam, int nsites, int c1, int c2);
int isseg(int start, int c, int *psg);
void pick2_chrom(int pop, int config[], int *pc1, int *pc2);
int links(int c);

#endif

//WCC:add
#include "R_ms.h"


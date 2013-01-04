To implement, first copy all the files under the same directory. The
commands for running the score test with the clustered haplotypes are
stored in "score.FDRD.public.r". Please read and run the commands in
"score.FDRD.public.r" line by line, and check if the same results that
are attached at the end of the file could be obtained.


===========================================================================
(1) genodata.mat : an example data matrix (genotypic data in Schaid et
                    al. 2002 format; use "help(haplo.score)" to see the
                    specific format). Note that this file is an R object.

(2) score.FDRD.public.r : the main codes to do data analysis, including
                    Schaid et al. (2002)'s full-dimensional (FD) analysis
                    and our reduced-dimensional (RD) analysis

(3) main.functions.forRDscoreTest.public.r : main functions for RD score test 

(4) functions.r: supporting functions, mainly for calculating the variance
                    of the RD score 

(5) code.getPIstar.getBigMatB : supporting functions, mainly for obtaining
                                 "matrix B" for clustering haplotypes


ps. I tend to use "FD" for full dimensional (original) haplotypes,
"RD" for reduced-dimensional (clustered) haplotypes. Please email me
(jytzeng@stat.ncsu.edu) when running into any questions.




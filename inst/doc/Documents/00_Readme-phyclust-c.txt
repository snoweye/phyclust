Library: phyclust-c
Version: 0.1-5
Date: 2012-04-21
Author: Wei-Chen Chen and Karin S. Dorman
Maintainer: Wei-Chen Chen <wccsnow@gmail.com>
Description: Model Based Phylogenetic Clustering
License: GPL (>= 2)
URL: http://thirteen-01.stat.iastate.edu/snoweye/phyclust/

Examples of main functions are in "test/" and short commands to compile
is in "test/make.test". A quick example is "test_toy.c" and one can
type the following. The "test_toy_se_update.c" is for sequencing error models.

> gcc -std=gnu99 -O3 -Wall -o test_toy test_toy.c ../*.c -I../ -lm
> test_toy Tt.009.anc.015.phy 4 0 0 3 2 0 2 3 10 0
> test_toy Tt.009.anc.015.gap.phy 4 0 0 3 2 0 2 3 10 0

> gcc -std=gnu99 -O3 -Wall -o test_toy_se_update test_toy_se_update.c ../*.c -I../ -lm
> test_toy_se_update Tt.009.anc.015.gap.phy 4 0 0 3 2 0 0 3 10 0

:Wei-Chen Chen

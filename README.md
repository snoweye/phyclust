# phyclust

* **License:** [![License](http://img.shields.io/badge/license-GPL%20v2-orange.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.en.html)
* **Download:** [![Download](http://cranlogs.r-pkg.org/badges/phyclust)](https://cran.r-project.org/package=phyclust)
* **Author:** Wei-Chen Chen and Karin Dorman


Phylogenetic clustering (phyloclustering) is an evolutionary
Continuous Time Markov Chain model-based approach to identify
population structure from molecular data without assuming
linkage equilibrium. The package phyclust (Chen 2011) provides a
convenient implementation of phyloclustering for DNA and SNP data,
capable of clustering individuals into subpopulations and identifying
molecular sequences representative of those subpopulations. It is
designed in C for performance, interfaced with R for visualization,
and incorporates other popular open source programs including
ms (Hudson 2002), seq-gen (Rambaut and Grassly 1997),
Hap-Clustering (Tzeng 2005) and PAML baseml (Yang 1997, 2007), for
simulating data, additional analyses, and searching the best tree.
See the phyclust website for more information, documentations and
examples.



## Installation

phyclust requires
* R version 3.0.0 or higher.
* R package ape.

The package can be installed from the CRAN via the usual
`install.packages("phyclust")`, or via the devtools package:

```r
library(devtools)
install_github("snoweye/phyclust")
```


## Copyright

See phyclust/inst/Documents/ for files in src/msdir/,
src/seq-gen/, src/paml_baseml, and R/ttzeng-*.r.


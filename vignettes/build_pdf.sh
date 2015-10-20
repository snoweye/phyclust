#!/bin/sh

R_HOME=`Rscript -e 'cat(R.home()[1])'`
JSS_BST=${R_HOME}/share/texmf/bibtex/bst/jss.bst
JSS_CLS=${R_HOME}/share/texmf/tex/latex/jss.cls

rm *.aux *.bbl *.blg *.log *.out *.toc *.bst *.cls
cp ${JSS_BST} ./
cp ${JSS_CLS} ./

pdflatex phyclust-guide.Rnw
bibtex phyclust-guide
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
Rscript -e "tools::compactPDF('.', gs_quality='ebook')"
rm *.aux *.bbl *.blg *.log *.out *.toc *.bst *.cls

qpdf phyclust-guide.pdf ../inst/doc/phyclust-guide.pdf
# mv -f *.pdf ../inst/doc/
rm -f *.pdf
cp -f *.Rnw ../inst/doc/
cp -f *.html ../inst/doc/

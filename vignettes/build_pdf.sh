#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex phyclust-guide.Rnw
bibtex phyclust-guide
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
Rscript -e "tools::compactPDF('.', gs_quality='ebook')"
rm *.aux *.bbl *.blg *.log *.out *.toc

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/
cp -f *.html ../inst/doc/

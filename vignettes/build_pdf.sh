#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex phyclust-guide.Rnw
bibtex phyclust-guide
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
pdflatex phyclust-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc

# mv phyclust-guide.pdf phyclust-guide.pdf.org
# qpdf phyclust-guide.pdf.org phyclust-guide.pdf
# rm *.org

mv -f *.pdf ../inst/doc/
cp -f *.Rnw ../inst/doc/
cp -f *.html ../inst/doc/

include ../make.inc
include ../common.inc

PRODUCTS = lmie.aux lmie.bbl lmie.blg lmie.brf lmie.log lmie.out lmie.toc \
           lmie_equations.aux lmie_equations.log lmie_equations.out lmie_equations.toc \
           # lmie.pdf lmie_equations.pdf *.html

all: lmie.pdf lmie_equations.pdf index.html

lmie.pdf: common.tex lmie.tex lmie.bib
	pdflatex -halt-on-error lmie.tex
	bibtex lmie
	pdflatex -halt-on-error lmie.tex
	pdflatex -halt-on-error lmie.tex

index.html: common.tex lmie.tex lmie.bib lmie.pdf
	$(HOME)/src/tth_C-4.03/tth -Llmie -s < lmie.tex | sed 's/\&times;</\&times; </g' | $(HOME)/src/tthgold/tthsplit
	sed 's|href="chap|href="sec|g' index.html > temp.out
	sed 's|Calllmie output format</a><br />|Calllmie output format</a><br /><a href="refs.html">Bibliography</a><br />|g' temp.out > index.html
	rm -f temp.out

lmie_equations.pdf: common.tex lmie_equations.tex
#	pdflatex -halt-on-error lmie_equations.tex
#	bibtex lmie
	pdflatex -halt-on-error lmie_equations.tex
	pdflatex -halt-on-error lmie_equations.tex

clean_doc:
	rm -f lmie.pdf lmie_equations.pdf *.html

include dep.inc

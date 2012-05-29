include ../make.inc

PRODUCTS = lmie.aux lmie.bbl lmie.blg lmie.brf lmie.log lmie.out lmie.toc # lmie.pdf *.html

all: lmie.pdf index.html

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

include dep.inc

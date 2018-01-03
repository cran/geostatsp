all: lgcp.pdf glgm.pdf

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz


%.tex: %.Rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('%<')" $(Rargs)

%.pdf: %.tex
	xelatex $<;
	xelatex $<;
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('$@.pdf', gs_quality='ebook')" 



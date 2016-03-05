all: lgcp.pdf

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz


lgcp.tex: lgcp.Rnw
	'$(R_HOME)/bin/Rscript' -e "knitr::knit('lgcp.Rnw')" $(Rargs)

lgcp.pdf: lgcp.tex
	xelatex lgcp;
	xelatex lgcp;
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('lgcp.pdf', gs_quality='ebook')" 



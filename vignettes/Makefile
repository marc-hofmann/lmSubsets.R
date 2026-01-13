all:
	pdflatex --shell-escape lmSubsets
	bibtex lmSubsets.aux
	pdflatex --shell-escape lmSubsets
	pdflatex --shell-escape lmSubsets
	pdflatex --shell-escape lmSubsets

clean:
	rm -f lmSubsets.tex
	rm -f lmSubsets.aux
	rm -f lmSubsets.bbl
	rm -f lmSubsets.blg
	rm -f lmSubsets.log
	rm -f lmSubsets-*.log
	rm -f lmSubsets.out
	rm -f lmSubsets-*.pdf
	rm -f lmSubsets.upa
	rm -f lmSubsets.upb

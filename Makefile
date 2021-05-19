cit.tar.gz: inst/NEWS 
	R CMD build . && mv cit_*.tar.gz cit.tar.gz
install: cit.tar.gz
	cd .. && R CMD INSTALL cit_.tar.gz
checkv: cit.tar.gz
	R CMD check --use-valgrind cit.tar.gz
check: cit.tar.gz
	R CMD check cit.tar.gz
clean:
	rm cit.tar.gz
debug:
	R -d valgrind --debugger-args='--leak-check=full'
inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
		head -n 80 inst/NEWS
.PHONY: install check clean checkv debug


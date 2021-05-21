VERSION:=$(shell Rscript -e 'x<-readLines("DESCRIPTION");cat(gsub(".+[:]\\s*", "", x[grepl("^Vers", x)]))')

cit_$(VERSION).tar.gz: inst/NEWS R/*.R src/*.cpp
	R CMD build . 

install: cit_$(VERSION).tar.gz
	R CMD INSTALL cit_$(VERSION).tar.gz

checkv: cit_$(VERSION).tar.gz
	R CMD check --use-valgrind cit_$(VERSION).tar.gz
	
check: cit_$(VERSION).tar.gz
	R CMD check --as-cran cit_$(VERSION).tar.gz
clean:
	rm -rf cit*.tar.gz; rm -rf cit*.Rcheck
debug:
	R -d valgrind --debugger-args='--leak-check=full'
inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
		head -n 80 inst/NEWS
version:
	echo $(VERSION)
.PHONY: install check clean checkv debug version


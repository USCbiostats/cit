../cit_.tar.gz: 
	cd .. && R CMD build cit && mv cit_*.tar.gz cit_.tar.gz
install: ../cit_.tar.gz
	cd .. && R CMD INSTALL cit_.tar.gz
checkv: ../cit_.tar.gz
	cd .. && R CMD check --use-valgrind cit_.tar.gz
check: ../cit_.tar.gz
	cd .. && R CMD check cit_.tar.gz
clean:
	rm ../cit_.tar.gz
.PHONY: install check clean checkv


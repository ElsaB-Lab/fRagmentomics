# prepare the package for release
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
R ?= R

build:
	$(R) CMD build .

install:
	$(R) -e 'install.packages("./$(PKGNAME)_$(PKGVERS).tar.gz")'

check:
	$(R) CMD check $(PKGNAME)_$(PKGVERS).tar.gz

manual:
	$(R) -e 'devtools::document();devtools::build_manual(path=".")'

test:
	$(R) -e 'if (any(as.data.frame(devtools::test())[["failed"]] > 0)) stop("Some tests failed.")'

tags: 
	ctags -R R

clean:
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) -f tags
	$(RM) -f $(PKGNAME)_$(PKGVERS).pdf
 
clean-all: clean
	$(RM) -r $(PKGNAME)_$(PKGVERS).tar.gz 

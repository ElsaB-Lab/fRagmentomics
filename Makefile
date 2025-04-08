# prepare the package for release
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
R ?= R

build:
	$(R) CMD build --no-manual .

install: build-cran
	$(R) CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

build-cran:
	$(R) CMD build .

check: build-cran
	$(R) CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran --no-timestamp

manual:
	$(R) -e 'devtools::document();devtools::build_manual(path=".")'

test:
	$(R) -e 'if (any(as.data.frame(devtools::test())[["failed"]] > 0)) stop("Some tests failed.")'

ctags: 
	ctags -R R

clean:
	$(RM) -r $(PKGNAME).Rcheck/
	$(RM) -f tags
	$(RM) -f $(PKGNAME)_$(PKGVERS).pdf
 
clean-all: clean
	$(RM) -r $(PKGNAME)_$(PKGVERS).tar.gz 
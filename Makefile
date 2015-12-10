# Main Makefile

export LIBS = -L/usr/local/lib -larmadillo -framework Accelerate -lgsl -lgslcblas -lm
export INCS = -I/usr/local/include
export FLAGS= -O2

RTDIR = $(shell pwd)

export CPP = g++
export SRCDIR = $(RTDIR)/src
export BINDIR = $(RTDIR)/bin

main: Makefile
	@cd $(SRCDIR) ; make main

prm_vars: Makefile
	@cd $(SRCDIR) ; make prm_vars

.PHONY: clean
clean:
	@cd $(SRCDIR) ; rm *.o *~
	@cd $(BINDIR) ; rm main prm_vars
	rm *~
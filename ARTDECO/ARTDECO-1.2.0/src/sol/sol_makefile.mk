
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -J$(DIRMOD) 
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -module $(DIRMOD)
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+=  -qmoddir=$(DIRMOD)
endif

ALL : $(DIROBJ)/sol_solrad.o

$(DIROBJ)/sol_solrad.o :: sol_solrad.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

clean:
	rm -f $(DIROBJ)/sol_*.o


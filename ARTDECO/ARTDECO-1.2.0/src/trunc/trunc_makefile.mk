
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -J$(DIRMOD) -I$(DIRMOD) 
LOCALFFLAGS+= -fimplicit-none
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
LOCALFFLAGS+= -implicitnone
LOCALFFLAGS+= -module $(DIRMOD)
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -I$(DIRMOD) -qmoddir=$(DIRMOD)
LOCALFFLAGS+= -qundef
LOCALFFLAGS+= -qrealsize=8
endif


ALL : $(DIROBJ)/trunc_deltafit.o $(DIROBJ)/trunc_potter.o $(DIROBJ)/trunc_DM.o

$(DIROBJ)/trunc_deltafit.o :: trunc_deltafit.f90 $(DIROBJ)/trunc_constants.o $(DIROBJ)/trunc_mathsub.o $(DIROBJ)/trunc_compBl.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/trunc_deltafit.o

$(DIROBJ)/trunc_potter.o :: trunc_potter.f90 $(DIROBJ)/trunc_constants.o $(DIROBJ)/trunc_mathsub.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/trunc_potter.o

$(DIROBJ)/trunc_DM.o :: trunc_DM.f90 $(DIROBJ)/trunc_constants.o $(DIROBJ)/trunc_mathsub.o $(DIROBJ)/trunc_compBl.o 
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/trunc_DM.o

$(DIROBJ)/trunc_compBl.o :: trunc_compBl.f90
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/trunc_mathsub.o :: trunc_mathsub.f
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/trunc_constants.o :: trunc_constants.f90
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

clean:
	rm -f $(DIROBJ)/trunc_*.o $(DIRMOD)/mtrunc*

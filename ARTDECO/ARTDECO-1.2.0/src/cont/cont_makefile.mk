
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= #-fdefault-double-8 -fdefault-real-8
LOCALFFLAGS+= -std=legacy 
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= #-real-size 64 -double-size 64
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= #-qrealsize=8
endif

ALL : $(DIROBJ)/cont_continuum.o

$(DIROBJ)/cont_continuum.o :: cont_continuum.f
	$(FC) $(LOCALFFLAGS) -c $< -o $@

clean:
	rm -f $(DIROBJ)/cont_*.o


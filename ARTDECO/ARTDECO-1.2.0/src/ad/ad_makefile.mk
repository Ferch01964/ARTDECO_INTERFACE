
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
LOCALFFLAGS+= -std=legacy 
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -qrealsize=8
endif

ALL : $(DIROBJ)/ad_addoub.o

$(DIROBJ)/ad_addoub.o :: ad_addoub.f $(DIROBJ)/ad_supermatrix.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/ad_addoub.o

$(DIROBJ)/ad_supermatrix.o :: ad_supermatrix.f    
	$(FC) $(LOCALFFLAGS) -c $< -o $@

clean:
	rm -f $(DIROBJ)/ad_*.o


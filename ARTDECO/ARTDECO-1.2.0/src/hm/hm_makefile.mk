
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -fimplicit-none 
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
LOCALFFLAGS+= -implicitnone
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -qrealsize=8
LOCALFFLAGS+= -qundef
endif


ALL : $(DIROBJ)/hm_rhm1p1.o $(DIROBJ)/hm_ihm2p0.o $(DIROBJ)/hm_phm2p4.o

$(DIROBJ)/hm_phm2p4.o :: hm_phm2p4.f90 $(DIROBJ)/hm_utility.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/hm_phm2p4.o

$(DIROBJ)/hm_rhm1p1.o :: hm_rhm1p1.f90 $(DIROBJ)/hm_utility.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/hm_rhm1p1.o

$(DIROBJ)/hm_ihm2p0.o :: hm_ihm2p0.f90 $(DIROBJ)/hm_utility.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/hm_ihm2p0.o

$(DIROBJ)/hm_utility.o :: hm_utility.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

clean:
	rm -f $(DIROBJ)/hm_phm2p4.o

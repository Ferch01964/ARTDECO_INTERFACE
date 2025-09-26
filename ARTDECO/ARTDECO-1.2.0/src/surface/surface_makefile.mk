
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -J$(DIRMOD) -I$(DIRMOD) 
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -module $(DIRMOD)
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -I$(DIRMOD) -qmoddir=$(DIRMOD)
endif

ALL : $(DIROBJ)/surface_brdf.o 

$(DIROBJ)/surface_brdf.o : surface_brdf.f90 $(DIROBJ)/brdf_ocean.o $(DIROBJ)/brdf_pavel.o  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/brdf_ocean.o : brdf_ocean.f90 $(DIROBJ)/artdeco_constants.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/brdf_pavel.o : brdf_pavel.f90 $(DIROBJ)/artdeco_constants.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/artdeco_constants.o :: $(MAINDIR)/artdeco_constants.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_common.o :: $(MAINDIR)/artdeco_common.f90  $(DIROBJ)/artdeco_constants.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@


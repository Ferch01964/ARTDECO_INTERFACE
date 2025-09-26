
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

ALL : $(DIROBJ)/mcrad1d_main.o 

$(DIROBJ)/mcrad1d_constants.o : mcrad1d_constants.f90
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/mcrad1d_constants.o 

$(DIROBJ)/mcrad1d_common.o : mcrad1d_common.f90 $(DIROBJ)/mcrad1d_constants.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/mcrad1d_common.o 

$(DIROBJ)/mcrad1d_ziggurat.o : mcrad1d_ziggurat.f90
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/mcrad1d_ziggurat.o 

$(DIROBJ)/mcrad1d_interpolate.o : mcrad1d_interpolate.f90
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/mcrad1d_interpolate.o 

$(DIROBJ)/mcrad1d_submain.o : mcrad1d_submain.f90 $(DIROBJ)/mcrad1d_interpolate.o $(DIROBJ)/mcrad1d_constants.o $(DIROBJ)/mcrad1d_ziggurat.o $(DIROBJ)/mcrad1d_common.o $(DIROBJ)/surface_brdf.o
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/mcrad1d_submain.o 

$(DIROBJ)/mcrad1d_main.o : mcrad1d_main.f90 $(DIROBJ)/mcrad1d_interpolate.o $(DIROBJ)/mcrad1d_constants.o $(DIROBJ)/mcrad1d_ziggurat.o $(DIROBJ)/mcrad1d_submain.o $(DIROBJ)/surface_brdf.o $(DIROBJ)/mcrad1d_common.o
	$(FC) $(LOCALFFLAGS) $(FLAGOMP) -c $< -o $(DIROBJ)/mcrad1d_main.o 

# other sub programs
$(DIROBJ)/surface_brdf.o ::
	cd $(DIR_SURF)       && $(MAKE) -f surface_makefile.mk FC='$(FC)' FFLAGS='$(FFLAGS)'   DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) MAINDIR=$(MAINDIR)  && cd ..

clean:
	rm -f $(DIROBJ)/mcrad1d_*.o $(DIRMOD)/mcrad1d*


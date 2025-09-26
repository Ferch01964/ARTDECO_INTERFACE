
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -I$(DIRMOD) 
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
LOCALFFLAGS+= -module $(DIRMOD)
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -I$(DIRMOD)
LOCALFFLAGS+= -qrealsize=8
endif

ALL : $(DIROBJ)/doad_doad.o

$(DIROBJ)/doad_doad.o :: doad_doad.f90 $(DIROBJ)/doad_base.o $(DIROBJ)/doad_rgrndm.o $(DIROBJ)/doad_ssurf.o $(DIROBJ)/doad_gau_xmu.o $(DIROBJ)/doad_get_surf.o $(DIROBJ)/doad_plkavg.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_get_surf.o :: doad_get_surf.f90 $(DIROBJ)/doad_base.o $(DIROBJ)/surface_brdf.o  
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_plkavg.o :: doad_plkavg.f90 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_base.o :: doad_base.f 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_ssurf.o :: doad_ssurf.f 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_rgrndm.o :: doad_rgrndm.f 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/doad_gau_xmu.o :: doad_gau_xmu.f 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

# other sub programs
$(DIROBJ)/surface_brdf.o ::
	cd $(DIR_SURF)       && $(MAKE) -f surface_makefile.mk FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) MAINDIR=$(MAINDIR)  && cd ..


clean:
	rm -f $(DIROBJ)/doad_*.o


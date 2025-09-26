
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -I$(DIRMOD) 
LOCALFFLAGS+= -fimplicit-none 
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
LOCALFFLAGS+= -implicitnone
LOCALFFLAGS+= -module $(DIRMOD)
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -I$(DIRMOD)
LOCALFFLAGS+= -qundef
LOCALFFLAGS+= -qrealsize=8
endif

ALL :: $(DIROBJ)/od_disort.o

$(DIROBJ)/od_disortutils.o :: od_disortutils.f
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/od_rdi1mach.o :: od_rdi1mach.f
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/od_bdref.o :: od_bdref.f90 $(DIROBJ)/surface_brdf.o 
	$(FC) $(LOCALFFLAGS)  -c $< -o $@

$(DIROBJ)/od_disort.o :: od_disort.f $(DIROBJ)/od_disortutils.o $(DIROBJ)/od_rdi1mach.o $(DIROBJ)/od_bdref.o 
	$(FC) $(LOCALFFLAGS)  -c $< -o $(DIROBJ)/od_disort.o

# other sub programs
$(DIROBJ)/surface_brdf.o ::
	cd $(DIR_SURF)       && $(MAKE) -f surface_makefile.mk FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) MAINDIR=$(MAINDIR) && cd ..

clean:
	rm -f $(DIROBJ)/od_*.o

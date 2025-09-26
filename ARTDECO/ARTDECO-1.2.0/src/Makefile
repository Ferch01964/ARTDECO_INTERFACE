
###################
# Choose compiler
FC = gfortran
#FC = ifort
##################
####  Line added on 24/03/14 for make tar
PROJECT_NAME=ARTDECO
VERSION=1.1.1
TAR_NAME=$(PROJECT_NAME)_V$(VERSION)

DIR_ARTDECO = $(CURDIR)/..
#/home/winiarek/ARTDECO

DIROBJ  = $(DIR_ARTDECO)/obj
DIRMOD  = $(DIR_ARTDECO)/mod
DIRF2PY = $(DIR_ARTDECO)/f2py

DIR_SRC         = $(DIR_ARTDECO)/src
DIR_AD          = $(DIR_SRC)/ad
DIR_OD          = $(DIR_SRC)/disort
DIR_MIE3        = $(DIR_SRC)/mie3
DIR_TRUNC       = $(DIR_SRC)/trunc
DIR_HM          = $(DIR_SRC)/hm
DIR_NM          = $(DIR_SRC)/numerical_methods
DIR_MCRAD1D     = $(DIR_SRC)/mcrad1d
DIR_SURF        = $(DIR_SRC)/surface
DIR_DOAD        = $(DIR_SRC)/doad
DIR_CONT        = $(DIR_SRC)/cont
DIR_SOL         = $(DIR_SRC)/sol

EXE    = artdeco

########## gfortran compiler options #############
ifeq ($(FC), gfortran)

#### PARANOID Last step before mental diseases hospital :-)
#FFLAGS =-g
#FFLAGS+=-Wall -Waliasing -Wampersand -Wcharacter-truncation -Wconversion
#FFLAGS+=-Wimplicit-interface -Wline-truncation -Wsurprising
#FFLAGS+=-Wunused-parameter 
#FFLAGS+=-fbounds-check
#FFLAGS+=-O0 -m64
#FFLAGS+=-fno-omit-frame-pointer

#### Development purposes (-g -O0)
# FFLAGS =-O0 -m64
# FFLAGS+=-g
# FFLAGS+=-Wall -Wextra -Wuninitialized
# FFLAGS+=-fbounds-check
# FFLAGS+=-fno-omit-frame-pointer 
# FFLAGS+=-frecursive

### Production mode 
FFLAGS = -m64 -O2 -fbounds-check
#FFLAGS = -m32 -O2 -fbounds-check
# Necessary for F2PY module creation
FFLAGS += -fPIC -cpp 


# artdeco_* files will be compiled with these options :
LOCALFFLAGS = $(FFLAGS)
LOCALFFLAGS += -J$(DIRMOD) -I$(DIRMOD) -fimplicit-none 

# flag tor OpenMP
#FLAGOMP= -fopenmp

endif


########## ifort compiler options #############
ifeq ($(FC), ifort)

# other 'library' (or non artdeco_ developped files)
# will be compiled with those options: 
### Production mode
FFLAGS = -m64 -O2

# Necessary for F2PY module creation
FFLAGS += -fpp  -assume protect_parens,minus0,byterecl

# artdeco_* files will be compiled with these options :
LOCALFFLAGS = $(FFLAGS)
LOCALFFLAGS += -module $(DIRMOD) -implicitnone

# flag tor OpenMP
FLAGOMP= -openmp

endif

########## xlf compiler options #############
ifeq ($(FC), xlf90)
#############################################
#                                           #
#             NOT WORKING YET               #
#                                           #
#############################################
# other 'library' (or non artdeco developped files)
# will be compiled with those options: 
### Production mode
FFLAGS = -q64 -qnosave -qsuffix=f=f90 

# artdeco_* files will be compiled with these options :
LOCALFFLAGS = $(FFLAGS)
LOCALFFLAGS += -I$(DIRMOD) -qmoddir=$(DIRMOD) -qundef

# flag tor OpenMP
FLAGOMP= -qreport -qsmp=noauto

endif

OBJ_SOL = $(DIROBJ)/sol_solrad.o

OBJ_CONT = \
      $(DIROBJ)/cont_continuum.o

OBJ_DOAD = \
      $(DIROBJ)/doad_doad.o \
      $(DIROBJ)/doad_base.o \
      $(DIROBJ)/doad_rgrndm.o \
      $(DIROBJ)/doad_ssurf.o \
      $(DIROBJ)/doad_gau_xmu.o \
      $(DIROBJ)/doad_plkavg.o \
      $(DIROBJ)/doad_get_surf.o

OBJ_SURF = \
      $(DIROBJ)/surface_brdf.o \
      $(DIROBJ)/brdf_ocean.o \
      $(DIROBJ)/brdf_pavel.o


OBJ_NM = \
      $(DIROBJ)/nm_get_index.o \
      $(DIROBJ)/nm_util.o \
      $(DIROBJ)/nm_sort.o \
      $(DIROBJ)/nm_sort2.o \
      $(DIROBJ)/nm_gauss_nodes.o \
      $(DIROBJ)/nm_nm.o \
      $(DIROBJ)/nm_type.o


OBJ_HM = \
      $(DIROBJ)/hm_phm2p4.o \
      $(DIROBJ)/hm_rhm1p1.o \
      $(DIROBJ)/hm_ihm2p0.o \
      $(DIROBJ)/hm_utility.o 

OBJ_TRUNC = \
      $(DIROBJ)/trunc_constants.o \
      $(DIROBJ)/trunc_compBl.o \
      $(DIROBJ)/trunc_mathsub.o \
      $(DIROBJ)/trunc_deltafit.o \
      $(DIROBJ)/trunc_potter.o \
      $(DIROBJ)/trunc_DM.o

OBJ_MIE3 = \
      $(DIROBJ)/mie3.o 

OBJ_AD =  \
      $(DIROBJ)/ad_addoub.o \
      $(DIROBJ)/ad_supermatrix.o 

OBJ_MCRAD1D = \
      $(DIROBJ)/mcrad1d_ziggurat.o \
      $(DIROBJ)/mcrad1d_common.o \
      $(DIROBJ)/mcrad1d_constants.o \
      $(DIROBJ)/mcrad1d_interpolate.o \
      $(DIROBJ)/mcrad1d_submain.o \
      $(DIROBJ)/mcrad1d_main.o 

OBJ_OD =  \
      $(DIROBJ)/od_disortutils.o \
      $(DIROBJ)/od_rdi1mach.o \
      $(DIROBJ)/od_bdref.o \
      $(DIROBJ)/od_disort.o

OBJ_ARTDECO = \
      $(DIROBJ)/artdeco_constants.o \
      $(DIROBJ)/artdeco_common.o \
      $(DIROBJ)/artdeco_utility.o \
      $(DIROBJ)/artdeco_inout.o \
      $(DIROBJ)/artdeco_get_opt.o \
      $(DIROBJ)/artdeco_get_betal.o \
      $(DIROBJ)/artdeco_solrad.o \
      $(DIROBJ)/artdeco_layers.o \
      $(DIROBJ)/artdeco_call_ad.o \
      $(DIROBJ)/artdeco_call_od.o \
      $(DIROBJ)/artdeco_call_mcrad1d.o \
      $(DIROBJ)/artdeco_sinsca.o \
      $(DIROBJ)/artdeco_call_doad.o \
      $(DIROBJ)/artdeco_band_int.o

OBJ_ALL = $(OBJ_ARTDECO) \
          $(OBJ_NM) \
          $(OBJ_MIE3) \
          $(OBJ_HM) \
          $(OBJ_TRUNC) \
          $(OBJ_AD) \
          $(OBJ_OD) \
          $(OBJ_MCRAD1D) \
          $(OBJ_SURF) \
          $(OBJ_DOAD) \
	  $(OBJ_CONT) \
          $(OBJ_SOL)

# This are the source files from ARTDECO
# containing fortran modules to be included into
# the F2PY module

OBJ_F2PY = $(OBJ_NM) \
          $(OBJ_MIE3) \
          $(OBJ_HM) \
          $(OBJ_TRUNC) \
          $(OBJ_AD) \
          $(OBJ_OD) \
          $(OBJ_MCRAD1D) \
          $(OBJ_SURF) \
          $(OBJ_DOAD) \
	  $(OBJ_CONT) \
          $(OBJ_SOL)

SRC_ARTDECO_F2PY = \
      $(DIR_SRC)/artdeco_constants.f90 \
      $(DIR_SRC)/artdeco_common.f90 \
      $(DIR_SRC)/artdeco_utility.f90 \
      $(DIR_SRC)/artdeco_f2py_utility.f90 \
      $(DIR_SRC)/artdeco_solrad.f90 \
      $(DIR_SRC)/artdeco_inout.f90 \
      $(DIR_SRC)/artdeco_layers.f90 \
      $(DIR_SRC)/artdeco_sinsca.f90 \
      $(DIR_SRC)/artdeco_call_doad.f90 \
      $(DIR_SRC)/artdeco_band_int.f90 \
      $(DIR_SRC)/artdeco_call_od.f90 \
      $(DIR_SRC)/artdeco_get_betal.f90 \
      $(DIR_SRC)/artdeco_call_mcrad1d.f90 



all: $(EXE)
	@echo " "
	@echo "That's it!"
	@echo " "

f2py :: $(OBJ_ARTDECO) $(DIROBJ)/artdeco_f2py_utility.o
	rm -rf $(DIRF2PY)/*
	cat $(SRC_ARTDECO_F2PY) > $(DIRF2PY)/artdeco_f2py.f90
	f2py $(DIRF2PY)/artdeco_f2py.f90 -h $(DIRF2PY)/artdeco_f2py.pyf
	f2py -I$(DIRMOD) --verbose --fcompiler=$(FC) --f90flags="-fbounds-check -O2 -m64"  -m artdeco -c $(OBJ_F2PY) $(DIRF2PY)/artdeco_f2py.f90
	mv artdeco.so  $(DIRF2PY)/
	@echo " "
	@echo "make f2py DONE!"
	@echo " "

$(EXE) :: artdeco.f90 $(OBJ_ARTDECO)
	$(FC) $(LOCALFFLAGS) $(FLAGOMP) $(OBJ_ALL) artdeco.f90 -o $(EXE)

$(DIROBJ)/artdeco_constants.o :: artdeco_constants.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_common.o :: artdeco_common.f90  $(DIROBJ)/artdeco_constants.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_utility.o :: artdeco_utility.f90  $(DIROBJ)/artdeco_constants.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_f2py_utility.o :: artdeco_f2py_utility.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@



$(DIROBJ)/artdeco_inout.o :: artdeco_inout.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/artdeco_utility.o \
	$(DIROBJ)/nm_sort2.o $(DIROBJ)/nm_sort.o $(DIROBJ)/surface_brdf.o $(DIROBJ)/artdeco_layers.o $(DIROBJ)/nm_get_index.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_solrad.o :: artdeco_solrad.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/artdeco_utility.o \
                               $(DIROBJ)/sol_solrad.o 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_layers.o :: artdeco_layers.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/artdeco_utility.o \
				$(DIROBJ)/nm_gauss_nodes.o $(DIROBJ)/cont_continuum.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_call_ad.o :: artdeco_call_ad.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/ad_addoub.o \
	                       $(DIROBJ)/artdeco_layers.o 
	$(FC) $(LOCALFFLAGS) $(FLAGOMP) -c $< -o $@

$(DIROBJ)/artdeco_call_mcrad1d.o :: artdeco_call_mcrad1d.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/mcrad1d_main.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_band_int.o :: artdeco_band_int.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/artdeco_utility.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_sinsca.o :: artdeco_sinsca.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/artdeco_layers.o $(DIROBJ)/surface_brdf.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_call_od.o :: artdeco_call_od.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/od_disort.o \
	                       $(DIROBJ)/artdeco_layers.o $(DIROBJ)/artdeco_sinsca.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_get_betal.o :: artdeco_get_betal.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/trunc_deltafit.o \
				$(DIROBJ)/trunc_DM.o $(DIROBJ)/trunc_potter.o $(DIROBJ)/nm_gauss_nodes.o 
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/artdeco_get_opt.o :: artdeco_get_opt.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o \
				 $(DIROBJ)/mie3.o $(DIROBJ)/hm_phm2p4.o $(DIROBJ)/hm_rhm1p1.o $(DIROBJ)/hm_ihm2p0.o
	$(FC) $(LOCALFFLAGS) $(FLAGOMP) -c $< -o $@

$(DIROBJ)/artdeco_call_doad.o :: artdeco_call_doad.f90  $(DIROBJ)/artdeco_constants.o $(DIROBJ)/artdeco_common.o $(DIROBJ)/doad_doad.o \
	                       $(DIROBJ)/artdeco_layers.o $(DIROBJ)/artdeco_sinsca.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

# used sub-program :

$(DIROBJ)/sol_solrad.o :: 
	cd $(DIR_SOL)      && $(MAKE) -f sol_makefile.mk     FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/surface_brdf.o ::
	cd $(DIR_SURF)     && $(MAKE) -f surface_makefile.mk FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) MAINDIR=$(DIR_SRC) && cd ..

$(DIROBJ)/cont_continuum.o :: 
	cd $(DIR_CONT)     && $(MAKE) -f cont_makefile.mk    FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/doad_doad.o ::
	cd $(DIR_DOAD)     && $(MAKE) -f doad_makefile.mk    FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) DIR_SURF=$(DIR_SURF) MAINDIR=$(DIR_SRC) && cd ..

$(DIROBJ)/mcrad1d_main.o ::
	cd $(DIR_MCRAD1D)  && $(MAKE) -f mcrad1d_makefile.mk FC='$(FC)' FFLAGS='$(FFLAGS)' FLAGOMP=$(FLAGOMP) DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) DIR_SURF=$(DIR_SURF) MAINDIR=$(DIR_SRC) && cd ..

$(DIROBJ)/mie3.o ::
	cd $(DIR_MIE3)     && $(MAKE) -f mie3_makefile.mk    FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD)  && cd ..

$(DIROBJ)/ad_addoub.o ::
	cd $(DIR_AD)       && $(MAKE) -f ad_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/od_disort.o :: 
	cd $(DIR_OD)       && $(MAKE) -f od_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) DIR_SURF=$(DIR_SURF) MAINDIR=$(DIR_SRC) && cd ..

$(DIROBJ)/trunc_deltafit.o ::
	cd $(DIR_TRUNC)    && $(MAKE) -f trunc_makefile.mk   FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/trunc_DM.o ::
	cd $(DIR_TRUNC)    && $(MAKE) -f trunc_makefile.mk   FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/trunc_potter.o ::
	cd $(DIR_TRUNC)    && $(MAKE) -f trunc_makefile.mk   FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/hm_phm2p4.o :: 
	cd $(DIR_HM)       && $(MAKE) -f hm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/hm_rhm1p1.o :: 
	cd $(DIR_HM)       && $(MAKE) -f hm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/hm_ihm2p0.o ::
	cd $(DIR_HM)       && $(MAKE) -f hm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)' DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD) && cd ..

$(DIROBJ)/nm_sort.o ::
	cd $(DIR_NM)       && $(MAKE) -f nm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD)  && cd ..

$(DIROBJ)/nm_sort2.o ::
	cd $(DIR_NM)       && $(MAKE) -f nm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD)  && cd ..

$(DIROBJ)/nm_gauss_nodes.o ::
	cd $(DIR_NM)       && $(MAKE) -f nm_makefile.mk      FC='$(FC)' FFLAGS='$(FFLAGS)'  DIROBJ=$(DIROBJ) DIRMOD=$(DIRMOD)  && cd ..

##### Line added on 24/03/14
DIR_TAR=${HOME}/tmp
tar : cleanall
	rm -rf $(DIR_TAR)
	mkdir -p  $(DIR_TAR)
	echo "Creation of the archive $(TAR_NAME).tar.gz in ${DIR_TAR}"
	rm -rf $(DIR_TAR)/$(TAR_NAME)
	mkdir -p $(DIR_TAR)/$(TAR_NAME)
	cp -rp ../* $(DIR_TAR)/$(TAR_NAME)
	rm -f $(DIR_TAR)/$(TAR_NAME)/*.tar.gz
	tar zcf $(DIR_TAR)/$(TAR_NAME).tar.gz -C $(DIR_TAR)/ $(TAR_NAME)
	rm -rf $(DIR_TAR)/$(TAR_NAME)    	

check: 
	 ../test/GAME/run_artdeco_vs_game.sh

.PHONY: check


clean:
	rm -rf $(DIROBJ)/*.o $(DIRMOD)/*.mod *.so artdeco.dSYM/  $(DIRF2PY)/*

cleanall:
	rm -rf $(DIROBJ)/*.o $(DIRMOD)/*.mod $(EXE) *.so artdeco.dSYM/  $(DIRF2PY)/* python/*.pyc

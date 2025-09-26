LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -J$(DIRMOD) -I$(DIRMOD) 
LOCALFFLAGS+= -fimplicit-none
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -module $(DIRMOD)
LOCALFFLAGS+= -implicitnone
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -I$(DIRMOD) -qmoddir=$(DIRMOD)
LOCALFFLAGS+= -qundef
endif

ALL :  $(DIROBJ)/nm_sort.o  $(DIROBJ)/nm_sort2.o $(DIROBJ)/nm_gauss_nodes.o



$(DIROBJ)/nm_nm.o :: nm_nm.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/nm_type.o :: nm_type.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/nm_util.o :: nm_util.f90
	$(FC) $(LOCALFFLAGS) -c $< -o $@

$(DIROBJ)/nm_gauss_nodes.o :: nm_gauss_nodes.f90 $(DIROBJ)/nm_type.o $(DIROBJ)/nm_util.o  $(DIROBJ)/nm_nm.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/nm_gauss_nodes.o

$(DIROBJ)/nm_sort.o :: nm_sort.f90 $(DIROBJ)/nm_type.o $(DIROBJ)/nm_util.o $(DIROBJ)/nm_nm.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/nm_sort.o

$(DIROBJ)/nm_sort2.o :: nm_sort2.f90 $(DIROBJ)/nm_type.o $(DIROBJ)/nm_util.o $(DIROBJ)/nm_nm.o  $(DIROBJ)/nm_get_index.o
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/nm_sort2.o    

$(DIROBJ)/nm_get_index.o :: nm_get_index.f90 $(DIROBJ)/nm_type.o $(DIROBJ)/nm_util.o $(DIROBJ)/nm_nm.o
	$(FC) $(LOCALFFLAGS) -c $< -o $@

clean:
	rm -f $(DIROBJ)/nm_*.o $(DIRMOD)/*nm_*

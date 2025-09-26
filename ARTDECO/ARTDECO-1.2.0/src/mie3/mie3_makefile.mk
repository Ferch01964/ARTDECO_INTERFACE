
LOCALFFLAGS = $(FFLAGS) 

ifeq ($(FC), gfortran)
LOCALFFLAGS+= -std=legacy 
LOCALFFLAGS+= -fdefault-double-8 -fdefault-real-8
endif
ifeq ($(FC), ifort)
LOCALFFLAGS+= -real-size 64 -double-size 64
endif
ifeq ($(FC), xlf90)
LOCALFFLAGS+= -qrealsize=8
endif


ALL : $(DIROBJ)/mie3.o 

$(DIROBJ)/mie3.o :: mie3.f 
	$(FC) $(LOCALFFLAGS) -c $< -o $(DIROBJ)/mie3.o

clean:
	rm -f $(DIROBJ)/mie3.o

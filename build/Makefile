#   lotka volterra
#------------------------------------
#             makefile
#------------------------------------
#   author :: perrin w. davidson
#   email  :: pwd@uchicago.edu
#------------------------------------

#	set compiler ::
FCOMP = gfortran

#	get base path ::
BASEPATH = $(realpath $(CURDIR)/..)

#	get build and source paths ::
BUILDDIR = $(BASEPATH)/build
SRCDIR = $(BASEPATH)/code

#	set main.F and main.exe path ::
TARGET = $(BUILDDIR)/lotka_volterra
SRC = $(SRCDIR)/lotka_volterra

#	set debug flags :: 
FCFLAGS = -g -ffpe-trap=zero,invalid,overflow,underflow 

#	set main.o path ::
OBJECTS = $(TARGET).o

#	compile main executable ::
$(TARGET): $(OBJECTS)
	@$(FCOMP) $(FCFLAGS) -o $(TARGET) $(OBJECTS)

#	compile main.o ::
$(TARGET).o: $(SRC).F90
	@$(FCOMP) -c $(SRC).F90

#	set clean ::
.PHONY: clean
clean:
	@rm -f $(TARGET) $(OBJECTS)

#	end makefile

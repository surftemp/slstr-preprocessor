
ifdef DEBUG
    COMPILE_FLAGS = -g -fbounds-check
else
    COMPILE_FLAGS = -O3
endif

# Customise for your system, you will need netcdf libraries installed
FORTRAN = gfortran
INCDIRS = -I/usr/local/include
LIBDIRS = -L/usr/local/lib

all: Preprocess_SLSTR TestNeighbourhoods

SLSTR_Preprocessor.o: SLSTR_Preprocessor.f90 GbcsPath.o
	$(FORTRAN) $(COMPILE_FLAGS) $(INCDIRS) -c -o SLSTR_Preprocessor.o SLSTR_Preprocessor.f90

Preprocess_SLSTR.o: Preprocess_SLSTR.f90 SLSTR_Preprocessor.o GbcsPath.o
	$(FORTRAN) $(COMPILE_FLAGS) $(INCDIRS) -c -o Preprocess_SLSTR.o Preprocess_SLSTR.f90

TestNeighbourhoods.o: TestNeighbourhoods.f90 SLSTR_Preprocessor.o GbcsPath.o
	$(FORTRAN) $(COMPILE_FLAGS) $(INCDIRS) -c -o TestNeighbourhoods.o TestNeighbourhoods.f90

GbcsPath.o: GbcsPath.f90
	$(FORTRAN) $(COMPILE_FLAGS) $(INCDIRS) -c -o GbcsPath.o GbcsPath.f90

Preprocess_SLSTR: SLSTR_Preprocessor.o Preprocess_SLSTR.o GbcsPath.o
	$(FORTRAN) $(LIBDIRS) $(LDFLAGS) -lnetcdff Preprocess_SLSTR.o SLSTR_Preprocessor.o GbcsPath.o -o Preprocess_SLSTR

TestNeighbourhoods: SLSTR_Preprocessor.o TestNeighbourhoods.o GbcsPath.o
	$(FORTRAN) $(LIBDIRS) $(LDFLAGS) -lnetcdff TestNeighbourhoods.o SLSTR_Preprocessor.o GbcsPath.o -o TestNeighbourhoods

clean:
	rm -f Preprocess_SLSTR
	rm -f TestNeighbourhoods
	rm -f *.o
	rm -f *.mod

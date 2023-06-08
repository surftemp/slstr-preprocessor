# Customise for your system, you will need netcdf libraries installed
FC = gfortran
INCLUDES = -I/usr/lib64/gfortran/modules
LDFLAGS = -L/usr/local/lib
LDLIBS = -lnetcdff

ifdef DEBUG
    FC_FLAGS = -Og -g -fcheck=all
else
    FC_FLAGS = -O2
endif


all: s3testpp s3regrid

clean:
	$(RM) s3testpp s3regrid
	$(RM) *.o *.mod

s3testpp: s3testpp.o SLSTR_Preprocessor.o GbcsPath.o
s3regrid: s3regrid.o GbcsKinds.o GbcsNetCDF.o GbcsPath.o SLSTR_Preprocessor.o

# Dependencies
SLSTR_Preprocessor.o: SLSTR_Preprocessor.f90 GbcsPath.o
s3testpp.o: s3testpp.f90 SLSTR_Preprocessor.o GbcsPath.o
GbcsPath.o: GbcsPath.f90
GbcsKinds.o: GbcsKinds.f90
GbcsNetCDF.o: GbcsNetCDF.f90
s3regrid.o: s3regrid.f90 GbcsKinds.o GbcsNetCDF.o GbcsPath.o SLSTR_Preprocessor.o


# Default rules
%.o : %.f90
	$(FC) $(FC_FLAGS) $(INCLUDES) -c $< -o $@

# Always link using Fortran
LINK.o = $(FC) $(FC_FLAGS) $(LDFLAGS)
# Alternatively we could use explicit: $(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@
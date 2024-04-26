# Compiler settings - Can change to your preferred compiler (e.g., gfortran, ifort)
FC = gfortran
# Compiler flags
FFLAGS = -O2

# Build target
TARGET = electron_fluxes

# Object files
OBJS = subroutines.o electron_fluxes.o

# Build rule
$(TARGET): $(OBJS)
	$(FC) -o $@ $(OBJS)

# Compilation rules
electron_fluxes.o: electron_fluxes.f90 subroutines.o
	$(FC) $(FFLAGS) -c electron_fluxes.f90

subroutines.o: subroutines.f90
	$(FC) $(FFLAGS) -c subroutines.f90

# Clean rule
clean:
	rm -f $(TARGET) $(OBJS)

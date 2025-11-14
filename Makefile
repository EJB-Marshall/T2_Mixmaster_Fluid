FC = gfortran
FFLAGS = -O3 -Wall -I/opt/homebrew/include #\
		-I$(shell brew --prefix openblas)/include

LDFLAGS = -L/opt/homebrew/lib \
		-lhdf5_fortran -lhdf5_f90cstub -lhdf5 #\
		-L$(shell brew --prefix openblas)/lib #\
		-lopenblas

SRC = real_type_mod.o read_parameters_mod.o \
	diff_op_mod.o finite_difference_mod.o boundary_conditions_mod.o \
	T2_fluid_mod.o time_stepping_mod.o hdf_io_mod.o \
	evolve_system.o

.PHONY: clean

main: $(SRC)
	$(FC) $(FFLAGS) -o $@ $(SRC) $(LDFLAGS)
%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
clean:
	@rm *.o *.mod main
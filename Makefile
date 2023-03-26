TARGET = lattice_density_functional_theory.out
OBJECTS = main_lattice_density_functional_theory.o sub_boundary_condition.o sub_update_density.o sub_check_convergence.o sub_calc_ave_density.o sub_calc_omega.o
F90 = ifort -qopenmp
FFLAGS = -fpp -Ddebug -check all -warn declarations -CB -fpe0 -traceback -integer-size 64

# FFLAGS = -fpp -check all -warn declarations -CB -fpe0 -traceback -integer-size 64
# FFLAGS = -check all -warn declarations -CB -fpe0 -traceback -integer-size 64
# FLAGS = -check -fpe0 -traceback
#COMMON_MOD = modules.f90 

#  -fpp : release mode,  -fpp -Ddebug : debug mode

.SUFFIXES :
.SUFFIXES : .o .f90 
.f90.o:
	${F90}	-c $< ${FFLAGS}

${TARGET} : ${OBJECTS}
	${F90}	-o $@ ${OBJECTS}
clean :
	rm  *.o ${TARGET}


# ${OBJECTS} : ${COMMON_MOD}

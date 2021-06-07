# MAKEFILE FOR EULER

# COMPILER
cc=gfortran

# LIBRARY LOCATION
lb_lc= -I/home/sugan/fftw/include

# LIBRARY FILE
lb_fftw=-L/home/sugan/fftw/lib -lfftw3 -lm

# PROGRAM
program=euler3D_system.f90

# MODULES
timer_mod=system_timer.f90
fft_mod=system_fftw.f90
pvdoutput_mod=system_pvdoutput.f90
constants_mod=system_constants.f90
auxilaries_mod=system_auxilaries.f90
variables_mod=system_variables.f90
initialcondition_mod=system_initialcondition.f90
functions_mod=system_functions.f90
asolver_mod=system_advectionsolver.f90
vsolver_mod=system_vorticitysolver.f90
test_mod=system_test.f90
output_mod=system_output.f90
main_mod=system_main.f90

# OBJECTS
obj=system_timer.o\
	system_fftw.o\
	system_pvdoutput.o\
	system_constants.o\
	system_auxilaries.o\
	system_variables.o\
	system_initialcondition.o\
	system_functions.o\
	system_advectionsolver.o\
	system_vorticitysolver.o\
	system_output.o\
	system_test.o\
	system_main.o

# EXECUTABLE
run=./ex

# CLEAN COMMANDS
rmex=rm ex

mkcl=make cl
#----------------------------end-------

# MAKEFILE
# ---------------------------start-----
ex:$(ob)
	$(cc) $(lb_lc) -c $(fft_mod) $(lb_fftw)
	$(cc) -c $(timer_mod)
	$(cc) -c $(pvdoutput_mod)
	$(cc) -c $(constants_mod)
	$(cc) -c $(auxilaries_mod)
	$(cc) -c $(variables_mod)
	$(cc) -c $(initialcondition_mod)
	$(cc) -c $(functions_mod)
	$(cc) -c $(output_mod)
	$(cc) -c $(asolver_mod)
	$(cc) -c $(vsolver_mod)
	$(cc) -c $(test_mod)
	$(cc) -c $(main_mod)
	$(cc) $(lb_lc) $(program) $(obj) $(lb_fftw) -o ex
	$(mkcl)
	$(run)

#----------------------------end-------

# CLEANING
# ---------------------------start-----
clean:
	rm ex
	clear
cl:
	rm *.mod
	rm *.o
#----------------------------end-------

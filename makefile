CMP     = gfortran
OPTFLAG = -O3
FREEFORM = -ffree-form
F132FORM = -ffixed-line-length-132
LAPACK =   -framework accelerate
ARPACK =  -L/opt/homebrew/lib/ -larpack
INCLUDE =  -I/opt/homebrew/include
WARN = #-Wall

4BodySVD.x:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.o
	${CMP} ${WARN} ${OPT} ${FREEFORM} matrix_stuff.o 4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90  ${LAPACK} ${ARPACK} ${INCLUDE} 

matrix_stuff.o: matrix_stuff.f
	${CMP} ${WARN} ${F132FORM} ${OPTFLAG} -c matrix_stuff.f

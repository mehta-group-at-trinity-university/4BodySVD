CMP     = ifort
OPTFLAG = -O2 -parallel
F132FORM = -extend_source
LAPACK =   -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lm
ARPACK =  -L/opt/ARPACK/ -larpack_Intel 
INCLUDE =  -I/usr/local/include

4BodySVD.x:	4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90
	${CMP} ${OPT} ${F132FORM} ${INCLUDE}  4BodySVD.f90 adiabaticSolver.f90 Bsplines.f90 FourBodyPOT.f90 matrix_stuff.f90 ${ARPACK} ${LAPACK} 


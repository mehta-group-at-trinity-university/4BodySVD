FourBody1D.x:	FourBody1D.f FourBody1D_Bsplines.f FourBody1D_matrix_stuff.f FourBody1D.o FourBody1D_Bsplines.o FourBody1D_matrix_stuff.o FourBodyPOT.o
	ifort -O4 -CB -debug extended FourBody1D.o FourBody1D_matrix_stuff.o FourBody1D_Bsplines.o FourBodyPOT.o -L/usr/local/lib -L/usr/local/intel/mkl/lib/em64t -larpack -L/usr/lib64 -lguide -lmkl_solver -lmkl_lapack -lmkl_em64t -lpthread -Vaxlib -L/usr/local/cernlib/lib/ -o FourBody1D.x

FourBodyPOT.o:	FourBodyPOT.f
	ifort -O4 -CB -c FourBodyPOT.f

FourBody1D.o:	FourBody1D.f
	ifort -O4 -CB -extend_source -debug extended -c FourBody1D.f

FourBody1D_matrix_stuff.o:	FourBody1D_matrix_stuff.f
	ifort -O4 -CB -extend_source -L/usr/local/lib/ -larpack -c FourBody1D_matrix_stuff.f

FourBody1D_Bsplines.o:	FourBody1D_Bsplines.f
	ifort -O4 -CB -extend_source -c FourBody1D_Bsplines.f

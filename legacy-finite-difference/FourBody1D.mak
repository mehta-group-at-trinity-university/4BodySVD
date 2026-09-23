FourBody1D.x:	FourBody1D.f FourBody1D_Bsplines.f FourBody1D_matrix_stuff.f FourBody1D.o FourBody1D_Bsplines.o FourBody1D_matrix_stuff.o FourBodyPOT.o
	ifort -O4 -CB -debug extended FourBody1D.o FourBody1D_matrix_stuff.o FourBody1D_Bsplines.o FourBodyPOT.o -i8 -I$(MKLROOT)/include/intel64/ilp64 -I$(MKLROOT)/include $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -L/opt/ARPACK/ -larpack_Intel -o FourBody1D.x

FourBodyPOT.o:	FourBodyPOT.f
	ifort -O4 -CB -c FourBodyPOT.f

FourBody1D.o:	FourBody1D.f
	ifort -O4 -CB -extend_source -c FourBody1D.f

FourBody1D_matrix_stuff.o:	FourBody1D_matrix_stuff.f
	ifort -O4 -CB -extend_source -L/opt/ARPACK/ -larpack_Intel -c FourBody1D_matrix_stuff.f

FourBody1D_Bsplines.o:	FourBody1D_Bsplines.f
	ifort -O4 -CB -extend_source -c FourBody1D_Bsplines.f

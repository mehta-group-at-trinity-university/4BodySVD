      
      subroutine sumpairwisepot(r12, r13, r14, r23, r24, r34, 
     >     potvalue) 
      
c     returns the value of the sum of pair-wise interactions 
      
      implicit none
      real *8 c0, cutoff, potvalue, r12, r13, r14, r23, r24, r34,dd,L
      
c      common /vindex/ index(67,3)
      
c      dimension vcof(67)
      
      potvalue=0.d0

c      c0=-8.694192e+00 ! a2=100
c      c0=-8.605163e+00 ! a2=-100
c      c0=-7.321700e+00 ! a2=-2
c      c0=-8.887028e+00 ! a2=20
c      cutoff=1.d0
c      potvalue=c0*cutoff*(
c     >     dexp(-(r12*cutoff)**2.d0) +
c     >     dexp(-(r13*cutoff)**2.d0) + 
c     >     dexp(-(r14*cutoff)**2.d0) + 
c     >     dexp(-(r23*cutoff)**2.d0) + 
c     >     dexp(-(r24*cutoff)**2.d0) + 
c     >     dexp(-(r34*cutoff)**2.d0))

      L=1.d0
c      dd=6.272844447382345d0    ! a2even=20, with one deep bound state L=1
c      dd=5.64556d0  !bound state at E0 = -3.72 and E1 = -0.86
c      dd=6.0d0 ! a2even=Infinity with one deep bound state L=1
c      dd=24.0d0 ! a2even=Infinity with one deep bound state L=1/2

      dd=0.853138729379683 ! a2even=2 B2=0.3028346165725004  L=1

!      dd=0.5d0 ! a2even=2, B2=0.25 L=2
!      L=2.0d0


c      dd=0.05260553185042266d0 ! a2even=20  B2=0.0025096021704961394  L=1
c      L=1.0d0
c      dd=0.11086487977899903d0 ! a2even=10 L=1 B2=0.010144578990854353
      

      potvalue=-dd*(
     >     1.001*dcosh(r12/L)**(-2.0d0) + 
     >     dcosh(r13/L)**(-2.0d0) + 
     >     dcosh(r14/L)**(-2.0d0) + 
     >     dcosh(r23/L)**(-2.0d0) + 
     >     dcosh(r24/L)**(-2.0d0) + 
     >     dcosh(r34/L)**(-2.0d0)) 
      return
      end 

      

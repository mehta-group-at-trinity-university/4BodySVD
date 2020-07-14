
      
      subroutine sumpairwisepot(r12, r13, r14, r23, r24, r34, V2Depth, potvalue,rvec,vvec) 
      
!     returns the value of the sum of pair-wise interactions 
      
      implicit none
      integer run
      double precision c0, cutoff, potvalue, r12, r13, r14, r23, r24, r34,dd,r0,a,L
!      double precision v1,v2,v3,v4,v5,v6
      double precision V2Depth
      double precision rvec(161),vvec(161)

!      common /vindex/ index(67,3)
      
!      dimension vcof(67)
      
      potvalue=0.d0

!      c0=-8.694192e+00 ! a2=100
!      c0=-8.605163e+00 ! a2=-100
!      c0=-7.321700e+00 ! a2=-2
!      c0=-8.887028e+00 ! a2=20
!      cutoff=1.d0
!      potvalue=c0*cutoff*(
!     >     dexp(-(r12*cutoff)**2.d0) +
!     >     dexp(-(r13*cutoff)**2.d0) + 
!     >     dexp(-(r14*cutoff)**2.d0) + 
!     >     dexp(-(r23*cutoff)**2.d0) + 
!     >     dexp(-(r24*cutoff)**2.d0) + 
!     >     dexp(-(r34*cutoff)**2.d0))

!      dd=1.0d0 ! a2even=20, with one deep bound state L=1

!      dd = 0.853138729379683d0
!      dd=dd*0.9d0
!      dd=6.0d0 ! a2even=Infinity with one deep bound state L=1
!      dd=24.0d0 ! a2even=Infinity with one deep bound state L=1/2

!      dd=0.853138729379683 ! a2even=2 B2=0.3028346165725004  L=1

!      dd=0.5d0 ! a2even=2, B2=0.25 L=2
!      L=2.0d0

      !write(6,*) "V2Depth", V2Depth, d
!      dd=0.05260553185042266d0 ! a2even=20  B2=0.0025096021704961394  L=1
!      L=1.0d0
!      dd=0.11086487977899903d0 ! a2even=10 L=1 B2=0.010144578990854353
      
!
! Morse Potential
      dd=V2Depth
      r0=1.0d0
      a=1.0d0
      potvalue = dd*((1-exp(-a*(r12-r0)))**(2.0d0)) - dd + &
          dd*((1d0 - exp(-a*(r13-r0)))**(2.0d0)) - dd + &
          dd*((1d0 - exp(-a*(r34-r0)))**(2.0d0)) - dd + &
          dd*((1d0 - exp(-a*(r23-r0)))**(2.0d0)) - dd + &
          dd*((1d0 - exp(-a*(r24-r0)))**(2.0d0)) - dd + &
          dd*((1d0 - exp(-a*(r14-r0)))**(2.0d0)) - dd

     !      dd=V2Depth
!      potvalue=-dd*(dcosh(r12/L)**(-2.0d0) + 
!     >     dcosh(r13/L)**(-2.0d0) + 
!     >     dcosh(r34/L)**(-2.0d0) + 
!     >     dcosh(r23/L)**(-2.0d0) + 
!     >     dcosh(r24/L)**(-2.0d0) + 
!     >     dcosh(r14/L)**(-2.0d0)) 

!      potvalue = 0.0d0
      !      write(6,*) "pot= ", potvalue
!      potvalue=-dd/(1+(r12/L)**(6.0d0)) -
!     >     dd/(1+(r13/L)**(6.0d0)) - 
!     >     dd/(1+(r14/L)**(6.0d0)) -
!     >     dd/(1+(r23/L)**(6.0d0)) -
!     >     dd/(1+(r24/L)**(6.0d0)) -
!     >     dd/(1+(r34/L)**(6.0d0))


!cc new potential ccc
!!$      call singlePot(0,vvec,r12,v1,run)
!!$      call singlePot(1,vvec,r13,v2,run)
!!$      call singlePot(0,vvec,r34,v3,run)
!!$      call singlePot(1,vvec,r23,v4,run)
!!$      call singlePot(1,vvec,r24,v5,run)
!!$      call singlePot(1,vvec,r14,v6,run)
!!$      potvalue = v1+v2+v3+v4+v5+v6

!      potvalue=4.0d0*4000.0d0*((4.0d0/r12)**12.0d0-(4.0d0/r12)**6.0d0+
!     >      (4.0d0/r13)**12.0d0-(4.0d0/r13)**6.0d0+
!     >      (4.0d0/r23)**12.0d0-(4.0d0/r23)**6.0d0+
!     >      (4.0d0/r14)**12.0d0-(4.0d0/r14)**6.0d0+
!     >      (4.0d0/r24)**12.0d0-(4.0d0/r24)**6.0d0+
!     >      (4.0d0/r34)**12.0d0-(4.0d0/r34)**6.0d0)


!      write(6,*) potvalue
      return
    end subroutine sumpairwisepot

      subroutine singlePot(rflag,vvec,R,potval,run)
      integer rflag,run
      double precision R, potval, depth
      double precision A,B
      double precision a0,littleb,Rm,C6,C8,C10,Aex,gamma,beta
      double precision rvec(161),vvec(161)
      double precision, allocatable :: avec(:)
      allocate(avec(36))
      

      ! Realistic atom-atom potential (not sure what atom!  Have to ask Kelly or Kaden, or look up the C6 to determine the species.
      A = -0.731001924d4;
      B = 0.590300435d6;
      avec = (/ 0.5530165860346669d1, 0.1402690940409950d5, &
             0.1048318435177453d5, -0.4087971898030580d4, &
             -0.1817659880541463d5, -0.2739889053154796d5,&
             -0.3930730951201573d5, 0.1755013192862486d6, &
             0.4366867152340260d6, -0.6031147948206772d7, &
              -0.7948028201366105d7, 0.1172293247025888d9, &
             0.9209455256556711d8, -0.1533597323741656d10,&
             -0.6982526182585652d9, 0.1418207140035536d11, &
             0.3514113318792982d10, -0.9529625435553510d11, &
             -0.1096758619400216d11, 0.4743826569413903d12, &
             0.1477139521952726d11, -0.1768149926253122d13, &
             0.3496665669778944d11, 0.4947896972615874d13, &
             -0.2373506540656795d12, -0.1033913349509970d14,& 
              0.6234843170080333d12, 0.1588626857834081d14, &
              -0.9740574795906707d12, -0.1741861156434701d14, &
             0.9410306699092036d12, 0.1289614354386887d14, &
             -0.5220986485346664d12, -0.5776754651606419d13,&
             0.1278431437418345d12, 0.1182411100178171d13 /)

      a0 = -4217.814754d0;
      littleb = -0.39d0;
      Rm = 4.06818602d0;
      C6 = 0.2072097d8;
      C8 = 0.6509487d9;
      C10 = 0.2575245d11;
      Aex = 0.1469387d5;
      gamma = 5.25669;
      beta = 2.11445;
      
      potval=0.0d0

      if (abs(R) .le. 1.5d0) then
          potval =109293.d0
          else if (abs(R) .le. 3.0d0) then
          potval = A + B/(R**4)
      else if (abs(R) .le. 11.0d0) then
          potval = a0
          do i = 1,36
              potval = potval+avec(i)*((R - Rm)/(R + b*Rm))**i
          enddo  
      else
          potval=-Aex*(R**gamma)*exp(-beta*R)-C6/R**6-C8/R**8-C10/R**10
          
       endif

       !6-12 potential
!      if (abs(R).le.3.7d0) then
!          potval=152.345d0
!      else
!          potval=4.0d0*40.0d0*((4.0d0/abs(R))**12.0d0-(4.0d0/abs(R))**6.0d0)
!      endif

!      if (rflag.eq.1d0) then
        !    potval=-20.0d0*(cosh(abs(R)/0.5d0))**(-2.0d0)
        !if (run==1) depth=50d0
        !if (run==2) depth=50d0
        !if (run==3) depth=100d0
        !if (run==4) depth=100d0
        !if (run==5) depth=300d0
      !depth=20d0
      !potval = depth*(((1-exp(-1.0d0/1.0d0*(abs(R)-1.0d0)))**2.0d0)-1.0d0)  
!      else
!          if (R.le.0.01d0) then
!              potval=10000.d0
!          else
!              potval=0.d0
!          endif
!      endif

!      potval = 20.0d0*(((1-exp(-1.0d0/1.0d0*(abs(R)-1.0d0)))**2.0d0)-1.0d0)
      
      return
    end subroutine singlePot
!!$
!!$     
!!$      
!!$!***********************************************************************
!!$!                                                                      *
!!$!     SUBROUTINE AKIMA (L,X,Y,N,U,V)                                   *
!!$!     INTERPOLATION OF A SINGLE-VALUED FUNCTION                        *
!!$!                                                                      *
!!$!     L = NUMBER OF INPUT DATA POINTS                                  *
!!$!          (MUST BE 2 OR GREATER)                                      *
!!$!     X = ARRAY OF DIMENSION L STORING THE X VALUES                    *
!!$!         (ABSCISSAS) OF INPUT DATA POINTS                             *
!!$!          (IN ASCENDING ORDER)                                        *
!!$!     Y = ARRAY OF DIMENSION L STORING THE Y VALUES                    *
!!$!         (ORDINATES) OF INPUT DATA POINTS                             *
!!$!     N = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE               *
!!$!         Y VALUE (ORDINATE) IS DESIRED                                *
!!$!          (MUST BE 1 OR GREATER)                                      *
!!$!     U = ARRAY OF DIMENSION N STRORING THE X VALUES                   *
!!$!         (ABSCISSAS) OF THE DESIRED POINTS                            *
!!$!     V = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y                *
!!$!         VALUES (ORDINATES) ARE TO BE DISPLAYED                       *
!!$!                                                                      *
!!$!***********************************************************************
!!$!
!!$      SUBROUTINE AKIMA (L,X,Y,N,U,V)
!!$      IMPLICIT REAL*8 (A-H,O-Z)
!!$      DIMENSION X(L),Y(L),U(N),V(N)
!!$      EQUIVALENCE (P0,X3),(Q0,Y3),(Q1,T3)
!!$      REAL*8 M1,M2,M3,M4,M5
!!$      EQUIVALENCE (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),(J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
!!$ 10   L0=L
!!$      LM1=L0-1
!!$      LM2=LM1-1
!!$      LP1=L0+1
!!$      N0=N
!!$      IF (LM2.LT.0) GO TO 90
!!$      IF (N0.LE.0) GO TO 91
!!$      DO 11 I=2,L0
!!$      IF (X(I-1)-X(I)) 11,95,96
!!$ 11   CONTINUE
!!$      IPV=0
!!$      DO 80 K=1,N0
!!$      UK=U(K)
!!$ 20   IF (LM2.EQ.0) GO TO 27
!!$      IF (UK.GE.X(L0)) GO TO 26
!!$      IF (UK.LT.X(1)) GO TO 25
!!$      IMN=2
!!$      IMX=L0
!!$ 21   I=(IMN+IMX)/2
!!$      IF (UK.GE.X(I)) GO TO 23
!!$ 22   IMX=I
!!$      GO TO 24
!!$ 23   IMN=I+1
!!$ 24   IF (IMX.GT.IMN) GO TO 21
!!$      I=IMX
!!$      GO TO 30
!!$ 25   I=1
!!$      GO TO 30
!!$ 26   I=LP1
!!$      GO TO 30
!!$ 27   I=2
!!$ 30   IF (I.EQ.IPV) GO TO 70
!!$      IPV=I
!!$ 40   J=I
!!$      IF (J.EQ.1) J=2
!!$      IF (J.EQ.LP1) J=L0
!!$      X3=X(J-1)
!!$      Y3=Y(J-1)
!!$      X4=X(J)
!!$      Y4=Y(J)
!!$      A3=X4-X3
!!$      M3=(Y4-Y3)/A3
!!$      IF (LM2.EQ.0) GO TO 43
!!$      IF (J.EQ.2) GO TO 41
!!$      X2=X(J-2)
!!$      Y2=Y(J-2)
!!$      A2=X3-X2
!!$      M2=(Y3-Y2)/A2
!!$      IF (J.EQ.L0) GO TO 42
!!$ 41   X5=X(J+1)
!!$      Y5=Y(J+1)
!!$      A4=X5-X4
!!$      M4=(Y5-Y4)/A4
!!$      IF (J.EQ.2) M2=M3+M3-M4
!!$      GO TO 45
!!$ 42   M4=M3+M3-M2
!!$      GO TO 45
!!$ 43   M2=M3
!!$      M4=M3
!!$ 45   IF (J.LE.3) GO TO 46
!!$      A1=X2-X(J-3)
!!$      M1=(Y2-Y(J-3))/A1
!!$      GO TO 47
!!$ 46   M1=M2+M2-M3
!!$ 47   IF (J.GE.LM1) GO TO 48
!!$      A5=X(J+2)-X5
!!$      M5=(Y(J+2)-Y5)/A5
!!$      GO TO 50
!!$ 48   M5=M4+M4-M3
!!$ 50   IF (I.EQ.LP1) GO TO 52
!!$      W2=DABS(M4-M3)
!!$      W3=DABS(M2-M1)
!!$      SW=W2+W3
!!$      IF (SW.NE.0.D0) GO TO 51
!!$      W2=.5D0
!!$      W3=.5D0
!!$      SW=1.D0
!!$ 51   T3=(W2*M2+W3*M3)/SW
!!$      IF (I.EQ.1) GO TO 54
!!$ 52   W3=DABS(M5-M4)
!!$      W4=DABS(M3-M2)
!!$      SW=W3+W4
!!$      IF (SW.NE.0.D0) GO TO 53
!!$      W3=.5D0
!!$      W4=.5D0
!!$      SW=1.D0
!!$ 53   T4=(W3*M3+W4*M4)/SW
!!$      IF (I.NE.LP1) GO TO 60
!!$      T3=T4
!!$      SA=A2+A3
!!$      T4=.5D0*(M4+M5-A2*(A2-A3)*(M2-M3)/(SA*SA))
!!$      X3=X4
!!$      Y3=Y4
!!$      A3=A2
!!$      M3=M4
!!$      GO TO 60
!!$ 54   T4=T3
!!$      SA=A3+A4
!!$      T3=.5D0*(M1+M2-A4*(A3-A4)*(M3-M4)/(SA*SA))
!!$      X3=X3-A4
!!$      Y3=Y3-M2*A4
!!$      A3=A4
!!$      M3=M2
!!$ 60   Q2=(2.D0*(M3-T3)+M3-T4)/A3
!!$      Q3=(-M3-M3+T3+T4)/(A3*A3)
!!$ 70   DX=UK-P0
!!$ 80   V(K)=Q0+DX*(Q1+DX*(Q2+DX*Q3))
!!$      RETURN
!!$ 90   WRITE (6,2090)
!!$      GO TO 99
!!$ 91   WRITE (6,2091)
!!$      GO TO 99
!!$ 95   WRITE (6,2095)
!!$      GO TO 97
!!$ 96   WRITE (6,2096)
!!$ 97   WRITE (6,2097) I,X(I)
!!$ 99   WRITE (6,2099) L0,N0
!!$      RETURN
!!$ 2090 FORMAT ('***   L = 1 OR LESS.')
!!$ 2091 FORMAT (' ***   N = 0 OR LESS.')
!!$ 2095 FORMAT (' ***   IDENTICAL X VALUES.')
!!$ 2096 FORMAT (' ***   X VALUES OUT OF SEQUENCE.')
!!$ 2097 FORMAT ('   I =',I7,10X,'X(I) =',D12.3)
!!$ 2099 FORMAT ('   L =',I7,10X,'N =',I7/'ERROR DETECTED IN ROUTINE     AKIMA')
!!$      END

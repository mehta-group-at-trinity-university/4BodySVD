      program FourBody

      integer LegPoints,xNumPoints,yNumPoints
      integer NumStates,PsiFlag,Order,Left,Right,Bottom,Top
      integer RSteps,CouplingFlag
      double precision alpha,m,Shift
      double precision RLeft,RRight,RDerivDelt
      DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX
      double precision xMin,xMax,yMin,yMax
      double precision, allocatable :: R(:)
      double precision, allocatable :: xPoints(:),yPoints(:)

      logical, allocatable :: Select(:)

      integer iparam(11),ncv,info
      integer i,j,k,iR
      integer LeadDim,MatrixDim,HalfBandWidth
      integer xDim,yDim
      integer, allocatable :: iwork(:)
      integer, allocatable :: xBounds(:),yBounds(:)
      double precision Tol
      double precision TotalMemory
      double precision xNUM,dNUM
      double precision mu, mu12, r0diatom, dDiatom, etaOVERpi, Pi
      double precision u1,v1,sys_ss_pot
      double precision, allocatable :: LUFac(:,:),workl(:)
      double precision, allocatable :: workd(:),Residuals(:)
      double precision, allocatable :: xLeg(:),wLeg(:)
      double precision, allocatable :: u(:,:,:),uxx(:,:,:),v(:,:,:)
      double precision, allocatable :: vy(:,:,:),vyy(:,:,:)
      double precision, allocatable :: S(:,:),H(:,:)
      double precision, allocatable :: lPsi(:,:),mPsi(:,:),rPsi(:,:),
     >     Energies(:,:)
      double precision, allocatable :: P(:,:),Q(:,:),dP(:,:)

      common/MassInfo/mu12,r0diatom,dDiatom

      character*64 LegendreFile

c     read in number of energies and states to print
      read(5,*)
      read(5,*) NumStates,PsiFlag,CouplingFlag
      write(6,*) NumStates,PsiFlag,CouplingFlag

c     read in Gauss-Legendre info
      read(5,*)
      read(5,*)
      read(5,1002) LegendreFile
      write(6,1002) LegendreFile
      read(5,*)
      read(5,*)
      read(5,*) LegPoints
      write(6,*) LegPoints,' LegPoints'

c     read in boundary conditions
      read(5,*)
      read(5,*)
      read(5,*) Shift,Order,Left,Right,Bottom,Top
      write(6,*) Shift,Order,Left,Right,Bottom,Top

c     read in potential parameters
      read(5,*)
      read(5,*)
      read(5,*) alpha,m, r0diatom, xNUM
      write(6,*) alpha,m,r0diatom,xNUM
      dNUM=2.d0*xNUM*xNUM-xNUM
c     mu12=m/2.d0
c     mu12=40.d0*m/2.d0
      Pi=dacos(-1.d0)
      write(6,*) 'Pi=',Pi
c     dDiatom=(dNUM)/(mu12*r0diatom**2)
c     write(6,*) 'dDiatom=',dDiatom,' xNUM=',xNUM,' dNUM',dNUM

c     read in grid information
      read(5,*)
      read(5,*)
      read(5,*) xNumPoints,xMin,xMax
      write(6,*) xNumPoints,xMin,xMax

      read(5,*)
      read(5,*)
      read(5,*) yNumPoints,yMin,yMax
      write(6,*) yNumPoints,yMin,yMax

      read(5,*)
      read(5,*)
      read(5,*) RSteps,RDerivDelt,RFirst,RLast
      write(6,*) RSteps,RDerivDelt,RFirst,RLast
c     c	XFirst=dsqrt(RFirst)
c     c	XLast=dsqrt(RLast)
      XFirst=RFirst**(1.d0/3.d0)
      XLast=RLast**(1.d0/3.d0)
c     XFirst = dlog10(RFirst)
c      XLast = dlog10(RLast)
      StepX=(XLast-XFirst)/(RSteps-1.d0)
      
      allocate(R(RSteps))
      do i = 1,RSteps
c     read(5,*) R(i)
         R(i)= (XFirst+(i-1)*StepX)**3
c     R(i)= 10.d0**(XFirst+(i-1)*StepX)
c         R(i) = RFirst + (i-1)*(RLast - RFirst)/dble(RSteps-1) ! linear radial grid
      enddo

      if (mod(xNumPoints,2) .ne. 0) then
         write(6,*) 'xNumPoints not divisible by 2'
         xNumPoints = (xNumPoints/2)*2
         write(6,*) '   truncated to ',xNumPoints
      endif

      if (mod(yNumPoints,2) .ne. 0) then
         write(6,*) 'yNumPoints not divisible by 2'
         yNumPoints = (yNumPoints/2)*2
         write(6,*) '   truncated to ',yNumPoints
      endif

      RLeft = 0.0d0
c     u1 = sys_ss_pot(RLeft,v1,2,.TRUE.)
c      u1 = 0.1d0*dexp(Rleft)
c     mu = m/dsqrt(3.0d0)
c     mu = m/(4.d0**(1.d0/3.d0))
c     mu=0.5
      mu=1.d0
c     mu=2.d0*mu12
      allocate(xLeg(LegPoints),wLeg(LegPoints))
      call GetGaussFactors(LegendreFile,LegPoints,xLeg,wLeg)

      xDim = xNumPoints+Order-3
      if (Left .eq. 2) xDim = xDim + 1
      if (Right .eq. 2) xDim = xDim + 1
      yDim = yNumPoints+Order-3
      if (Top .eq. 2) yDim = yDim + 1
      if (Bottom .eq. 2) yDim = yDim + 1

      MatrixDim = xDim*yDim
      HalfBandWidth = yDim*Order+Order

      TotalMemory = (4*(HalfBandWidth+1)+(3*HalfBandWidth+1)+
     >     6*NumStates)*8.0d0*MatrixDim
      TotalMemory = TotalMemory/(1024.0d0*1024.0d0)
      
      write(6,*)
      write(6,*) 'MatrixDim ',MatrixDim
      write(6,*) 'HalfBandWidth ',HalfBandWidth
      write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
      write(6,*)

      allocate(xPoints(xNumPoints),yPoints(yNumPoints))
      allocate(xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order))
      allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,
     >     xDim))
      allocate(v(LegPoints,yNumPoints,yDim),vy(LegPoints,yNumPoints,
     >     yDim),vyy(LegPoints,yNumPoints,yDim))
      allocate(S(HalfBandWidth+1,MatrixDim),H(HalfBandWidth+1,
     >     MatrixDim))
      allocate(P(NumStates,NumStates),Q(NumStates,NumStates),
     >     dP(NumStates,NumStates))
      
      ncv = 2*NumStates
      LeadDim = 3*HalfBandWidth+1
      allocate(iwork(MatrixDim))
      allocate(Select(ncv))
      allocate(LUFac(LeadDim,MatrixDim))
      allocate(workl(ncv*ncv+8*ncv))
      allocate(workd(3*MatrixDim))
      allocate(lPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv),
     >     rPsi(MatrixDim,ncv))
      allocate(Residuals(MatrixDim))
      allocate(Energies(ncv,2))
      info=0
         call GridMaker(m,mu,R(iR),11.65d0,xNumPoints,xMin,xMax,
     >        yNumPoints,yMin,yMax,xPoints,yPoints)

         call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >        xDim,xBounds,xNumPoints,0,u)
         call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,
     >        xDim,xBounds,xNumPoints,2,uxx)
         call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,
     >        yDim,yBounds,yNumPoints,0,v)
         call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,
     >        yDim,yBounds,yNumPoints,1,vy)
         call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,
     >        yDim,yBounds,yNumPoints,2,vyy)
         
         call CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,
     >        xDim,yDim,
     >        xNumPoints,yNumPoints,u,v,
     >        xBounds,yBounds,HalfBandWidth,S)

      do iR = 1,RSteps


         if (CouplingFlag .ne. 0) then

            RLeft = R(iR)-RDerivDelt
            call CalcHamiltonian(alpha,RLeft,mu,Order,xPoints,yPoints,
     >           LegPoints,xLeg,wLeg,xDim,yDim,
     >           xNumPoints,yNumPoints,u,v,vy,
     >           uxx,vyy,xBounds,yBounds,HalfBandWidth,H)
            call MyDsband(Select,Energies,lPsi,MatrixDim,Shift,MatrixDim,
     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
     >           NumStates,Tol,Residuals,ncv,lPsi,MatrixDim,iparam,workd,workl,
     >           ncv*ncv+8*ncv,iwork,info)
            if (info.ne.0) write(6,*) 'Error in MyDsband.  info = ',info
            if (iR .gt. 1) then 
               call FixPhase(NumStates,HalfBandWidth,
     >              MatrixDim,S,ncv,mPsi,lPsi)
!               write(6,*) 'Finished with FixPhase at RLeft.'
            endif

            call CalcEigenErrors(info,iparam,MatrixDim,H,
     >           HalfBandWidth+1,
     >           S,HalfBandWidth,NumStates,lPsi,Energies,ncv)
!            IF(R(iR).GT. 2.d0) Shift = 1.05d0*Energies(1,1)
!            do i = 1,min(NumStates,iparam(5))
c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*RLeft*RLeft)
!            enddo
            write(100,10) RLeft,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) RLeft
            do i = 1,min(NumStates,iparam(5))
               write(6,*) i,Energies(i,1),Energies(i,2)
            enddo

            RRight = R(iR)+RDerivDelt
            call CalcHamiltonian(alpha,RRight,mu,Order,xPoints,yPoints,
     >           LegPoints,xLeg,wLeg,xDim,yDim,
     >           xNumPoints,yNumPoints,u,v,vy,uxx,vyy,xBounds,yBounds,
     >           HalfBandWidth,H)
            call MyDsband(Select,Energies,rPsi,MatrixDim,Shift,
     >           MatrixDim,
     >           H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,
     >           NumStates,Tol,
     >           Residuals,ncv,rPsi,MatrixDim,iparam,workd,workl,
     >           ncv*ncv+8*ncv,iwork,info)
            if (info.ne.0) write(6,*) 'Error in MyDsband.  info = ',info
            call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,
     >           lPsi,rPsi)
!            write(6,*) 'Finished with FixPhase at RRight.'
            call CalcEigenErrors(info,iparam,MatrixDim,H,
     >           HalfBandWidth+1,
     >           S,HalfBandWidth,NumStates,rPsi,Energies,ncv)
!            IF(R(iR).GT. 2.2d0) Shift = 0.93d0*Energies(1,1)
!            do i = 1,min(NumStates,iparam(5))
c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*RRight*RRight)
!            enddo
            write(100,10) RRight,(Energies(i,1), i = 1,NumStates)
            write(6,*)
            write(6,*) RRight
            do i = 1,min(NumStates,iparam(5))
               write(6,*) i,Energies(i,1),Energies(i,2)
            enddo

         endif
         
         call CalcHamiltonian(alpha,R(iR),mu,Order,xPoints,yPoints,
     >        LegPoints,xLeg,wLeg,xDim,yDim,
     >        xNumPoints,yNumPoints,u,v,vy,uxx,vyy,xBounds,yBounds,
     >        HalfBandWidth,H)
         call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,
     >        H,S,HalfBandWidth+1,LUFac,LeadDim,HalfBandWidth,NumStates,
     >        Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,
     >        ncv*ncv+8*ncv,iwork,info)
         if (info.ne.0) write(6,*) 'Error in MyDsband.  info = ',info
         
         if (CouplingFlag .ne. 0) then 
            call FixPhase(NumStates,HalfBandWidth,
     >           MatrixDim,S,ncv,rPsi,mPsi)
!            write(6,*) 'Finished with FixPhase at RRight.'
         endif

         call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,
     >        HalfBandWidth,NumStates,mPsi,Energies,ncv)
!         IF(R(iR).GT. 2.2d0) Shift = 0.93d0*Energies(1,1)
         do i = 1,min(NumStates,iparam(5))

c     Energies(i,1) = Energies(i,1) + 1.875d0/(mu*R(iR)*R(iR))

         enddo
         write(200,20) R(iR),(Energies(i,1), i = 1,min(NumStates,
     >        iparam(5)))
         write(6,*)
         write(6,*) R(iR)
         do i = 1,min(NumStates,iparam(5))
            write(6,*) i,Energies(i,1),Energies(i,2)
         enddo

         if (CouplingFlag .ne. 0) then
c     call CalcPMatrix(min(NumStates,iparam(5)),HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,P)
c     call CalcQMatrix(min(NumStates,iparam(5)),HalfBandWidth,MatrixDim,RDerivDelt,lPsi,mPsi,rPsi,S,Q)
            call CalcCoupling(NumStates,HalfBandWidth,MatrixDim,
     >           RDerivDelt,lPsi,mPsi,rPsi,S,P,Q,dP)
            
            write(101,*) R(iR)
            write(102,*) R(iR)
            write(103,*) R(iR)
            do i = 1,min(NumStates,iparam(5))
               write(101,20) (P(i,j), j = 1,min(NumStates,iparam(5)))
               write(102,20) (Q(i,j), j = 1,min(NumStates,iparam(5)))
               write(103,20) (dP(i,j), j = 1,min(NumStates,iparam(5)))
            enddo
         endif
         
         if (PsiFlag .ne. 0) then
            do i = 1,xNumPoints
               write(97,*) xPoints(i)
            enddo
            do i = 1,yNumPoints
               write(98,*) yPoints(i)
            enddo
            do i = 1,MatrixDim
               write(999+iR,20) (mPsi(i,j), j = 1,NumStates)
            enddo
            close(unit=999+iR)
         endif

      enddo

      deallocate(S,H)
      deallocate(Energies)
      deallocate(iwork)
      deallocate(Select)
      deallocate(LUFac)
      deallocate(workl)
      deallocate(workd)
      deallocate(lPsi,mPsi,rPsi)
      deallocate(Residuals)
      deallocate(P,Q,dP)
      deallocate(xPoints,yPoints)
      deallocate(xLeg,wLeg)
      deallocate(xBounds,yBounds)
      deallocate(u,uxx)
      deallocate(v,vy,vyy)
      deallocate(R)

 10   format(1P,100e25.15)
 20   format(1P,100e16.8)
 1002 format(a64)

      stop
      end

      subroutine CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,
     >     xNumPoints,yNumPoints,u,v,xBounds,yBounds,HalfBandWidth,S)
      implicit none
      integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints,
     >     xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order),HalfBandWidth
      double precision xPoints(*),yPoints(*),xLeg(*),wLeg(*)
      double precision S(HalfBandWidth+1,xDim*yDim)
      double precision u(LegPoints,xNumPoints,xDim)
      double precision v(LegPoints,yNumPoints,yDim)

      integer ix,iy,ixp,iyp,kx,ky,lx,ly
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)
      double precision a,b,m
      double precision xTempS
      double precision ax,bx
      double precision y,ay,by,yScaledZero,yTempS
      double precision, allocatable :: xIntScale(:),xS(:,:)
      double precision, allocatable :: siny(:,:),yIntScale(:),yS(:,:)

      allocate(xIntScale(xNumPoints),xS(xDim,xDim))
      allocate(siny(LegPoints,yNumPoints),yIntScale(yNumPoints),yS(yDim,yDim))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(kyMin(yDim,yDim),kyMax(yDim,yDim))

      S = 0.0d0

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
      enddo
      do ky = 1,yNumPoints-1
         ay = yPoints(ky)
         by = yPoints(ky+1)
         yIntScale(ky) = 0.5d0*(by-ay)
         yScaledZero = 0.5d0*(by+ay)
         do ly = 1,LegPoints
            y = yIntScale(ky)*xLeg(ly)+yScaledZero
c     siny(ly,ky) = dsin(4.0d0*y)
            siny(ly,ky) = dsin(y)
         enddo
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo
      do iy = 1,yDim
         do iyp = 1,yDim
            kyMin(iyp,iy) = max(yBounds(iy),yBounds(iyp))
            kyMax(iyp,iy) = min(yBounds(iy+Order+1),yBounds(iyp+Order+1))-1
         enddo
      enddo

      do ix = 1,xDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            xS(ixp,ix) = 0.0d0
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               xTempS = 0.0d0
               do lx = 1,LegPoints
                  a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                  b = a*u(lx,kx,ixp)
                  xTempS = xTempS + b
               enddo
               xS(ixp,ix) = xS(ixp,ix) +   xTempS
            enddo
         enddo
      enddo

      do iy = 1,yDim
         do iyp = max(1,iy-Order),min(yDim,iy+Order)
            yS(iyp,iy) = 0.0d0
            do ky = kyMin(iyp,iy),kyMax(iyp,iy)
               yTempS = 0.0d0
               do ly = 1,LegPoints
                  a = wLeg(ly)*yIntScale(ky)*v(ly,ky,iy)
                  b = a*v(ly,ky,iyp)
                  yTempS = yTempS + b*siny(ly,ky)
               enddo
               yS(iyp,iy) = yS(iyp,iy) +   yTempS
            enddo
         enddo
      enddo

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            do iy = 1,yDim
               Row = i1+iy
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  Col = i1p+iyp
                  if (Col .ge. Row) then
                     NewRow = HalfBandWidth+1+Row-Col
                     S(NewRow,Col) = xS(ixp,ix)*yS(iyp,iy)
                  endif
               enddo
            enddo
         enddo
      enddo

      deallocate(xIntScale,xS)
      deallocate(siny,yIntScale,yS)
      deallocate(kxMin,kxMax)
      deallocate(kyMin,kyMax)

      return
      end

      subroutine CalcHamiltonian(alpha,R,mu,
     >     Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,
     >     xNumPoints,yNumPoints,u,v,vy,uxx,vyy,xBounds,yBounds,HalfBandWidth,H)

      integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints,xBounds(*),yBounds(*),HalfBandWidth
      double precision alpha,R,mu
      double precision xPoints(*),yPoints(*),xLeg(*),wLeg(*)
      double precision H(HalfBandWidth+1,xDim*yDim)
      double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
      double precision v(LegPoints,yNumPoints,yDim),vy(LegPoints,yNumPoints,yDim),vyy(LegPoints,yNumPoints,yDim)

      integer ix,iy,ixp,iyp,kx,ky,lx,ly
      integer i1,i1p
      integer Row,NewRow,Col
      integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)
      double precision a,b,m,Pi
      double precision Rall,r12,r12a,r23,r23a,r23b,r23c,r13,r13a,r13b,r13c,r14,r24,r34
      double precision u1,sys_ss_pot,V12,V23,V31
      double precision VInt,VTempInt,potvalue
c     double precision TempPot,VInt,VTempInt
      double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
      double precision y,ay,by,yScaledZero,yTempT,yTempF,yInt
      double precision, allocatable :: Pot(:,:,:,:)
      double precision, allocatable :: xIntScale(:),xTempV(:),xS(:,:),xT(:,:)
      double precision, allocatable :: siny(:,:),cosy(:,:),tany(:,:),yIntScale(:),yTempV(:,:),yF(:,:),yT(:,:)
      double precision, allocatable :: cos2x0(:,:),cos2xp(:,:),cos2xm(:,:),cos2y(:,:),cosx(:,:),sinx(:,:)

      double precision mu12,r0diatom,dDiatom
      common/MassInfo/mu12,r0diatom,dDiatom

      allocate(xIntScale(xNumPoints),xTempV(LegPoints),xS(xDim,xDim),xT(xDim,xDim))
      allocate(siny(LegPoints,yNumPoints),cosy(LegPoints,yNumPoints),tany(LegPoints,yNumPoints),yIntScale(yNumPoints))
      allocate(cos2x0(LegPoints,xNumPoints),cos2xp(LegPoints,xNumPoints),cos2xm(LegPoints,xNumPoints),cos2y(LegPoints,yNumPoints))
      allocate(sinx(LegPoints,xNumPoints),cosx(LegPoints,xNumPoints))
      allocate(yTempV(LegPoints,yNumPoints),yF(yDim,yDim),yT(yDim,yDim))
      allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
      allocate(kyMin(yDim,yDim),kyMax(yDim,yDim))
      allocate(Pot(LegPoints,LegPoints,yNumPoints,xNumPoints))

      Pi = 3.1415926535897932385d0

      H = 0.0d0

      m = -0.5d0/(mu*R*R)

c     m = -0.5d0/(mu12*R*R)

      do kx = 1,xNumPoints-1
         ax = xPoints(kx)
         bx = xPoints(kx+1)
         xIntScale(kx) = 0.5d0*(bx-ax)
         xScaledZero = 0.5d0*(bx+ax)
         do lx = 1,LegPoints
            x = xIntScale(kx)*xLeg(lx)+xScaledZero
c     cos2x0(lx,kx) = dcos(2.0d0*x)
            cosx(lx,kx) = dcos(x)
            sinx(lx,kx) = dsin(x)
c     cos2xp(lx,kx) = dcos(2.0d0*(x+Pi/3.0d0))
c     cos2xm(lx,kx) = dcos(2.0d0*(x-Pi/3.0d0))
         enddo
      enddo
      do ky = 1,yNumPoints-1
         ay = yPoints(ky)
         by = yPoints(ky+1)
         yIntScale(ky) = 0.5d0*(by-ay)
         yScaledZero = 0.5d0*(by+ay)
         do ly = 1,LegPoints
            y = yIntScale(ky)*xLeg(ly)+yScaledZero
c     cosy(ly,ky) = 4.0d0*dcos(4.0d0*y)
            cosy(ly,ky) = dcos(y)
            cos2y(ly,ky) = dcos(2.0d0*y)
c     siny(ly,ky) = dsin(4.0d0*y)
            siny(ly,ky) = dsin(y)
            tany(ly,ky) = 2.0d0*dtan(2.0d0*y)
         enddo
      enddo

      do ix = 1,xDim
         do ixp = 1,xDim
            kxMin(ixp,ix) = max(xBounds(ix),xBounds(ixp))
            kxMax(ixp,ix) = min(xBounds(ix+Order+1),xBounds(ixp+Order+1))-1
         enddo
      enddo
      do iy = 1,yDim
         do iyp = 1,yDim
            kyMin(iyp,iy) = max(yBounds(iy),yBounds(iyp))
            kyMax(iyp,iy) = min(yBounds(iy+Order+1),yBounds(iyp+Order+1))-1
         enddo
      enddo

      do ix = 1,xDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            xS(ixp,ix) = 0.0d0
            xT(ixp,ix) = 0.0d0
            do kx = kxMin(ixp,ix),kxMax(ixp,ix)
               xTempS = 0.0d0
               xTempT = 0.0d0
               do lx = 1,LegPoints
                  a = wLeg(lx)*xIntScale(kx)*u(lx,kx,ix)
                  b = a*u(lx,kx,ixp)
                  xTempS = xTempS + b
                  xTempT = xTempT + a*uxx(lx,kx,ixp)
               enddo
               xS(ixp,ix) = xS(ixp,ix) + xTempS
               xT(ixp,ix) = xT(ixp,ix) + xTempT
            enddo
         enddo
      enddo

      do iy = 1,yDim
         do iyp = max(1,iy-Order),min(yDim,iy+Order)
            yF(iyp,iy) = 0.0d0
            yT(iyp,iy) = 0.0d0
            do ky = kyMin(iyp,iy),kyMax(iyp,iy)
               yTempF = 0.0d0
               yTempT = 0.0d0
               do ly = 1,LegPoints
                  a = wLeg(ly)*yIntScale(ky)*v(ly,ky,iy)
                  ap = wLeg(ly)*yIntScale(ky)*vy(ly,ky,iy)
                  b = a*v(ly,ky,iyp)
                  bp = ap*vy(ly,ky,iyp)
                  yTempF = yTempF + b/siny(ly,ky)
                  yTempT = yTempT + a*(siny(ly,ky)*vyy(ly,ky,iyp)+cosy(ly,ky)*vy(ly,ky,iyp))
               enddo
               yF(iyp,iy) = yF(iyp,iy) + yTempF
               yT(iyp,iy) = yT(iyp,iy) + yTempT
            enddo
         enddo
      enddo

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            do iy = 1,yDim
               Row = i1+iy
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  Col = i1p+iyp
                  if (Col .ge. Row) then
                     NewRow = HalfBandWidth+1+Row-Col
                     H(NewRow,Col) = m*(xT(ixp,ix)*yF(iyp,iy) + xS(ixp,ix)*yT(iyp,iy))
                  endif
               enddo
            enddo
         enddo
      enddo

c     if potential integral is not separable, use the following code section
c     to do 2D integrals

      do kx = 1,xNumPoints-1
         do ky = 1,yNumPoints-1
            do lx = 1,LegPoints
               do ly = 1,LegPoints
                  r12 = (2.d0*R*cosx(lx,kx)*siny(ly,ky))/dsqrt(2.d0)
                  r34 = (2.d0*R*sinx(lx,kx)*siny(ly,ky))/dsqrt(2.d0)
                  r13 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(cosx(lx,kx)-sinx(lx,kx)))/dsqrt(2.d0)
                  r14 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)
                  r23 = (dsqrt(2.d0)*R*cosy(ly,ky) - R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)
                  r24 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(-cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)
                  call  sumpairwisepot(r12, r13, r14, r23, r24, r34, potvalue)
                  Pot(ly,lx,ky,kx) = alpha*potvalue
               enddo
            enddo
         enddo
      enddo

      do ix = 1,xDim
         i1 = (ix-1)*yDim
         do ixp = max(1,ix-Order),min(xDim,ix+Order)
            i1p = (ixp-1)*yDim
            do iy = 1,yDim
               Row = i1+iy
               do iyp = max(1,iy-Order),min(yDim,iy+Order)
                  Col = i1p+iyp
                  if (Col .ge. Row) then
                     
                     NewRow = HalfBandWidth+1+Row-Col

                     VInt = 0.0d0
                     do ky = kyMin(iyp,iy),kyMax(iyp,iy)
                        do ly = 1,LegPoints
                           yTempV(ly,ky) = wLeg(ly)*siny(ly,ky)*v(ly,ky,iy)*v(ly,ky,iyp)
                        enddo
                     enddo
                     do kx = kxMin(ixp,ix),kxMax(ixp,ix)
                        do lx = 1,LegPoints
                           xTempV(lx) = wLeg(lx)*u(lx,kx,ix)*u(lx,kx,ixp)
                        enddo
                        do ky = kyMin(iyp,iy),kyMax(iyp,iy)
                           VTempInt = 0.0d0
                           do lx = 1,LegPoints
                              do ly = 1,LegPoints
                                 VTempInt = VTempInt + xTempV(lx)*yTempV(ly,ky)*Pot(ly,lx,ky,kx)
                              enddo
                           enddo
                           VInt = VInt + xIntScale(kx)*yIntScale(ky)*VTempInt
                        enddo
                     enddo

                     H(NewRow,Col) = H(NewRow,Col)+VInt
c                     write(25,*) ix,ixp,H(NewRow,Col)
                  endif
               enddo
            enddo
         enddo
      enddo

      deallocate(Pot)
      deallocate(xIntScale,xTempV,xS,xT)
      deallocate(cosy,siny,tany,yIntScale,yTempV,yF,yT)
      deallocate(cos2y,cos2x0,cos2xp,cos2xm)
      deallocate(kxMin,kxMax)
      deallocate(kyMin,kyMax)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P)

      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision P(NumStates,NumStates)

      integer i,j,k
      double precision a,ddot
      double precision, allocatable :: TempPsi1(:),TempPsi2(:)

      allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))

      a = 0.5d0/RDelt

      do j = 1,NumStates
         do k = 1,MatrixDim
            TempPsi1(k) = rPsi(k,j)-lPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
         do i = 1,NumStates
            P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
         enddo
      enddo

      deallocate(TempPsi1,TempPsi2)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CalcQMatrix(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,Q)
      
      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision Q(NumStates,NumStates)
      
      integer i,j,k
      double precision a,ddot
      double precision, allocatable :: TempPsi1(:),TempPsi2(:)
      
      allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))
      
      a = 1.0d0/(RDelt**2)
      
      do j = 1,NumStates
         do k = 1,MatrixDim
            TempPsi1(k) = lPsi(k,j)+rPsi(k,j)-2.0d0*mPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
         do i = 1,NumStates
            Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
         enddo
      enddo
      
      deallocate(TempPsi1,TempPsi2)
      
      return
      end
c      cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)
      
      integer NumStates,HalfBandWidth,MatrixDim,ncv
      double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
      double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)

      integer i,j
      double precision Phase,ddot
      double precision, allocatable :: TempPsi(:)
!      write(6,*) 'in FixPhase: allocating memory for TempPsi'
      allocate(TempPsi(MatrixDim))

      do i = 1,NumStates
!         write(6,*) 'in FixPhase: i = ', i
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
         Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
         if (Phase .lt. 0.0d0) then
            do j = 1,MatrixDim
               rPsi(j,i) = -rPsi(j,i)
            enddo
         endif
      enddo

      deallocate(TempPsi)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GridMaker222(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)

      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      integer i,j,k
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2
      double precision yDelt,y0,y1,y2

      Pi = 3.1415926535897932385d0

      r0New = 3.0d0*r0

      xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0))))

      if (R .gt. xRswitch) then
         x0 = xMin
         x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
         x2 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
            k = k + 1
         enddo
      else
         x0 = xMin
         x1 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
      endif

      yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))

      if (R .gt. yRswitch) then
         y0 = yMin
         y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
         y2 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints/2)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
         yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y1
            k = k + 1
         enddo
      else
         y0 = yMin
         y1 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints-1)
         do i = 1,yNumPoints
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
      endif

      return
      end


      
      subroutine GridMaker111(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      
      integer i,j,k
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2
      double precision yDelt,y0,y1,y2
      
      Pi = 3.1415926535897932385d0
      

      r0New = 3.0d0*r0
      
      xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0)) ))
      
      
      if (R .gt. xRswitch) then
         x0 = xMin
         x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
         x2 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
            k = k + 1
         enddo
      else
         x0 = xMin
         x1 = xMax
         k = 1
         xDelt = (x1-x0)/dfloat(xNumPoints-1)
         do i = 1,xNumPoints
            xPoints(k) = (i-1)*xDelt + x0
            k = k + 1
         enddo
      endif
      
      yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))
      
      if (R .gt. yRswitch) then
         y0 = yMin
         y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
         y2 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints/2)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
         yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
         do i = 1,yNumPoints/2
            yPoints(k) = (i-1)*yDelt + y1
            k = k + 1
         enddo
      else
         y0 = yMin
         y1 = yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints-1)
         do i = 1,yNumPoints
            yPoints(k) = (i-1)*yDelt + y0
            k = k + 1
         enddo
      endif
      
      return
      end
      

      subroutine GridMaker(m,mu,R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)

      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      integer i,j,k
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2
      double precision yDelt,y0,y1,y2

      Pi = 3.1415926535897932385d0

c     r0New = 3.0d0*r0

c     xRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0+dcos(2.0d0*Pi*(1/12.0d0+1/3.0d0))))

c     write(96,*) 'xRswitch,R=',xRswitch,R
c     if (R .gt. xRswitch) then
c     x0 = xMin
c     x1 = 0.5d0*dacos(dsqrt(3.0d0)*r0New**2/R**2-1.0d0) - Pi/3.0d0
c     x2 = xMax
c     write(96,*) 'x0,x1,x2=',x0,x1,x2
c     k = 1
c     xDelt = (x1-x0)/dfloat(xNumPoints/2)
c     do i = 1,xNumPoints
c     xPoints(k) = (i-1)*xDelt + x0
c     k = k + 1
c     enddo
c     c       xDelt = (x2-x1)/dfloat(xNumPoints/2-1)
c     c       do i = 1,xNumPoints/2
c     c        xPoints(k) = (i-1)*xDelt + x1
c     c        k = k + 1
c     c       enddo
c     c 
c     xDelt = dsqrt(x2-x1)/dfloat(xNumPoints/2-1)
c     do i = 1,xNumPoints/2
c     xPoints(k) = x2-(xNumPoints/2-i)**2*xDelt*xDelt
c     k = k + 1
c     enddo

c     else
      x0 = xMin
      x1 = xMax
c     write(96,*) 'x0,x1=',x0,x1
      k = 1
      xDelt = (x1-x0)/dfloat(xNumPoints-1)
      do i = 1,xNumPoints
         xPoints(k) = (i-1)*xDelt + x0
c     xPoints(i) = (i-1)*xDelt + x0
         k = k + 1
      enddo
c     endif
c     write(96,15) (xPoints(k),k=1,xNumPoints)
 15   format(6(1x,1pd12.5))
      

c     yRswitch = dsqrt(dsqrt(3.0d0)*r0New**2/(1.0d0-dcos(Pi/4.0d0)))

c     if (R .gt. yRswitch) then
c     y0 = yMin
c     y1 = 0.5d0*dacos(1.0d0-dsqrt(3.0d0)*r0New**2/R**2)
c     y2 = yMax
c     write(96,*) 'y0,y1,y2=',y0,y1,y2
c     k = 1
c     yDelt = dsqrt(y1-y0)/dfloat(yNumPoints/2)
c     do i = 1,yNumPoints/2
c     yPoints(k) = ((i-1)*yDelt)**2 + y0
c     k = k + 1
c     enddo
c     yDelt = (y2-y1)/dfloat(yNumPoints/2-1)
c     do i = 1,yNumPoints/2
c     yPoints(k) = (i-1)*yDelt + y1
c     k = k + 1
c     enddo
c     else
      y0 = yMin
      y1 = yMax
c     write(96,*) 'y0,y1=',y0,y1
      k = 1
      yDelt = (y1-y0)/dfloat(yNumPoints-1)
      do i = 1,yNumPoints
         yPoints(k) = (i-1)*yDelt + y0
         k = k + 1
      enddo
c     endif
c     write(96,15) (yPoints(k),k=1,yNumPoints)

      return
      end

      subroutine GridMakerBetter(m,mu,R,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)

      integer xNumPoints,yNumPoints
      double precision m,mu,R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

      integer i,j,OPGRID
      double precision Pi
      double precision r0New
      double precision xRswitch,yRswitch
      double precision xDelt,x0,x1,x2,x3
      double precision yDelt,y0,y1,y2,y3

      Pi = 3.1415926535897932385d0

      x0 = xMin
      x1 = xMax
c     write(96,*) 'x0,x1=',x0,x1
      k = 1
      xDelt = (x1-x0)/dfloat(xNumPoints-1)
      do i = 1,xNumPoints
         xPoints(k) = (i-1)*xDelt + x0
c     xPoints(i) = (i-1)*xDelt + x0
         k = k + 1
      enddo
      y0 = yMin
      y1 = yMax
c     write(96,*) 'y0,y1=',y0,y1
      k = 1
      yDelt = (y1-y0)/dfloat(yNumPoints-1)
      do i = 1,yNumPoints
         yPoints(k) = (i-1)*yDelt + y0
         k = k + 1
      enddo
      
c     write(96,15) (yPoints(k),k=1,yNumPoints)

c
c     write(96,15) (xPoints(k),k=1,xNumPoints)
      OPGRID=1
      if(OPGRID.eq.1) then
c     print*, 'R>xRswitch!! using modified grid!!'
         x0 = xMin
         x1=x0+Pi/(2.0d0*8.0d0)
         x2 = xMax-Pi/(2.0d0*8.0d0)
         x3=xMax
         k = 1
c         write(6,*) ' The phi grid:'
         xDelt = (x1-x0)/dfloat(xNumPoints/4)
         do i = 1,xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x0
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
         xDelt = (x2-x1)/dfloat(xNumPoints/2)
         do i = 1,xNumPoints/2
            xPoints(k) = (i-1)*xDelt + x1
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
         xDelt = (x3-x2)/dfloat(xNumPoints/4-1)
         do i = 1,xNumPoints/4
            xPoints(k) = (i-1)*xDelt + x2
c            print*, k, xPoints(k), xDelt
            k = k + 1
         enddo
         
         do k=1,2
            do i=2,xNumPoints-1
               xPoints(i)=(xPoints(i-1)+2.0d0*xPoints(i)+xPoints(i+1))/4.0d0
            enddo
         enddo

c         write(6,*) 'the theta grid:'
         y0 = yMin
         y1=y0+Pi/(2.0d0*3.0d0)
         y2=yMax
         k = 1
         yDelt = (y1-y0)/dfloat(yNumPoints/4)
         do i = 1,yNumPoints/4
            yPoints(k) = (i-1)*yDelt + y0
c            print*, k, yPoints(k), yDelt
            k = k + 1
         enddo
         yDelt = (y2-y1)/dfloat(3*yNumPoints/4-1)
         do i = 1,3*yNumPoints/4
            yPoints(k) = (i-1)*yDelt + y1
c            print*, k, yPoints(k), yDelt
            k = k + 1
         enddo
         do k=1,2
            do i=2,yNumPoints-1
               yPoints(i)=(yPoints(i-1)+2.0d0*yPoints(i)+yPoints(i+1))/4.0d0
            enddo
         enddo

      endif
      
      
      
 15   format(6(1x,1pd12.5))
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CalcCoupling(NumStates,HalfBandWidth,MatrixDim,RDelt,lPsi,mPsi,rPsi,S,P,Q,dP)

      integer NumStates,HalfBandWidth,MatrixDim
      double precision RDelt
      double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
      double precision S(HalfBandWidth+1,MatrixDim)
      double precision P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates)

      integer i,j,k
      double precision aP,aQ,ddot
      double precision, allocatable :: lDiffPsi(:),rDiffPsi(:),TempPsi(:),TempPsiB(:),rSumPsi(:)

      allocate(lDiffPsi(MatrixDim),rDiffPsi(MatrixDim),TempPsi(MatrixDim),TempPsiB(MatrixDim),rSumPsi(MatrixDim))

      aP = 0.5d0/RDelt
      aQ = aP*aP

      do j = 1,NumStates
         do k = 1,MatrixDim
            rDiffPsi(k) = rPsi(k,j)-lPsi(k,j)
            rSumPsi(k)  = lPsi(k,j)+mPsi(k,j)+rPsi(k,j)
         enddo
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
         call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)
         do i = 1,NumStates
            P(i,j) = aP*ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
            dP(i,j)= ddot(MatrixDim,mPsi(1,i),1,TempPsiB,1)
            do k = 1,MatrixDim
               lDiffPsi(k) = rPsi(k,i)-lPsi(k,i)
            enddo
            Q(i,j) = -aQ*ddot(MatrixDim,lDiffPsi,1,TempPsi,1)
         enddo
      enddo

      do j=1,NumStates
	 do i=j,NumStates
            dP(i,j)=2.d0*aQ*(dP(i,j)-dP(j,i))
            dP(j,i)=-dP(i,j)
	 enddo
      enddo

      deallocate(lDiffPsi,rDiffPsi,TempPsi,rSumPsi,TempPsiB)

      return
      end
c      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Double Precision Function Vpot95(r)
      implicit real*8(a-h,o-z)
      double precision eps,D,c6,c8,c10,Astar,alpha,beta
      double precision rmin,Fscale,vscale,vlr

c     convert from angstroms to a.u.
      rmin = 2.96830d0/0.529177249d0

c     convert from Kelvin to a.u.
      eps = 10.956d0*3.166829d-6

      astar = 1.86924404d5
      alpha = 10.5717543
      beta = -2.07758779d0
      c6 = 1.35186623d0
      c8 = 0.41495143d0
      c10 = 0.17151143d0
      D = 1.438d0

      xscale = r/rmin

      vscale = Astar*dexp(-alpha*xscale + beta*xscale*xscale)

      vlr = c6/(xscale**6) + c8/(xscale**8) + c10/(xscale**10)

      if (xscale .gt. D)then
         Fscale = 1.0d0
      else
         Fscale = dexp(-(D/xscale-1.0d0)*(D/xscale-1.0d0))
      endif 

      Vpot95 = eps*(vscale - Fscale*vlr)

      return
      end


c     SUBROUTINE potential3(dim, DIMEN, DIMMAX, wlkdmc, vx, vy, vz, 
c     *                      r_ij_old, mass, e_pot)
      DOUBLE PRECISION FUNCTION Vpot91(r)
      IMPLICIT REAL*8(A-H,O-Z)

CCCCC LM2M2 with add-on potential: R. A. Aziz, and M. J. Slaman
CCCCC J. Chem. Phys. 94, 8047 (1991)

      INTEGER dim, DIMEN, DIMMAX, wlkdmc
      DOUBLE PRECISION mass, e_pot
      INTEGER i, j, k, position

      DOUBLE PRECISION B, alpha, beta, A, epsilo, r_m, r
      DOUBLE PRECISION c_6, c_8, c_10, d, x, fhelp, sum, x_0, x_1

      DOUBLE PRECISION PI, add_on

      PARAMETER(B = 0.0026d0, A = 1.89635353d0 * 10**5)
      PARAMETER(alpha = 10.70203539d0, beta = -1.90740649d0)
      PARAMETER(c_6 = 1.34687065d0, c_8 = 0.41308398d0, 
     *     c_10 = 0.17060159d0)
      PARAMETER(epsilo = 3.4739577d0 / 10**5, D = 1.4088d0, 
     *     r_m = 5.6115d0)
      PARAMETER(x_0 = 1.003535949d0, x_1 = 1.454790369d0)
      PARAMETER(PI = 3.14159265358979324d0)
C     PI = 2.0d0 * ASIN(1.0d0)

      e_pot = 0.0d0      
C     C      DO j = 1, wlkdmc
C     C         DO i = 1, dim
C     C            DO k = i + 1, dim
C     C               position = (i - 1) * dim - (i - 1) * i / 2 + k - i
      x = r / r_m
      IF(x .LT. D) THEN
         fhelp = exp(-(D / x - 1.0d0) * (D / x - 1.0d0))
      ELSE
         fhelp = 1.0d0
      ENDIF
      sum = c_6 / x**6 + c_8 / x**8 + c_10 / x**(10)
      IF ((x .LE. x_1) .AND. (x .GE. x_0)) THEN
         add_on = B * 
     *        (SIN(2.0d0 * PI * (x - x_0) / (x_1 - x_0) - 
     *        0.5d0 * PI) + 1.0d0) 
      ELSE
         add_on = 0.0d0
      ENDIF
      e_pot = e_pot + 
     *     (A * exp(-alpha * x + beta * x * x) - 
     *     fhelp * sum + add_on) * epsilo 
C     C            END DO
C     C         END DO
C     C      END DO
      Vpot91=e_pot
      RETURN
      END

      DOUBLE PRECISION FUNCTION Vpot(r)
      IMPLICIT REAL*8(A-H,O-Z)
      double precision kelvinPERau, angstromPERau, eps, sigma
      data kelvinPERau,angstromPERau/315933.16d0,0.529177249d0/
      data eps,sigma/35.6d0,2.749d0/
      ratio = (sigma/(angstromPERau*r))**6
      Vpot = 4.d0 * (eps/kelvinPERau) * (ratio**2-ratio)
      return
      end

Subroutine adiabaticSolver(NumStates,PsiFlag,CouplingFlag,LegendreFile,LegPoints,Shift,Order,Left,Right,Top,Bottom,alpha,massarray,&
     xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,RSteps,RDerivDelt,RFirst,RLast,V2Depth,&
     R,Uad,Psi,eDim,psiDim,S,sDim)
  implicit none

  integer LegPoints,xNumPoints,yNumPoints
  integer NumStates,PsiFlag,Order,Left,Right,Bottom,Top
  integer RSteps,CouplingFlag,CalcNewBasisFunc
  double precision alpha,m,Shift,V2Depth
  double precision RLeft,RRight,RDerivDelt
  DOUBLE PRECISION RFirst,RLast,XFirst,XLast,StepX
  double precision xMin,xMax,yMin,yMax
  double precision :: R(RSteps)
  !double precision, allocatable :: R(:)
  double precision, allocatable :: xPoints(:),yPoints(:)

  logical, allocatable :: Select(:)

  integer iparam(11),ncv,info
  integer i,j,k,iR
  integer LeadDim,MatrixDim,HalfBandWidth
  integer xDim,yDim,eDim,psiDim,sDim
  integer, allocatable :: iwork(:)
  integer, allocatable :: xBounds(:),yBounds(:)
  double precision Tol
  double precision TotalMemory
  double precision xNUM,dNUM
  double precision mu, mu12,mu34,mu1234,massarray(4),r0diatom, dDiatom, etaOVERpi, Pi
  double precision u1,v1,sys_ss_pot
  double precision, allocatable :: LUFac(:,:),workl(:)
  double precision, allocatable :: workd(:),Residuals(:)
  double precision, allocatable :: xLeg(:),wLeg(:)
  double precision, allocatable :: u(:,:,:),uxx(:,:,:),v(:,:,:)
  double precision, allocatable :: vy(:,:,:),vyy(:,:,:)
  double precision, allocatable :: H(:,:),Vmid(:,:,:),Vexp(:,:,:)
!  double precision, allocatable :: mPsi(:,:),lPsi(:,:),rPsi(:,:)
  double precision, allocatable :: mPsi(:,:)
  double precision :: Psi(RSteps,psiDim,eDim),Uad(RSteps,eDim,2),S(sDim,psiDim)
  double precision, allocatable ::P(:,:),Q(:,:),dP(:,:),Energies(:,:)
  double precision r0
  common/MassInfo/r0diatom,dDiatom

  character*64 LegendreFile
  character*64 filename
  filename="AdiabaticCurves.dat"
  
  open(unit=100,file=filename)
  write(*,*) NumStates,PsiFlag,CouplingFlag

  !     read in Gauss-Legendre info
  write(*,1002) LegendreFile
  write(*,*) LegPoints,' LegPoints'
  write(*,*) Shift,Order,Left,Right,Bottom,Top

  ! end if
  write(*,*) alpha,massarray(1),massarray(2),massarray(3),massarray(4)
  Pi=dacos(-1.d0)
  write(*,*) 'Pi=',Pi
  write(*,*) xNumPoints,xMin,xMax
  write(*,*) yNumPoints,yMin,yMax
  write(*,*) RSteps,RDerivDelt,RFirst,RLast
  write(*,*) V2Depth


  mu = (massarray(1)*massarray(2)*massarray(3)*massarray(4))
  !mu = mu/(massarray(1)+massarray(2)+massarray(3)+massarray(4))
  mu = mu**0.25d0
!  mu = (m1*m2*m3*m4/(m1 + m2 + m3 + m4))**(1.0d0/3.0d0)
  mu12 = massarray(1)*massarray(2)/(massarray(1) + massarray(2))
  mu34 = massarray(3)*massarray(4)/(massarray(3) + massarray(4))
  mu1234 = (massarray(1) + massarray(2))*(massarray(3) + massarray(4))/(massarray(1) + massarray(2) + massarray(3) + massarray(4))

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

  TotalMemory = (4*(HalfBandWidth+1)+(3*HalfBandWidth+1)+6*NumStates)*8.0d0*MatrixDim
  TotalMemory = TotalMemory/(1024.0d0*1024.0d0)

  write(6,*)
  write(6,*) 'MatrixDim ',MatrixDim
  write(6,*) 'HalfBandWidth ',HalfBandWidth
  write(6,*) 'Approximate peak memory usage (in Mb) ',TotalMemory
  write(6,*)

  allocate(xPoints(xNumPoints),yPoints(yNumPoints))
  allocate(xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order))
  allocate(u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim))
  allocate(v(LegPoints,yNumPoints,yDim),vy(LegPoints,yNumPoints,yDim),vyy(LegPoints,yNumPoints,yDim))
  allocate(H(HalfBandWidth+1,MatrixDim))

  allocate(Vmid(HalfBandWidth+1,MatrixDim,6))
  allocate(Vexp(NumStates,NumStates,6))

  allocate(P(NumStates,NumStates),Q(NumStates,NumStates),dP(NumStates,NumStates))

  ncv = NumStates*2d0
  LeadDim = 3*HalfBandWidth+1
  allocate(iwork(MatrixDim))
  allocate(Select(ncv))
  allocate(LUFac(LeadDim,MatrixDim))
  allocate(workl(ncv*ncv+8*ncv))
  allocate(workd(3*MatrixDim))
!  allocate(lPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv),mPsi(MatrixDim,ncv))
  allocate(mPsi(MatrixDim,ncv))
  allocate(Residuals(MatrixDim))
  allocate(Energies(ncv,2))
  info=0

  call GridMaker(R(1),11.65d0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)

  print*, "after gridmaker..."
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,0,u)
  call CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,xDim,xBounds,xNumPoints,2,uxx)
  call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,yDim,yBounds,yNumPoints,0,v)
  call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,yDim,yBounds,yNumPoints,1,vy)
  call CalcBasisFuncs(Bottom,Top,Order,yPoints,LegPoints,xLeg,yDim,yBounds,yNumPoints,2,vyy)

  call CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,&
       xNumPoints,yNumPoints,u,v,xBounds,yBounds,HalfBandWidth,S)

  do iR = 1,RSteps
     
     call CalcHamiltonian(alpha,R(iR),mu,mu12,mu34,mu1234,&
          massarray(1),massarray(2),massarray(3),massarray(4),Order,xPoints,yPoints,LegPoints,xLeg,&
          wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,vy,uxx,vyy,xBounds,yBounds,HalfBandWidth,V2Depth,H)
     call MyDsband(Select,Energies,mPsi,MatrixDim,Shift,MatrixDim,H,S,HalfBandWidth+1,LUFac,LeadDim,&
          HalfBandWidth,NumStates,Tol,Residuals,ncv,mPsi,MatrixDim,iparam,workd,workl,ncv*ncv+8*ncv,iwork,info)
     if (info.ne.0) write(6,*) 'Error in MyDsband.  info = ',info
     
     
     if (iR .gt. 1) then
        call FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,Psi(iR-1,:,:),mPsi)
     endif
     
     if(iparam(5).lt.numStates) write(*,*) "LESS THAN NUMSTATES"

     call CalcEigenErrors(info,iparam,MatrixDim,H,HalfBandWidth+1,S,HalfBandWidth,NumStates,mPsi,Energies,ncv)
     IF(R(iR).GT. 2.2d0) Shift = 1.1d0*Energies(1,1)
     
     write(6,*)
     write(6,*) R(iR)
     do i = 1,NumStates
        !        write(6,*) i,(2d0*mu*R(iR)**2)*Energies(i,1),Energies(i,2)
        write(6,*) i,Energies(i,1),Energies(i,2)
     enddo

     !     write(100,10) R(iR),(2d0*mu*R(iR)*R(iR)*Energies(j,1),j=1,numStates)
     write(100,10) R(iR),(Energies(j,1),j=1,numStates)

     do i=1,eDim
        Uad(iR,i,1)=Energies(i,1)
        Uad(iR,i,2)=Energies(i,2)
     end do
     do i=1,MatrixDim
        do j=1,eDim
           Psi(iR,i,j)=mPsi(i,j)
        end do
     end do


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
  close(100)
  !close(200)
  !close(101)
  !close(102)
  !close(103)
  deallocate(H)
  deallocate(iwork)
  deallocate(Select)
  deallocate(LUFac)
  deallocate(workl)
  deallocate(workd)
!  deallocate(lPsi,rPsi)
  deallocate(Residuals)
  deallocate(P,Q,dP)
  deallocate(xPoints,yPoints)
  deallocate(xLeg,wLeg)
  deallocate(xBounds,yBounds)
  deallocate(u,uxx)
  deallocate(v,vy,vyy)
  deallocate(Vmid,Vexp)

10 format(1P,1000e25.15)
20 format(1P,2000e16.8)
1002 format(a64)


end Subroutine adiabaticSolver

subroutine CalcOverlap(Order,xPoints,yPoints,LegPoints,xLeg,wLeg,xDim,yDim,xNumPoints,yNumPoints,&
     u,v,xBounds,yBounds,HalfBandWidth,S)

  integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints,xBounds(xNumPoints+2*Order),yBounds(yNumPoints+2*Order),HalfBandWidth
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
        !     siny(ly,ky) = dsin(4.0d0*y)
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
end subroutine CalcOverlap

subroutine CalcHamiltonian(alpha,R,mu,mu12,mu34,mu1234,m1,m2,m3,m4,Order,xPoints,yPoints,LegPoints,&
     xLeg,wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,vy,uxx,vyy,xBounds,yBounds,HalfBandWidth,V2Depth,H)
  implicit none
  integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints,xBounds(*),yBounds(*),HalfBandWidth

  double precision V2Depth
  character(len=74) :: PotFile
  double precision alpha,R,mu,mu12,mu34,mu1234,m1,m2,m3,m4,kxy,kxz,kyz
  double precision xPoints(*),yPoints(*),xLeg(*),wLeg(*)
  double precision H(HalfBandWidth+1,xDim*yDim)
  double precision u(LegPoints,xNumPoints,xDim),uxx(LegPoints,xNumPoints,xDim)
  double precision v(LegPoints,yNumPoints,yDim),vy(LegPoints,yNumPoints,yDim),vyy(LegPoints,yNumPoints,yDim)

  integer ix,iy,ixp,iyp,kx,ky,lx,ly
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)
  double precision a,b,mr,Pi, ap, bp
  double precision Rall,r12,r12a,r23,r23a,r23b,r23c,r13,r13a,r13b,r13c,r14,r24,r34,rp1,rp2,rp3
  double precision u1,sys_ss_pot,V12,V23,V31
  double precision VInt,VTempInt,potvalue
  !     double precision TempPot,VInt,VTempInt
  double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt, y1, y2, y3, x1, x2, x3, x4
  double precision y,ay,by,yScaledZero,yTempT,yTempF,yInt
  double precision, allocatable :: Pot(:,:,:,:)
  double precision, allocatable :: xIntScale(:),xTempV(:),xS(:,:),xT(:,:)
  double precision, allocatable :: siny(:,:),cosy(:,:),tany(:,:),yIntScale(:),yTempV(:,:),yF(:,:),yT(:,:)
  double precision, allocatable :: cos2x0(:,:),cos2xp(:,:),cos2xm(:,:),cos2y(:,:),cosx(:,:),sinx(:,:)
  double precision, allocatable :: rvec(:),vvec(:)


  double precision r0diatom,dDiatom
  common/MassInfo/r0diatom,dDiatom

  kxy=0.0d0
  kxz=0.0d0
  kyz=0.0d0
  
  allocate(xIntScale(xNumPoints),xTempV(LegPoints),xS(xDim,xDim),xT(xDim,xDim))
  allocate(siny(LegPoints,yNumPoints),cosy(LegPoints,yNumPoints),tany(LegPoints,yNumPoints),yIntScale(yNumPoints))
  allocate(cos2x0(LegPoints,xNumPoints),cos2xp(LegPoints,xNumPoints),cos2xm(LegPoints,xNumPoints),cos2y(LegPoints,yNumPoints))
  allocate(sinx(LegPoints,xNumPoints),cosx(LegPoints,xNumPoints))
  allocate(yTempV(LegPoints,yNumPoints),yF(yDim,yDim),yT(yDim,yDim))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(kyMin(yDim,yDim),kyMax(yDim,yDim))
  allocate(Pot(LegPoints,LegPoints,yNumPoints,xNumPoints))
  allocate(rvec(161),vvec(161))

  Pi = 3.1415926535897932385d0

  H = 0.0d0
  
  mr = -0.5d0/(mu*R*R)

  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x = xIntScale(kx)*xLeg(lx)+xScaledZero
        !     cos2x0(lx,kx) = dcos(2.0d0*x)
        cosx(lx,kx) = dcos(x)
        sinx(lx,kx) = dsin(x)
        !     cos2xp(lx,kx) = dcos(2.0d0*(x+Pi/3.0d0))
        !     cos2xm(lx,kx) = dcos(2.0d0*(x-Pi/3.0d0))
     enddo
  enddo
  do ky = 1,yNumPoints-1
     ay = yPoints(ky)
     by = yPoints(ky+1)
     yIntScale(ky) = 0.5d0*(by-ay)
     yScaledZero = 0.5d0*(by+ay)
     do ly = 1,LegPoints
        y = yIntScale(ky)*xLeg(ly)+yScaledZero
        !     cosy(ly,ky) = 4.0d0*dcos(4.0d0*y)
        cosy(ly,ky) = dcos(y)
        cos2y(ly,ky) = dcos(2.0d0*y)
        !     siny(ly,ky) = dsin(4.0d0*y)
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
              !     yTempF = yTempF + b*tany(ly,ky)
              yTempF = yTempF + b/siny(ly,ky)
              yTempT = yTempT + a*(siny(ly,ky)*vyy(ly,ky,iyp)+cosy(ly,ky)*vy(ly,ky,iyp))
              !                  yTempT = yTempT + bp*siny(ly,ky)
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
                 H(NewRow,Col) = mr*(xT(ixp,ix)*yF(iyp,iy)+xS(ixp,ix)*yT(iyp,iy)) !! previously was minus
              endif
           enddo
        enddo
     enddo
  enddo

  !     if potential integral is not separable, use the following code section
  !     to do 2D integrals

  !      Rall = R/dsqrt(dsqrt(3.0d0))

  !      write(6,*) 'r0diatom=',r0diatom,' Rall=',Rall,' dDiatom=',dDiatom

  do kx = 1,xNumPoints-1
     do ky = 1,yNumPoints-1
        do lx = 1,LegPoints
           do ly = 1,LegPoints
              !     r12 = Rall*dsqrt(1.0d0+cos2y(ly,ky)*cos2x0(lx,kx))
              !cc                 r12 = (2.d0*R*cosx(lx,kx)*siny(ly,ky))/dsqrt(2.d0)
              !cc                  r34 = (2.d0*R*sinx(lx,kx)*siny(ly,ky))/dsqrt(2.d0)
              !cc                  r13 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(cosx(lx,kx)-sinx(lx,kx)))/dsqrt(2.d0)
              !cc                  r14 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)
              !cc                  r23 = (dsqrt(2.d0)*R*cosy(ly,ky) - R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)
              !cc                  r24 = (dsqrt(2.d0)*R*cosy(ly,ky) + R*siny(ly,ky)*(-cosx(lx,kx)+sinx(lx,kx)))/dsqrt(2.d0)

              !cc the correct version for equal masses                   
              !r12 = (2.d0**(1.d0/6.d0))*R*cosx(lx,kx)*siny(ly,ky)
              !r34 = (2.d0**(1.d0/6.d0))*R*sinx(lx,kx)*siny(ly,ky)
              !r13 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)-sinx(lx,kx))/2.d0
              !r14 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx))/2.d0
              !r23 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) - (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx))/2.d0
              !r24 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(-cosx(lx,kx)+sinx(lx,kx))/2.d0

              !cc the correct version for different masses
!!$              r12 = (mu/mu12)**(1.d0/2.d0)*R*cosx(lx,kx)*siny(ly,ky)
!!$              r34 = (mu/mu34)**(1.d0/2.d0)*R*sinx(lx,kx)*siny(ly,ky)
!!$              rp1 = (mu/mu1234)**(1.d0/2.d0)*R*cosy(ly,ky)
!!$              rp2 = (mu/mu12)**(1.d0/2.d0)/(m1+m2)*R*cosx(lx,kx)*siny(ly,ky)
!!$              rp3 = (mu/mu34)**(1.d0/2.d0)/(m3+m4)*R*sinx(lx,kx)*siny(ly,ky)
!!$              r13 = rp1+rp2*m2-rp3*m4
!!$              r14 = rp1+rp2*m2+rp3*m3
!!$              r23 = rp1-rp2*m1-rp3*m4
!!$              r24 = rp1-rp2*m1+rp3*m3

              !NPM Version for different masses in H-type mass-scaled Jacobi coordinates
              y1=R*cosx(lx,kx)*siny(ly,ky)
              y2=R*sinx(lx,kx)*siny(ly,ky)
              y3=R*cosy(ly,ky)
              x1=Sqrt(mu)*((m2*y1)/((m1 + m2)*Sqrt(mu12)) - ((m3 + m4)*y3)/((m1 + m2 + m3 + m4)*Sqrt(mu1234)))
              x2=Sqrt(mu)*(-((m1*y1)/((m1 + m2)*Sqrt(mu12))) - ((m3 + m4)*y3)/((m1 + m2 + m3 + m4)*Sqrt(mu1234)))
              x3=Sqrt(mu)*((m4*y2)/((m3 + m4)*Sqrt(mu34)) + ((m1 + m2)*y3)/((m1 + m2 + m3 + m4)*Sqrt(mu1234)))
              x4=Sqrt(mu)*(-((m3*y2)/((m3 + m4)*Sqrt(mu34))) + ((m1 + m2)*y3)/((m1 + m2 + m3 + m4)*Sqrt(mu1234)))
              r12=abs(x1-x2)
              r13=abs(x1-x3)
              r14=abs(x1-x4)
              r23=abs(x2-x3)
              r24=abs(x2-x4)
              r34=abs(x3-x4)
              !write(6,*) V2Depth
              !     c          u1 = sys_ss_pot(r12,v12,2,.FALSE.)
              !ccc  v12 = dexp(-r12**2/200.d0)
              !ccc  v12 = -dDiatom/(dcosh(r12/r0diatom))**2
              !ccc  v12 = Vpot(r12)
              !     r23 = Rall*dsqrt(1.0d0+cos2y(ly,ky)*cos2xp(lx,kx))
              !ccc  u1 = sys_ss_pot(r23,v23,2,.FALSE.)
              !ccc  v23 = dexp(-r23**2/200.d0)
              !ccc  v23 = -dDiatom/(dcosh(r23/r0diatom))**2
              !ccc  v23 = Vpot(r23)
              !     r13 = Rall*dsqrt(1.0d0+cos2y(ly,ky)*cos2xm(lx,kx))
              !ccc  u1 = sys_ss_pot(r13,v31,2,.FALSE.)
              !ccc  v31 = dexp(-r13**2/200.d0)
              !ccc  v31 = -dDiatom/(dcosh(r13/r0diatom))**2
              !ccc  v31 = Vpot(r13)
              !     call  h3ppot(r12, r13, r23, potvalue)


              call  sumpairwisepot(r12,r13,r14,r23,r24,r34,V2Depth,potvalue,rvec,vvec)
              Pot(ly,lx,ky,kx) = alpha*potvalue
              !Pot(ly,lx,ky,kx) = alpha*(0.5d0*mu*(y1**2 + y2**2 + y3**2) + kxy*y1*y2 + kxz*y1*y3 + kyz*y2*y3)

              !write(6,*) potvalue
              !     Pot(ly,lx,ky,kx) = 0.d0
              !cc   Pot(ly,lx,ky,kx) = alpha*(V12+V23+V31)
              !cc   Pot(ly,lx,ky,kx) = alpha*(V12*V23*V31)
              !cc   Pot(ly,lx,ky,kx) = alpha*(1.d0/r12**2+1.d0/r23**2+1.d0/r13**2)
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

                 H(NewRow,Col) = H(NewRow,Col) + VInt
                 !                     write(25,*) ix,ixp,H(NewRow,Col)
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
  deallocate(rvec,vvec)

  return
end subroutine CalcHamiltonian

subroutine CalcVMatrix(NumStates,HalfBandWidth,MatrixDim,mPsi,S,Vexp,R,Order,xPoints,yPoints,&
     LegPoints,xLeg,wLeg,xDim,yDim,xNumPoints,yNumPoints,u,v,xBounds,yBounds,Vmid)

  integer Order,LegPoints,xDim,yDim,xNumPoints,yNumPoints,xBounds(*),yBounds(*)


  integer NumStates,HalfBandWidth,MatrixDim
  !double precision lPsi(MatrixDim,NumStates),mPsi(MatrixDim,NumStates),rPsi(MatrixDim,NumStates)
  double precision mPsi(MatrixDim,NumStates)
  double precision S(HalfBandWidth+1,MatrixDim)
  double precision Vexp(NumStates,NumStates,6)

  integer i,j,k,idv
  double precision ddot
  double precision, allocatable :: TempPsi1(:),TempPsi2(:)



  double precision R,L,dd
  double precision xPoints(*),yPoints(*),xLeg(*),wLeg(*)

  double precision u(LegPoints,xNumPoints,xDim)
  double precision v(LegPoints,yNumPoints,yDim)


  integer ix,iy,ixp,iyp,kx,ky,lx,ly
  integer i1,i1p
  integer Row,NewRow,Col
  integer, allocatable :: kxMin(:,:),kxMax(:,:),kyMin(:,:),kyMax(:,:)

  double precision r12,r34,r13,r14,r24,r23
  double precision VInt,VTempInt
  !     double precision TempPot,VInt,VTempInt
  double precision x,ax,bx,xScaledZero,xTempT,xTempS,xInt
  double precision y,ay,by,yScaledZero,yTempT,yTempF,yInt
  double precision, allocatable :: Vr(:,:,:,:,:)
  double precision, allocatable :: xIntScale(:),xTempV(:)
  double precision, allocatable :: siny(:,:),cosy(:,:),yIntScale(:),yTempV(:,:)
  double precision, allocatable :: cosx(:,:),sinx(:,:)

  double precision Vmid(HalfBandWidth+1,xDim*yDim,6)

  allocate(xIntScale(xNumPoints),xTempV(LegPoints))
  allocate(siny(LegPoints,yNumPoints),cosy(LegPoints,yNumPoints),yIntScale(yNumPoints))
  allocate(sinx(LegPoints,xNumPoints),cosx(LegPoints,xNumPoints))
  allocate(yTempV(LegPoints,yNumPoints))
  allocate(kxMin(xDim,xDim),kxMax(xDim,xDim))
  allocate(kyMin(yDim,yDim),kyMax(yDim,yDim))
  allocate(Vr(LegPoints,LegPoints,yNumPoints,xNumPoints,6))


  allocate(TempPsi1(MatrixDim),TempPsi2(MatrixDim))


  dd=6.272844447382345d0 ! a2even=20, with one deep bound state L=1
  dd=dd*0.9d0
  L=1.0d0
  do kx = 1,xNumPoints-1
     ax = xPoints(kx)
     bx = xPoints(kx+1)
     xIntScale(kx) = 0.5d0*(bx-ax)
     xScaledZero = 0.5d0*(bx+ax)
     do lx = 1,LegPoints
        x = xIntScale(kx)*xLeg(lx)+xScaledZero
        cosx(lx,kx) = dcos(x)
        sinx(lx,kx) = dsin(x)
     enddo
  enddo
  do ky = 1,yNumPoints-1
     ay = yPoints(ky)
     by = yPoints(ky+1)
     yIntScale(ky) = 0.5d0*(by-ay)
     yScaledZero = 0.5d0*(by+ay)
     do ly = 1,LegPoints
        y = yIntScale(ky)*xLeg(ly)+yScaledZero
        cosy(ly,ky) = dcos(y)
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

  do kx = 1,xNumPoints-1
     do ky = 1,yNumPoints-1
        do lx = 1,LegPoints
           do ly = 1,LegPoints                             
              r12 = (2.d0**(1.d0/6.d0))*R*cosx(lx,kx)*siny(ly,ky)
              r34 = (2.d0**(1.d0/6.d0))*R*sinx(lx,kx)*siny(ly,ky)
              r13 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)-sinx(lx,kx))/2.d0
              r14 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx))/2.d0
              r23 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) - (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(cosx(lx,kx)+sinx(lx,kx))/2.d0                     
              r24 = (2.d0**(-1.d0/3.d0))*R*cosy(ly,ky) + (2.d0**(1.d0/6.d0))*R*siny(ly,ky)*(-cosx(lx,kx)+sinx(lx,kx))/2.d0
              !Vr(ly,lx,ky,kx,1)=-dd*(dcosh(r12/L)**(-2.0d0))
              !Vr(ly,lx,ky,kx,2)=-dd*(dcosh(r34/L)**(-2.0d0))
              !Vr(ly,lx,ky,kx,3)=-dd*(dcosh(r13/L)**(-2.0d0))
              !Vr(ly,lx,ky,kx,4)=-dd*(dcosh(r14/L)**(-2.0d0))
              !Vr(ly,lx,ky,kx,5)=-dd*(dcosh(r23/L)**(-2.0d0))
              !Vr(ly,lx,ky,kx,6)=-dd*(dcosh(r24/L)**(-2.0d0)) 
              !cc Square well potential
              !                      Vr(ly,lx,ky,kx,1)=0.0d0
              !                      Vr(ly,lx,ky,kx,2)=0.0d0                      
              !                      Vr(ly,lx,ky,kx,3)=0.0d0
              !                      Vr(ly,lx,ky,kx,4)=0.0d0
              !                      Vr(ly,lx,ky,kx,5)=0.0d0
              !                      Vr(ly,lx,ky,kx,6)=0.0d0    
              !                      if (1.0d0 .ge. r12) then
              !                          Vr(ly,lx,ky,kx,1)=-5.0d0
              !                      endif
              !                      if (1.0d0 .ge. r34) then
              !                          Vr(ly,lx,ky,kx,2)=-5.0d0
              !                      endif
              !                      if (1.0d0 .ge. r13) then
              !                          Vr(ly,lx,ky,kx,3)=-5.0d0
              !                      endif
              !                      if (1.0d0 .ge. r14) then
              !                          Vr(ly,lx,ky,kx,4)=-5.0d0
              !                      endif
              !                      if (1.0d0 .ge. r23) then
              !                          Vr(ly,lx,ky,kx,5)=-5.0d0
              !                      endif
              !                      if (1.0d0 .ge. r24) then
              !                          Vr(ly,lx,ky,kx,6)=-5.0d0
              !                      endif                      
           enddo
        enddo
     enddo
  enddo

  do idv = 1,6
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
                    Vmid(NewRow,Col,idv) = 0.0d0
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo

  do idv = 1,6
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
                                VTempInt = VTempInt + xTempV(lx)*yTempV(ly,ky)*Vr(ly,lx,ky,kx,idv)
                             enddo
                          enddo
                          VInt = VInt + xIntScale(kx)*yIntScale(ky)*VTempInt
                       enddo
                    enddo
                    Vmid(NewRow,Col,idv) = Vmid(NewRow,Col,idv)+VInt
                    !                     write(25,*) ix,ixp,H(NewRow,Col)
                 endif
              enddo
           enddo
        enddo
     enddo
  enddo

  do idv = 1,6
     do j = 1,NumStates
        do k = 1,MatrixDim
           TempPsi1(k) = mPsi(k,j)
        enddo
        !call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,Vmid(:,:,idv),HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
        do i = 1,NumStates
           Vexp(i,j,idv) = ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
        enddo
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  deallocate(Vr)
  deallocate(xIntScale,xTempV)
  deallocate(cosy,siny,yIntScale,yTempV)
  deallocate(sinx,cosx)
  deallocate(kxMin,kxMax)
  deallocate(kyMin,kyMax)

  return
end subroutine CalcVMatrix

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
     !call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        P(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcPMatrix

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
     !call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,TempPsi1,1,0.0d0,TempPsi2,1)
     do i = 1,NumStates
        Q(i,j) = a*ddot(MatrixDim,TempPsi2,1,mPsi(1,i),1)
     enddo
  enddo

  deallocate(TempPsi1,TempPsi2)

  return
end subroutine CalcQMatrix

subroutine FixPhase(NumStates,HalfBandWidth,MatrixDim,S,ncv,mPsi,rPsi)

  integer NumStates,HalfBandWidth,MatrixDim,ncv
  double precision S(HalfBandWidth+1,MatrixDim),Psi(MatrixDim,ncv)
  double precision mPsi(MatrixDim,ncv),rPsi(MatrixDim,ncv)
  integer i,j
  double precision Phase,ddot
  double precision, allocatable :: TempPsi(:)
!  write(6,*) "in FixPhase: allocating memory for TempPsi"
  allocate(TempPsi(MatrixDim))

  do i = 1,NumStates
!     write(6,*) 'in FixPhase: i = ', i
     call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rPsi(1,i),1,0.0d0,TempPsi,1)
     Phase = ddot(MatrixDim,mPsi(1,i),1,TempPsi,1)
     if (Phase .lt. 0.0d0) then
        !print *, "needed to fix phase"
        do j = 1,MatrixDim
           rPsi(j,i) = -rPsi(j,i)
        enddo
     endif
  enddo

  deallocate(TempPsi)

  return
end subroutine FixPhase
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
end subroutine GridMaker222



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
end subroutine GridMaker111


subroutine GridMaker(R,r0,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)
  implicit none
  integer xNumPoints,yNumPoints
  double precision R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

  integer i,j,k
  double precision Pi
  double precision r0New
  double precision xRswitch,yRswitch
  double precision xDelt,x0,x1,x2
  double precision yDelt,y0,y1,y2

  Pi = 3.1415926535897932385d0

  x0 = xMin
  x1 = xMax
  k = 1
  xDelt = (x1-x0)/dfloat(xNumPoints-1)

  do i = 1,xNumPoints
     xPoints(k) = (i-1)*xDelt + x0
     k = k + 1
  enddo

15 format(6(1x,1pd12.5))

  y0 = yMin
  y1 = yMax

  k = 1
  yDelt = (y1-y0)/dfloat(yNumPoints-1)
  do i = 1,yNumPoints
     yPoints(k) = (i-1)*yDelt + y0
     k = k + 1
  enddo
  return
end subroutine GridMaker

subroutine GridMakerBetter(R,xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,xPoints,yPoints)

  integer xNumPoints,yNumPoints
  double precision R,r0,xMin,xMax,yMin,yMax,xPoints(xNumPoints),yPoints(yNumPoints)

  integer i,j,OPGRID
  double precision Pi
  double precision r0New
  double precision xRswitch,yRswitch
  double precision xDelt,x0,x1,x2,x3
  double precision yDelt,y0,y1,y2,y3

  Pi = 3.1415926535897932385d0

  x0 = xMin
  x1 = xMax
  !     write(96,*) 'x0,x1=',x0,x1
  k = 1
  xDelt = (x1-x0)/dfloat(xNumPoints-1)
  do i = 1,xNumPoints
     xPoints(k) = (i-1)*xDelt + x0
     !     xPoints(i) = (i-1)*xDelt + x0
     k = k + 1
  enddo
  y0 = yMin
  y1 = yMax
  !     write(96,*) 'y0,y1=',y0,y1
  k = 1
  yDelt = (y1-y0)/dfloat(yNumPoints-1)
  do i = 1,yNumPoints
     yPoints(k) = (i-1)*yDelt + y0
     k = k + 1
  enddo

  !     write(96,15) (yPoints(k),k=1,yNumPoints)

  !
  !     write(96,15) (xPoints(k),k=1,xNumPoints)
  OPGRID=1
  if(OPGRID.eq.1) then
     !     print*, 'R>xRswitch!! using modified grid!!'
     x0 = xMin
     x1=x0+Pi/(2.0d0*8.0d0)
     x2 = xMax-Pi/(2.0d0*8.0d0)
     x3=xMax
     k = 1
     !         write(6,*) ' The phi grid:'
     xDelt = (x1-x0)/dfloat(xNumPoints/4)
     do i = 1,xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x0
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
     xDelt = (x2-x1)/dfloat(xNumPoints/2)
     do i = 1,xNumPoints/2
        xPoints(k) = (i-1)*xDelt + x1
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo
     xDelt = (x3-x2)/dfloat(xNumPoints/4-1)
     do i = 1,xNumPoints/4
        xPoints(k) = (i-1)*xDelt + x2
        !            print*, k, xPoints(k), xDelt
        k = k + 1
     enddo

     do k=1,2
        do i=2,xNumPoints-1
           xPoints(i)=(xPoints(i-1)+2.0d0*xPoints(i)+xPoints(i+1))/4.0d0
        enddo
     enddo

     !         write(6,*) 'the theta grid:'
     y0 = yMin
     y1=y0+Pi/(2.0d0*3.0d0)
     y2=yMax
     k = 1
     yDelt = (y1-y0)/dfloat(yNumPoints/4)
     do i = 1,yNumPoints/4
        yPoints(k) = (i-1)*yDelt + y0
        !            print*, k, yPoints(k), yDelt
        k = k + 1
     enddo
     yDelt = (y2-y1)/dfloat(3*yNumPoints/4-1)
     do i = 1,3*yNumPoints/4
        yPoints(k) = (i-1)*yDelt + y1
        !            print*, k, yPoints(k), yDelt
        k = k + 1
     enddo
     do k=1,2
        do i=2,yNumPoints-1
           yPoints(i)=(yPoints(i-1)+2.0d0*yPoints(i)+yPoints(i+1))/4.0d0
        enddo
     enddo

  endif



15 format(6(1x,1pd12.5))
  return
end subroutine GridMakerBetter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
     !call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rDiffPsi,1,0.0d0,TempPsi,1)
     !call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,rSumPsi,1,0.0d0,TempPsiB,1)
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
end subroutine CalcCoupling
Double Precision Function Vpot95(r)
  implicit real*8(a-h,o-z)
  double precision eps,D,c6,c8,c10,Astar,alpha,beta
  double precision rmin,Fscale,vscale,vlr

  !     convert from angstroms to a.u.
  rmin = 2.96830d0/0.529177249d0

  !     convert from Kelvin to a.u.
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
end Function Vpot95


!     SUBROUTINE potential3(dim, DIMEN, DIMMAX, wlkdmc, vx, vy, vz, 
!     *                      r_ij_old, mass, e_pot)
DOUBLE PRECISION FUNCTION Vpot91(r)
  IMPLICIT REAL*8(A-H,O-Z)

  !CCCC LM2M2 with add-on potential: R. A. Aziz, and M. J. Slaman
  !CCCC J. Chem. Phys. 94, 8047 (1991)

  INTEGER dim, DIMEN, DIMMAX, wlkdmc
  DOUBLE PRECISION mass, e_pot
  INTEGER i, j, k, position

  DOUBLE PRECISION B, alpha, beta, A, epsilo, r_m, r
  DOUBLE PRECISION c_6, c_8, c_10, d, x, fhelp, sum, x_0, x_1

  DOUBLE PRECISION PI, add_on

  PARAMETER(B = 0.0026d0, A = 1.89635353d0 * 10**5)
  PARAMETER(alpha = 10.70203539d0, beta = -1.90740649d0)
  PARAMETER(c_6 = 1.34687065d0, c_8 = 0.41308398d0, c_10 = 0.17060159d0)
  PARAMETER(epsilo = 3.4739577d0 / 10**5, D = 1.4088d0, r_m = 5.6115d0)
  PARAMETER(x_0 = 1.003535949d0, x_1 = 1.454790369d0)
  PARAMETER(PI = 3.14159265358979324d0)
  !     PI = 2.0d0 * ASIN(1.0d0)

  e_pot = 0.0d0      
  !     C      DO j = 1, wlkdmc
  !     C         DO i = 1, dim
  !     C            DO k = i + 1, dim
  !     C               position = (i - 1) * dim - (i - 1) * i / 2 + k - i
  x = r / r_m
  IF(x .LT. D) THEN
     fhelp = exp(-(D / x - 1.0d0) * (D / x - 1.0d0))
  ELSE
     fhelp = 1.0d0
  ENDIF
  sum = c_6 / x**6 + c_8 / x**8 + c_10 / x**(10)
  IF ((x .LE. x_1) .AND. (x .GE. x_0)) THEN
     add_on = B * (SIN(2.0d0 * PI * (x - x_0) / (x_1 - x_0) - 0.5d0 * PI) + 1.0d0) 
  ELSE
     add_on = 0.0d0
  ENDIF
  e_pot = e_pot + (A * exp(-alpha * x + beta * x * x) - fhelp * sum + add_on) * epsilo 
  !     C            END DO
  !     C         END DO
  !     C      END DO
  Vpot91=e_pot
  RETURN
END FUNCTION Vpot91

DOUBLE PRECISION FUNCTION Vpot(r)
  IMPLICIT REAL*8(A-H,O-Z)
  double precision kelvinPERau, angstromPERau, eps, sigma
  data kelvinPERau,angstromPERau/315933.16d0,0.529177249d0/
  data eps,sigma/35.6d0,2.749d0/
  ratio = (sigma/(angstromPERau*r))**6
  Vpot = 4.d0 * (eps/kelvinPERau) * (ratio**2-ratio)
  return
end FUNCTION Vpot

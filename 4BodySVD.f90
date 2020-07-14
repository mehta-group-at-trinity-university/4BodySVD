program SVD
  implicit none
  real*8 n,w,m,Rmin,Rmax,res,temp,overlap,tempA,tempB,epsout
  real*8 Shift,alpha,r0diatom,etaOVERpi,xMin,xMax,yMin,yMax,RDerivDelt,RFirst,RLast,V2Depth
  real*8 e,eMin,eMax,m1,m2,m3,m4,mu4,massarray(4)
  real*8, allocatable :: nodesR(:),weightsR(:)
  real*8, allocatable :: Bsparse(:),eigenVals(:),ress(:)
  real*8, allocatable :: EVals(:),EVecs(:,:),beta(:)
  real*8, allocatable :: dXr(:,:),TmatR(:,:),weightsDMR(:,:),S(:,:),eigenVecs(:,:)
  real*8, allocatable :: iPsi(:,:),jPsi(:,:),H(:,:),B(:,:)
  real*8, allocatable :: Uad(:,:,:),Psi(:,:,:),O(:,:)
  integer LobattoPoints
  integer i,j,k,ivec,ifilenodes,ix,sR,sRc,row,mu,nu,inu,jmu,rCountA,numN0A,rCountB,numN0B,probSize
  integer xDim,yDim,ncv,matrixdim,HalfBandWidth,set
  integer NumStates,PsiFlag,CouplingFlag,LegPoints,order,Left,Right,Bottom,Top,xNumPoints,yNumPoints
  integer fparam(64),info,loop,Tstart,Tend,rate
  integer, allocatable :: Hrow(:),Hcol(:),Brow(:),Bcol(:),indexOf(:,:)
  character*64 LegendreFile
  character*64 Outputfile
  CHARACTER*64 LobattoFile
  character*3 wfnfile
  character*3 wfnsuffix
  character*64 wfnout
  integer iwfn
  double precision, parameter :: Pi = 3.1415926535897932384626433832795
  read(5,*)
  read(5,*) xNumPoints, yNumPoints
  read(5,*)
  read(5,*)
  read(5,*) LobattoPoints, NumStates, Order
  read(5,*)
  read(5,*)
  read(5,*) V2Depth, Rmin, Rmax, alpha
  read(5,*)
  read(5,*)
  read(5,*) m1,m2,m3,m4
  read(5,*)
  read(5,*)
  read(5,*) xMin, xMax, yMin, yMax
  read(5,*)
  read(5,*)
  read(5,*) Left, Right, Bottom, Top

  xMin = xMin*Pi
  xMax = xMax*Pi
  yMin = yMin*Pi
  yMax = yMax*Pi
  write(6,*) "Bounds = ", xMin, xMax, yMin, yMax
  write(6,*) "V2Depth = ", V2Depth

  !mu4=(m1*m2*m3*m4/(m1+m2+m3+m4))**(1.0d0/3.0d0)
  mu4=(m1*m2*m3*m4)**(0.25d0)
  write(6,*) "mu4 = ", mu4
  massarray = (/m1,m2,m3,m4/)
  Outputfile = "Eigenvals.dat"!"eigenVals-100x100x100x20x5-HHHHOdd.dat"
  open(unit=9,file=Outputfile)

  call system_clock(Tstart)

  sR=LobattoPoints
  sRc=sR-2

  allocate(nodesR(sR),weightsR(sR),dXr(sRc,sR),weightsDMR(sR,sR),TmatR(src,src))

  ! read in the Gauss-Lobatto nodes and weights.  Be sure LobattoPoints is equal to one of
  ! the options in the input following file:
  LobattoFile = 'GaussLobatto.dat'
  call GetGaussLobattoFactors(LobattoFile,LobattoPoints,nodesR,weightsR)

  call rescale(sR,Rmin,Rmax,nodesR,weightsR)
 ! write(6,*) "scaled nodes:"
 ! call printmatrix(nodesR,1,sr,6)
 ! write(6,*) "scaled weights:"
 ! call printmatrix(weightsR,1,sr,6)
  do i=1,sRc
     do j=1,sR
        call calcDX(nodesR,i+1,j,sR,weightsR,res)
        dXr(i,j)=res
     end do
  end do

  weightsDMR=0d0
  do i=1,sR
     weightsDMR(i,i)=weightsR(i)!*nodesR(i)**2
  end do
  TmatR=0d0
  TmatR = (0.5d0/mu4)*MATMUL(dXr,MATMUL(weightsDMR,TRANSPOSE(dXr)))
!  write(6,*) "Kinetic energy matrix:"
!  call printmatrix(TmatR,src,src,6)
  
!  call testbasis(nodesR, weightsR, sr)
!  call test1D(TmatR,sr,weightsR,nodesR)
!  stop
  !obtain the adiabatic eigenvalues and eigenvectors 
  LegendreFile='Legendre.dat'
  LegPoints=10

  ! Set the "shift" for the DSBAND solver to -5*depth
  Shift = -5d0*V2Depth

  RDerivDelt=0.0001d0
  RFirst=nodesR(2)
  RLast=nodesR(sR-1)

  !ncv=2*NumStates
  xDim=xNumPoints+Order-3
  yDim=yNumPoints+Order-3+1
  MatrixDim=xDim*yDim
  HalfBandWidth=yDim*Order+Order
  print *, "MatDIM: ", MatrixDim
  probSize=sRc*numStates
  PsiFlag=0

  write(6,*) NUmStates
  write(6,*) HalfBandWidth, MatrixDim
  allocate(Uad(sRc,numStates,2),Psi(sRc,MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim),&
       indexOf(sRc,NumStates),O(probSize,probSize),&
       iPsi(MatrixDim,numStates),jPsi(MatrixDim,numStates))
  call adiabaticSolver(NumStates,PsiFlag,0,LegendreFile,LegPoints,Shift,Order,Left,Right,top,Bottom,alpha,massarray,&
       xNumPoints,xMin,xMax,yNumPoints,yMin,yMax,sRc,RDerivDelt,RFirst,RLast,V2Depth,&
       nodesR(2:sR-1),Uad,Psi,numstates,MatrixDim,S,HalfBandWidth+1)

  write(6,*)"About to calculate the overlap matrix..."
  ! create an array for the collective index
  row=1
  do i=1,sRc
     do j=1,numStates
        indexOf(i,j)=row
        row=row+1
     end do
  end do
  ! Begin calculation of overlap matrix
  O=0d0

  do i=1,sRc
     write(6,*) i
     do j=1,sRc
        iPsi(:,:) = Psi(i,:,:)
        jPsi(:,:) = Psi(j,:,:)
        call CalcOMatrix(Numstates,HalfBandWidth,MatrixDim,iPsi,jPsi,&
             S,O,sRc,i,j,indexOf,probSize)
     end do
  end do

  allocate(H(probsize,probsize),EVals(probsize),EVecs(probsize,probsize))
  H=0d0

  write(6,*) "Constructing the Full Hamiltonian..."
  do i=1,sRc
     do nu=1,numStates
        do j=1,sRc
           do mu=1,numStates
              inu=indexOf(i,nu)
              jmu=indexOf(j,mu)

              !H(inu,jmu) = H(inu,jmu) + TmatR(i,j)*O(inu,jmu)
              H(inu,jmu) = TmatR(i,j)*O(inu,jmu)
              if(inu.eq.jmu) then
                 H(inu,jmu) = H(inu,jmu) + Uad(i,nu,1)
              endif

           end do
        end do
     end do
  end do

!!!!!!
  !Hamiltonian symmetry check
!!$     k=0
!!$     do i=1,sRc
!!$        do nu=1,numStates
!!$           do mu=1,numStates
!!$              inu=indexOf(i,nu)
!!$              jmu=indexOf(i,mu)
!!$              if(O(inu,jmu) .ne. O(inu,jmu) )  then
!!$                 k=k+1
!!$                 write(6,*) O(inu,jmu)
!!$              endif
!!$           end do
!!$        end do
!!$     end do
!!$     write(6,*) "DIFFS"
!!$     do i=1,sRc
!!$        do nu=1,numStates
!!$           !do j=1,sRc
!!$           !do mu=1,numStates
!!$           inu=indexOf(i,nu)
!!$           jmu=indexOf(i+1,nu)
!!$           write(6,*) O(inu,jmu)-O(jmu,inu)
!!$           !if (O(inu,jmu) .ne. O(jmu,inu)) k=k+1
!!$           !end do
!!$           !end do
!!$        end do
!!$     end do
!!$     write(6,*) k
!!$!!!!!

  write(6,*) "About to diagonalize the Hamiltonian. problem size = ", probsize
  call Mydsyev(H,probsize,EVals,EVecs)

  call system_clock(Tend,rate)
  write(6,*) info
  write(6,*) (Tend-Tstart)/rate
  write(9,*) '#',(Tend-Tstart)/rate, numstates, Uad(sRc,1,1)

  do i=1,probsize
     write(9,*) i, EVals(i)
    !write(301,*) (EVecs(i,j),j=1,probsize)
  end do
  wfnfile = "wfn"
  do ivec = 1, 20
     write(wfnsuffix,'(I3)') 100+ivec
     wfnout = wfnfile // "-" // wfnsuffix // ".dat"
     write(6,*) wfnout
     !        open(9999,"wfnout")
     !       do iR = 1, nR
     !          write(file,10) nodesR(iR+1), (evecs(indexOf(iR,nu),ivec)/sqrt(weightsR(iR+1)), nu = 1, NumStates)
     !      enddo
     !        close(9999+ifile)
     call Plotwfn(ivec,indexOf,nodesR,weightsR,sr,src,NumStates,probsize,EVecs,wfnout)

  enddo
10 format(1P,2000e16.8)     
  close(1)
  !close(2)
  !close(3)
  !close(4)
  !close(5)
  !close(6)
  !close(8)
  close(10)
  close(20)
  !close(15)
  !close(100)
  !close(101)
  !close(102)
  !close(200)
  !close(201)
  !close(300)
  !close(301)
  !close(400)
  !close(401)


  deallocate(nodesR,weightsR,dXr,weightsDMR,TmatR,Uad,Psi,S,&
       H,EVals,EVecs,indexOf,O,iPsi,jPsi)



end program SVD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Plotwfn(ivec,indexOf,nodesR,weightsR,ngl,nR,NumStates,probsize,evecs,wfnfile)
  implicit none
  integer ivec,nu,iR,ngl,NumStates,nR,probsize,inu,indexOf(nR,NumStates)
  real*8 nodesR(ngl),weightsR(ngl),evecs(probsize,probsize)
  character*64 wfnfile
  !  write(6,*) file
  open(unit=8,file=wfnfile)
  do iR = 1, nR
     write(8,10) nodesR(iR+1), (evecs(indexOf(iR,nu),ivec)*sqrt(weightsR(iR+1)), nu = 1, NumStates)
  enddo
  close(8)
10 format(1P,2000e16.8)

end subroutine Plotwfn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rescale(n,a,b,nodes,weights)
  implicit none
  integer  n,i
  real*8  nodes(n), tempx(n), tempw(n), weights(n)
  real*8  x1,xN,a,b
  x1= nodes(1)
  xN= nodes(n)
  do i=1,n
     !nodes(i)=a+(nodes(i)-x1)/(xN-x1)*(b-a)
     tempx(i) = 0.5d0*(a+b) + nodes(i)*0.5d0*(b-a)
     tempw(i) = 0.5d0*(b-a)*weights(i)
  end do
  nodes(:)=tempx(:)
  weights(:)=tempw(:)
end subroutine rescale

!find derivitave of the DVR basis at the DVR points
subroutine CalcDX(nodes,i,l,sz,weights,res)
  implicit none
  integer  i,sz,j,k,l
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res
  
  res=0

  do k=1,sz
     holder=0d0
     if(k .ne. i) then
        do j=1,sz
           if(j.ne.i) holder(j)=(nodes(l)-nodes(j))/(nodes(i)-nodes(j))
        end do
        holder(i)=1d0
        holder(k)=1d0

        val=PRODUCT(holder)/((nodes(i)-nodes(k)))
        res = res + val
     end if
  end do
  res=res*(weights(i)**(-.5d0))
end subroutine CalcDX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function CalcXfun(nodes,i,R,sz,weights)
  implicit none
  integer  i,sz,j,k
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res,R
  res=0
  do j=1,i-1
     holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
  end do
  holder(i)=1d0
  do j=i+1,sz
     holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
  end do
  res=PRODUCT(holder)/sqrt(weights(i))
  CalcXfun=res
end function CalcXfun
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision function CalcDXfun(nodes,i,R,sz,weights)
  implicit none
  integer i,sz,j,k
  real*8 nodes(sz),holder(sz),weights(sz)
  real*8 val,res,R
  res=0

  do k=1,sz
     holder=0d0
     if(k .ne. i) then
        do j=1,sz
           if(j.ne.i) holder(j)=(R-nodes(j))/(nodes(i)-nodes(j))
        end do
        holder(i)=1d0
        holder(k)=1d0

        val=PRODUCT(holder)/((nodes(i)-nodes(k)))
        res = res + val
     end if
  end do
  res=res*(weights(i)**(-.5d0))
  CalcDXfun=res
  
end function CalcDXfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine testbasis(nodes, weights, sr)
  implicit none
  integer i,j,sr,src,np
  double precision, allocatable :: R(:)
  double precision nodes(sr), weights(sr)
  double precision, external :: CalcXfun, CalcDXfun
  np=1000
  allocate(R(np))
  do i = 1, np
     R(i) = nodes(1) + (dble(i)-1d0)/(dble(np)-1d0) * (nodes(sr) - nodes(1))
  enddo
  src = sr-2
  do i = 1, np
     write(2000,*) R(i), (CalcXfun(nodes,j,R(i),sr,weights), j=1,sr)
  enddo

  do i = 1, np
     write(2001,*) R(i), (CalcDXfun(nodes,j,R(i),sr,weights), j=1,sr)
  enddo

end subroutine testbasis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test1D(Tmat,sr,weights,nodes)
  implicit none
  integer sr, src, i, j
  double precision Tmat(sr-2,sr-2),V(sr-2,sr-2)
  double precision weights(sr),nodes(sr)
  double precision evals(sr-2), evecs(sr-2,sr-2)

  src = sr-2
  V=0d0
  do i=1,src
     V(i,i) = 6d0/(2d0*nodes(i+1)**2)
  enddo
  call Mydsyev(Tmat+V,src,evals,evecs)
  open(unit = 9, file = "SphericalWell-Evals.dat")
  do i = 1, src
     write(9,*) evals(i)
  enddo
  close(9)
  open(unit=10, file = "SphericalWell-Evecs.dat")
  do i = 1, src
     write(10,*) nodes(i+1), (evecs(i,j)/sqrt(weights(i)), j=1,src)
  enddo
  close(10)
end subroutine test1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcOMatrix(NumStates,HalfBandWidth,MatrixDim,iPsi,jPsi,S,O,sRc,i,j,indexOf,probSize)
  implicit none
  integer NumStates,HalfBandWidth,MatrixDIm,i,j,k,R,sRc,nu,mu
  integer indexOf(sRc,numStates),inu,jmu,probSize
  real*8 ddot
  double precision iPsi(MatrixDim,NumStates),jPsi(MatrixDim,numStates),S(HalfBandWidth+1,MatrixDim)
  double precision O(probSize,probSize),TempPsi(MatrixDim)
  
  do nu=1,numStates
     call dsbmv('U',MatrixDim,Halfbandwidth,1.0d0,S,Halfbandwidth+1,iPsi(:,nu),1,0.0d0,tempPsi,1)
     do mu=1,numStates
        inu=indexOf(i,nu)
        jmu=indexOf(j,mu)
        O(inu,jmu) = ddot(MatrixDim,TempPsi,1,jPsi(1,mu),1)
        !               write(401,*) O(inu,jmu)
     end do
  end do
end subroutine CalcOMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE printmatrix(M,nr,nc,file)
  IMPLICIT NONE
  INTEGER nr,nc,file,j,k
  DOUBLE PRECISION M(nr,nc)
  
  DO j = 1,nr
     WRITE(file,30) (M(j,k), k = 1,nc)
  ENDDO
  
20 FORMAT(1P,100D16.8)
30 FORMAT(100D14.4)
END SUBROUTINE printmatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

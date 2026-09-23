      subroutine CalcEigenErrors(info,iparam,MatrixDim,H,LeadDim,S,HalfBandWidth,NumStates,Psi,Energies,MaxNumStates)

      integer info,iparam(*),MatrixDim,LeadDim,HalfBandWidth,NumStates,MaxNumStates
      double precision H(LeadDim,MatrixDim),S(HalfBandWidth+1,MatrixDim)
      double precision Psi(MatrixDim,MaxNumStates),Energies(MaxNumStates,2)

      integer j
      double precision dnrm2
      double precision, allocatable :: HPsi(:),SPsi(:)

      if ( info .eq. 0) then

       if (iparam(5) .lt. NumStates) write(6,*) 'Not all states found'

c Compute the residual norm: ||  A*x - lambda*x ||

       allocate(HPsi(MatrixDim))
       allocate(SPsi(MatrixDim))
       do j = 1,NumStates
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,H,LeadDim,Psi(1,j),1,0.0d0,HPsi,1)
        call dsbmv('U',MatrixDim,HalfBandWidth,1.0d0,S,HalfBandWidth+1,Psi(1,j),1,0.0d0,SPsi,1)
        call daxpy(MatrixDim,-Energies(j,1),SPsi,1,HPsi,1)
        Energies(j,2) = dnrm2(MatrixDim,HPsi,1)
        Energies(j,2) = Energies(j,2)/dabs(Energies(j,1))
       enddo
       deallocate(HPsi)
       deallocate(SPsi)
      else
       write(6,*) ' '
       write(6,*) ' Error with _sband, info= ', info
       write(6,*) ' Check the documentation of _sband '
       write(6,*) ' '
      end if

      return
      end
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  Matrices A and B are stored in LAPACK-style band form.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are called Ritz values and Ritz vectors respectively.  They are referred 
c  to as such in the comments that follow.  The computed orthonormal basis 
c  for the invariant subspace corresponding to these Ritz values is referred 
c  to as a Lanczos basis.
c
c \Usage
c   call dsband
c      ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, N, AB, MB, LDA, 
c        RFAC, KL, KU, WHICH, BMAT, NEV, TOL, RESID, NCV, V, 
c        LDV, IPARAM, WORKD, WORKL, LWORKL, IWORK, INFO )
c
c \Arguments
c
c  RVEC    Logical (INPUT)
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute the associated Ritz vectors. 
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute all Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is not referenced.
c
c  D       Double precision array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT) 
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c 
c  N       Integer.  (INPUT) 
c          Dimension of the eigenproblem.  
c 
c  AB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix A in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set.
c          The j-th column of A is stored in the j-th column of the
c          array AB as follows:
c          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c
c  MB      Double precision array of dimension LDA by N. (INPUT)
c          The matrix M in band storage, in rows KL+1 to
c          2*KL+KU+1; rows 1 to KL of the array need not be set. 
c          The j-th column of M is stored in the j-th column of the
c          array AB as follows:
c          MB(kl+ku+1+i-j,j) = M(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
c          Not referenced if IPARAM(7) = 1
c
c  LDA     Integer. (INPUT)
c          Leading dimension of AB, MB, RFAC.
c
c  RFAC    Double precision array of LDA by N. (WORKSPACE/OUTPUT)
c          RFAC is used to store the LU factors of MB when IPARAM(7) = 2 
c          is invoked.  It is used to store the LU factors of
c          (A-sigma*M) when IPARAM(7) = 3,4,5 is invoked.
c          It is not referenced when IPARAM(7) = 1.
c
c  KL      Integer. (INPUT)
c          Max(number of subdiagonals of A, number of subdiagonals of M)
c
c  KU      Integer. (OUTPUT)
c          Max(number of superdiagonals of A, number of superdiagonals of M)
c
c  WHICH   Character*2.  (INPUT)
c          When IPARAM(7)= 1 or 2,  WHICH can be set to any one of
c          the following.
c  
c            'LM' -> want the NEV eigenvalues of largest magnitude.
c            'SM' -> want the NEV eigenvalues of smallest magnitude.
c            'LA' -> want the NEV eigenvalues of largest REAL part.
c            'SA' -> want the NEV eigenvalues of smallest REAL part.
c            'BE' -> Compute NEV eigenvalues, half from each end of the 
c                    spectrum.  When NEV is odd, compute one more from 
c                    the high end than from the low end. 
c
c          When IPARAM(7) = 3, 4, or 5,  WHICH should be set to 'LM' only. 
c          
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

c  NEV     Integer. (INPUT)
c          Number of eigenvalues of OP to be computed.
c   
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          Represents the dimension of the Lanczos basis constructed
c          by dsaupd for OP.
c
c  V       Double precision array N by NCV.  (OUTPUT)
c          Upon INPUT: the NCV columns of V contain the Lanczos basis 
c                      vectors as constructed by dsaupd for OP.
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns 
c                       represent the Ritz vectors that span the desired 
c                       invariant subspace.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Lanczos basis vector array V computed by dsaupd. In this case
c          if RVEC=.TRUE., the first NCONV=IPARAM(5) columns of V contain
c          the desired Ritz vectors. 
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: 
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          It is set to 1 in this subroutine.  The user do not need
c          to set this parameter.
c          ------------------------------------------------------------
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced. 
c
c          IPARAM(3) = MXITER
c          On INPUT:  max number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" eigenvalues.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsband for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          Not referenced.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c WORKD    Double precision work array of length at least 3*n. (WORKSPACE)
c
c WORKL    Double precision work array of length LWORKL.  (WORKSPACE)
c
c LWORKL   Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV.
c
c IWORK    Integer array of dimension at least N. (WORKSPACE)
c          Used when IPARAM(7)=2,3,4,5 to store the pivot information in the 
c          factorization of M or (A-SIGMA*M).
c            
c INFO     Integer.  (INPUT/OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 in DSAUPD.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informational error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -13: HOWMNY must be one of 'A' or 'P'
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current
c                   Arnoldi factorization.
c
c\Routines called:
c     dsaupd  ARPACK reverse communication interface routine.
c     dseupd  ARPACK routine that returns Ritz values and (optionally)
c             Ritz vectors.
c     dgbtrf  LAPACK band matrix factorization routine.
c     dgbtrs  LAPACK band linear system solve routine. 
c     dlacpy  LAPACK matrix copy routine.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     ddot    Level 1 BLAS that computes the dot product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dgbmv   Level 2 BLAS that computes the band matrix vector product.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it.
c
c\Author    
c     Danny Sorensen
c     Richard Lehoucq
c     Chao Yang
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
      subroutine MyDsband(select,d,z,ldz,sigma,n,ab,mb,lda,rfac,ldrfac,k,nev,
     >                  tol,resid,ncv,v,ldv,iparam,workd,workl,lworkl,iwork,info)

      character        which*2, bmat, howmny
      integer          n, lda, ldrfac, k, nev, ncv, ldv, ldz, lworkl, info  
      Double precision tol, sigma
      logical          rvec

      integer          iparam(*), iwork(*)
      logical          select(*)
      Double precision d(*), resid(*), v(ldv,*), z(ldz,*), ab(lda,n), mb(lda,n), rfac(ldrfac,n), workd(*), workl(*)

      integer          ipntr(14)

      integer          ido, i, j, Row, Col, type, ierr

      Double precision one, zero
      parameter        (one = 1.0, zero = 0.0)

      Double precision ddot, dnrm2, dlapy2
      external         ddot, dcopy, dgbmv, dgbtrf, dgbtrs, dnrm2, dlapy2, dlacpy

c iparam(3) : Max number of Arnoldi iterations
      iparam(3) = 300
      iparam(7) = 3
      rvec = .TRUE.
      howmny = 'A'
      which = 'LM'
      bmat = 'G'
      type = 4 
      ido = 0
      iparam(1) = 1

      rfac = 0.0d0
      do i = 1,n
       do j = i,min(i+k,n)
        Row = k+1+i-j
        Col = j
        rfac(k+Row,Col) = ab(Row,Col) - sigma*mb(Row,Col)
       enddo
       do j = max(1,i-k),i-1
        Row = 2*k+1
        Col = j
        rfac(Row+i-j,j) = rfac(Row+j-i,i)
       enddo
      enddo

      call dgbtrf(n,n,k,k,rfac,ldrfac,iwork,ierr)
      if ( ierr .ne. 0 )  then
       print*, ' '
       print*, '_SBAND: Error with _gbtrf:',ierr
       print*, ' '
       go to 9000
      end if

  90  continue 

      call dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)

      if (ido .eq. -1) then
       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
       if (ierr .ne. 0) then
        print*, ' ' 
        print*, '_SBAND: Error with _gbtrs.'
        print*, ' ' 
        go to 9000
       end if
      else if (ido .eq. 1) then
       call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)
       call dgbtrs('Notranspose',n,k,k,1,rfac,ldrfac,iwork,workd(ipntr(2)),n,ierr)
       if (ierr .ne. 0) then 
        print*, ' '
        print*, '_SBAND: Error with _gbtrs.' 
        print*, ' '
        go to 9000
       end if
      else if (ido .eq. 2) then
       call dsbmv('U',n,k,1.0d0,mb,lda,workd(ipntr(1)),1,0.0d0,workd(ipntr(2)),1)
      else 
       if ( info .lt. 0) then
        print *, ' '
        print *, ' Error with _saupd info = ',info
        print *, ' Check the documentation of _saupd '
        print *, ' '
        go to 9000
       else 
        if ( info .eq. 1) then
         print *, ' '
         print *, ' Maximum number of iterations reached.'
         print *, ' '
        else if ( info .eq. 3) then
         print *, ' '
         print *, ' No shifts could be applied during 
     >        implicit Arnoldi update, try increasing NCV.'
         print *, ' '
        end if
        if (iparam(5) .gt. 0) then
         call dseupd(rvec,'A',select,d,z,ldz,sigma,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
         if ( info .ne. 0) then
          print *, ' ' 
          print *, ' Error with _neupd = ', info
          print *, ' Check the documentation of _neupd '
          print *, ' ' 
          go to 9000
         endif
        endif
       endif
       go to 9000
      endif

      go to 90 

 9000 continue

      end

   SUBROUTINE DCBSV( N, KU, NRHS, AB, LDAB, B, LDB, INFO ) 
   IMPLICIT NONE
!   Copyright (c) 2021   Anthony M de Beus
!   PURPOSE solves the cyclic/periodic general banded system, see LAPACK routine DGBSV by contrast
!   using an O(N/KU+KU)xKUxKU algorithm 
!   THIS VERSION DOES NOT REQUIRE LINKING WITH LAPACK AND BLAS

    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision 

!  Arguments copied and modified from dgbsv.f *  -- LAPACK routine (version 3.1) --
!      Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!      November 2006
!
!  AB      (input) array, dimension (LDAB,N)
!          On entry, the matrix A in band storage,
!          AB(KU+1+i-j,j) = A(i,j); KU+1+i-j understood "modulo N" so that the
!          subdiagonals wrap around to the last columns and first rows and the 
!          superdiagonals wrap around to the first columns and last rows i.e.
!          AB(2*KU+2-mod(N+KU+1+i-j,N),i) = A(i,j) (transpose of dgbsv.f)
!          See further details below.
!          
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB. LDAB >= 2*KU+1 
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 2*KU+2.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  B       (input/output) array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!  
!  KU      (input) INTEGER
!          Number of superdiagonals not including the central diagonal
!          The number of superdiagonals is equal to the number of subdiagonals.
!          The bandwidth is therefore 2*KU+1  KU>=1
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, a determinant is exactly zero, and the solution
!               has not been computed.
!
!  Further Details
!  ===============
!
!  The band diagonal storage scheme is illustrated by the following example, when
!  N = 9, KU = 2
!
!     AB(2*KU+2-mod(N+KU+1+i-j,N),i)                  A(i,j)=aij  
!     a18  a29  a31  a42  a53 a64 a75 a86 a97         a11 a12 a13  0   0   0   0 a18  a19
!     a19  a21  a32  a43  a54 a65 a76 a87 a98         a21 a22 a23 a24  0   0   0  0   a29
!     a11  a22  a33  a44  a55 a66 a77 a88 a99         a31 a32 a33 a34 a35  0   0  0    0
!     a12  a23  a34  a45  a56 a67 a78 a89 a91         0   a42 a43 a44 a45 a46  0  0    0
!     a13  a24  a35  a46  a57 a68 a79 a81 a92         0   0   a53 a54 a55 a56 a57 0    0
!                                                     0   0   0   a64 a65 a66 a67 a68  0
!                                                     0   0   0   0   a75 a76 a77 a78 a79
!                                                     a81 0   0   0    0  a86 a87 a88 a89
!                                                     a91 a92 0   0    0   0  a97 a98 a99
!
!  .. Scalar Arguments ..
   Integer, Intent(IN) ::  KU, LDAB, LDB, N, NRHS
   INTEGER, INTENT(OUT) :: INFO
!  .. Array Arguments ..
   Real(wp), Intent(IN) :: AB( ldab, * )
   Real(wp), Intent(INOUT) ::  B( ldb, * )

!  .. Work space ..
   REAL(wp) :: Bj(2*KU,N/(2*KU)+2*KU)  
   REAL(wp) :: Cj(2*KU,2*KU,N/(2*KU)+2*KU)
   REAL(wp) :: Pj(2*KU,2*KU,N/(2*KU)+2*KU)
   REAL(wp) :: Sj(2*KU,2*KU,N/(2*KU)+2*KU)   
   REAL(wp) :: BjL(2*KU+mod(N,2*KU))   
   REAL(wp) :: CjL(2*KU+mod(N,2*KU),2*KU+mod(N,2*KU)) 
   REAL(wp) :: UB(2*KU,2*KU,N/(2*KU)+2*KU)   ! Inverse(Cj+Sj*A(:,:,j-1))
   REAL(wp) :: UD(2*KU,2*KU,0:N/(2*KU)+2*KU) ! ud is my set of matrices Aj
   REAL(wp) :: UE(2*KU,0:N/(2*KU)+2*KU)      ! ue is my vectors vj 
   REAL(wp) :: A(2*KU,2*KU),BB(2*KU,2*KU) ! working copies
   REAL(wp) :: IDENTS(2*KU+mod(N,2*KU),2*KU),AL(2*KU,2*KU),BL(2*KU,2*KU+mod(N,2*KU))
   REAL(wp) :: AAL(2*KU+mod(N,2*KU),2*KU+mod(N,2*KU)),CCL(2*KU+mod(N,2*KU))
   INTEGER :: i,j,k,p,ii,jj,kk
    
   p=mod(N,2*KU)
      
!     INFO handling copied/modified from dgbsv.f *  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
      INFO = 0
      IF( N.LT.(2*KU+2) ) THEN
         INFO = -1
      ELSE IF( KU.LT.1 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
!     end info handling
!
!  Initialize
   Bj=0
   Cj=0
   Pj=0
   CjL=0
   BjL=0
   AAL=0
   CCL=0   
   Sj=0   
   jj=0  !index of number of arrays
       
    do j=1,(N-p)/2+1-KU,KU   
       jj=jj+1
         do i=1,KU
          do k=1,KU
           Cj(i,k,jj)=AB(KU+k-i+1,j+i-1)
          end do 
         end do       
         do i=KU+1,2*KU
          do k=KU+1,2*KU
           Cj(i,k,jj)=AB(KU+k-i+1,n-2*KU+i-j+1)
          end do 
         end do
         do i=1,KU
          do k=1,KU
           if ( k <= i ) then
           Pj(i,k,jj)=AB(2*KU+1+k-i,j+i-1)
           endif
          end do 
         end do
         do i=KU+1,2*KU
          do k=KU+1,2*KU
           if ( k >= i ) then
           Pj(i,k,jj)=AB(1+k-i,n-2*KU+i-j+1)
           endif
          end do 
         end do    
         do i=1,KU
          do k=1,KU
           if ( k >= i ) then
           Sj(i,k,jj)=AB(1+k-i,j+i-1)
           endif
          end do 
         end do       
         do i=KU+1,2*KU
          do k=KU+1,2*KU
           if ( k <= i ) then
           Sj(i,k,jj)=AB(2*KU+1+k-i,n-2*KU+i-j+1)
           endif
          end do 
         end do   
   end do

!  LAST ARRAYS needs extra p rows in middle of Cj & Pj now 2KU+p square
         j=(N-p)/2+1-KU ! is found at middle 
         do i=1,2*KU+p
          do k=1,2*KU+p
          if (ABS(k-i) <=  KU) then
           CjL(i,k)=AB(KU+k-i+1,j+i-1)
          endif
          end do 
         end do              
         
!  FIRST EQUATION 
!  initial values
!  A(0) = permuted identity
   UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1
      endif
     end do 
   end do

!  IDENTITY MATRIX WITH p extra rows of zeroes in the middle
   IDENTS=0
    do i=1,KU
     do k=1,KU
      if ( i == k ) then
       IDENTS(i,k)=1
      endif
     end do
    end do
    do i=1+KU+p,2*KU+p
     do k=KU,2*KU
      if ( i == k+p ) then
       IDENTS(i,k)=1
      endif
     end do 
    end do

!  ALL BUT THE LAST EQUATION, GENERATE UD & UB  !   (Cj+Sj*A(:,:,j-1))*A(:,:,j)=-Pj
   jj=0  !index of number of arrays
   do j=1,(N-p)/2-KU,KU     ! j is not used in this loop, it's just a counter 
    jj=jj+1     
    A=Cj(:,:,jj)+matmul(Sj(:,:,jj),UD(:,:,jj-1))
    BB=-Pj(:,:,jj)
    call GaussJordan( 2*KU, 2*KU, A, 2*KU, BB, 2*KU, INFO )
    if (info /= 0) then
     RETURN
    else
     UB(:,:,jj)=A
     UD(:,:,jj)=BB
    endif  
  end do
! DONE WITH ALL BUT LAST INVERSION
! GENERATE SOLUTIONS
  do kk=1,NRHS  
!  FIRST EQUATION V(0)=0
   UE(:,0)=0
!  ALL BUT THE LAST EQUATION GENERATE UE
   jj=0  !index of number of arrays
   do j=1,(N-p)/2-KU,KU          ! j is used here to generate Bj
    jj=jj+1
    do i=1,KU
      Bj(i,jj)=B(j+i-1,kk)   ! write Bj with RHS
    end do     
    do i=KU+1,2*KU
      Bj(i,jj)=B(N-2*KU+i-j+1,kk)
    end do     
!   (Cj+Sj*UD(:,:,j-1))*v(i,j)=bj-Sj*v(:,j-1) using UB=Inverse(Cj+Sj*UD(:,:,j-1)) from above  
    UE(:,jj)=matmul(UB(:,:,jj),Bj(:,jj)-matmul(Sj(:,:,jj),UE(:,jj-1)))    
   end do
!  LAST EQUATION  (last j from above+KU)
    jj=jj+1
    j=(N-p)/2+1-KU
!   remaining values of B centrally
    do i=1,2*KU+p
      BjL(i)=B(j+i-1,kk)   ! write BjL with RHS
    end do      
    CCL=BjL-matmul(IDENTS,matmul(Sj(:,:,jj),UE(:,jj-1)))
    AAL=CjL+matmul(IDENTS,matmul(matmul(Sj(:,:,jj),UD(:,:,jj-1)),Transpose(IDENTS)))    
    call GaussJordan( 2*KU+p, 1, AAL, 2*KU+p, CCL, 2*KU+p, INFO )
    if (info /= 0) then        
     RETURN
    else    
    do i=1,2*KU+p
      B(j+i-1,kk)=CCL(i)         ! overwrite RHS with solution for inner values (reverse of above)
    end do        
    do i=1,KU
      Bj(i,jj)=CCL(i)            ! write Bj with RHS
      end do                     ! but exclude middle p values      
    do i=KU+1,2*KU
      Bj(i,jj)=CCL(i+p)
    end do 
    endif
!   BACKSUBSTITUTION z(j-1)=UE(j-1)+UD(:,:,j-1)*z(j)
    ii=jj ! save jj 
    do while (jj > 1)
      Bj(:,jj-1)=UE(:,jj-1)+matmul(UD(:,:,jj-1),Bj(:,jj))
      jj=jj-1
    end do 
!    overwrite RHS with solution Bj      
    do j=(N-p)/2+1-KU,1,-KU            
     do i=1,KU
      B(j+i-1,kk)=Bj(i,ii)    
     end do   
     do i=KU+1,2*KU
      B(N-2*KU+i-j+1,kk)=Bj(i,ii)     
     end do 
     ii=ii-1      
    end do 
! next kk for NRHS  
  end do  

END SUBROUTINE dcbsv

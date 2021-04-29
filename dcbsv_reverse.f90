   SUBROUTINE DCBSV( N, KU, NRHS, AB, LDAB, B, LDB, INFO )
   Use lapackinterface  
   IMPLICIT NONE
!   Copyright (c) 2021   Anthony M de Beus
!   PURPOSE solves the cyclic/periodic general banded system, see LAPACK routine DGBSV by contrast
!   using an O(N/KU+KU)xKUxKU algorithm 

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
!          Rows of A are Columns of AB; Diagonals of A are Rows of AB.
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
   REAL(wp) :: Bj(2*KU,N/(2*KU)+2*KU,NRHS)  
   REAL(wp) :: Cj(2*KU,2*KU,N/(2*KU)+2*KU)
   REAL(wp) :: Pj(2*KU,2*KU,N/(2*KU)+2*KU)
   REAL(wp) :: Sj(2*KU,2*KU,N/(2*KU)+2*KU)       
   REAL(wp) :: UD(2*KU,2*KU,0:N/(2*KU)+2*KU) ! ud is my set of matrices Aj
   REAL(wp) :: UE(2*KU,0:N/(2*KU)+2*KU,NRHS)      ! ue is my vectors vj 
   REAL(wp) :: A(2*KU,2*KU),AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS) ! working copies
   REAL(wp) :: CCL(2*KU,NRHS)
   REAL(wp), ALLOCATABLE :: CJL(:,:),SPJ(:,:),EEK(:,:) ! pxp and px2*KU, and sometimes p=0
   INTEGER ::  i,j,k,kk,hh,p,ii,jj,ipiv(2*KU)
    
   p=mod(N,2*KU)
   
   if (p /= 0) then 
   write(*,*) 'currently only works for p=mod(N,2*KU)=0'
!   stop
   endif
   
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
         CALL XERBLA( 'DCBSV ', INFO )
         RETURN
      END IF
!     end info handling
!
!  Initialize
   ipiv=0
   Bj=0
   Cj=0
   Pj=0
   CCL=0
   EE=0   
   Sj=0   
   jj=0  !index of number of arrays

!   Setup the arrays             
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
         do i=1,KU
          Bj(i,jj,1:NRHS)=B(j+i-1,1:NRHS)
         end do     
         do i=KU+1,2*KU
          Bj(i,jj,1:NRHS)=B(N-2*KU+i-j+1,1:NRHS)             
         end do                              
   end do

!  LAST ARRAY in reverse uses p=0 formula
   UE(:,0,:)=0
   UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1
      endif
     end do 
   end do        
               

!  FIRST EQUATION if p=0 V(N/(2*KU)=0 & A(N/(2*KU) = permuted identity
   if (p == 0) then
    UE(:,N/(2*KU)+1,:)=0
    UD(:,:,N/(2*KU)+1)=0
    do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,N/(2*KU)+1)=1
      endif
     end do 
    end do
   else
!  if p /= 0 have to compute UD(:,:,N/(2*KU)+1), UE(:,N/(2*KU)+1,:)
    allocate (CJL(p,p),SPJ(p,2*KU),EEK(p,2*KU+NRHS))
!   FIRST ARRAYS now p square and px2*KU
    j=(N-p)/2+1    ! is the start of the middle p values 
    do i=1,p
     do k=1,p
      if (ABS(k-i) <=  KU) then
       CjL(i,k)=AB(KU+k-i+1,j+i-1)
      endif
     end do
!    j=(N-p)/2+1-KU    ! is the start of zj in the middle
     do k=1,KU
      if ( k >= i ) then
       SPJ(i,k)=AB(KU+k-i+1,j-KU+i-1)
      endif
     end do
     do k=KU+1,2*KU
      if ( k >= i ) then
       SPJ(i,k)=AB(2*KU+k-i+1,j-KU+i-1)
      endif
     end do 
    end do
!   solve for UE and precursor of UD
    CJL zpj = -SPJ zj + bj
!   concatenate UD precursor and UD solutions onto EEK        
     EEK(:,1:2*KU)=-SPj(:,:)
     do hh=1,NRHS
      do i=1,2*KU
       EEK(:,2*KU+hh)=B(j+i-1,hh)
      end do     
     end do
!   
    call DGESV(2*KU,2*KU+NRHS,A,2*KU,IPIV,EE,2*KU,INFO) ! overwrites EE into solution 
!    call GaussJordan( 2*KU, 2*KU+NRHS ,A ,2*KU , EE, 2*KU, INFO )   ! overwrites EE into solution                    
    if (info /= 0) then
     CALL XERBLA( 'DGETRS/DCBSV ', INFO )
     RETURN
    else
     UD(:,:,jj)=EE(:,1:2*KU)
    do hh=1,NRHS     
     UE(:,jj,hh)=EE(:,2*KU+hh)                
    end do    
    endif

!   Two cases, p < KU or p >= KU (p <=2*KU always since p=mod(N,2*KU)) 
    if (p < KU) then
!    take px2KU
!    p rows on top and bottom
!    needs A0 in center 2*KU-2*p in size
!    result is 2KUx2KU
     UE(:,(N-p)/(2*KU)+1,:)=0
     UD(:,:,(N-p)/(2*KU)+1)=0
     do i=1,2*KU
      do k=1,2*KU
       if ( (ABS(k - i) - 0  ) == 0 ) then
        UD(i,k,(N-p)/(2*KU)+1)=1
       endif
      end do 
     end do 

     UD(:,:,N/(2*KU)+1)=

    else
!    take px2KU
!    duplicate 2*KU-p middle rows
!    result is 2KUx2KU
     UD(:,:,N/(2*KU)+1)=
    endif
    deallocate (CJL,SPJ)                        
   endif

write(*,*) 'p',p
write(*,*) Transpose(UD(:,:,N/(2*KU)+1))  

       
!  ALL BUT THE LAST EQUATION, GENERATE UD & UE  !   (Cj+Pj*A(:,:,j+1))*A(:,:,j)=-Sj
   jj=N/(2*KU)+1  !index of number of arrays
    do j=N/2,2*KU,-KU            
    jj=jj-1     
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,Pj(:,:,jj),2*KU,UD(:,:,jj+1),2*KU,0.0_wp,AA,2*KU)
    A=Cj(:,:,jj)+AA
!    A=Cj(:,:,jj)+matmul(Pj(:,:,jj),UD(:,:,jj+1))        
!   concatenate Aj and vj solutions onto EE        
     EE(:,1:2*KU)=-Sj(:,:,jj)
     do hh=1,NRHS
      call DGEMV('N',2*KU,2*KU,1.0_wp,Pj(:,:,jj),2*KU,UE(:,jj+1,hh),1,0.0_wp,CC(:,hh),1) 
      EE(:,2*KU+hh)=Bj(:,jj,hh)-CC(:,hh)     
!      EE(:,2*KU+hh)=Bj(:,jj,hh)-matmul(Pj(:,:,jj),UE(:,jj+1,hh))
     end do
!    compute next UD,UE using factored A
    call DGESV(2*KU,2*KU+NRHS,A,2*KU,IPIV,EE,2*KU,INFO) ! overwrites EE into solution 
!    call GaussJordan( 2*KU, 2*KU+NRHS ,A ,2*KU , EE, 2*KU, INFO )   ! overwrites EE into solution                    
    if (info /= 0) then
     CALL XERBLA( 'DGETRS/DCBSV ', INFO )
     RETURN
    else
     UD(:,:,jj)=EE(:,1:2*KU)
    do hh=1,NRHS     
     UE(:,jj,hh)=EE(:,2*KU+hh)                
    end do    
    endif                          
   end do 
        
!  LAST EQUATION in reverse jj=1  so BjL=Bj(:,1,1:NRHS); CjL=Cj(:,:,1)+matmul(Sj(:,:,1),ud(:,:,0))
    call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,Pj(:,:,1),2*KU,UE(:,2,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)       
    CCL(:,1:NRHS)=Bj(:,1,1:NRHS)-CC(:,1:NRHS)
!    CCL(:,1:NRHS)=Bj(:,1,1:NRHS)-matmul(Pj(:,:,1),UE(:,2,1:NRHS))    
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,Pj(:,:,1),2*KU,UD(:,:,2),2*KU,0.0_wp,A,2*KU)
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,Sj(:,:,1),2*KU,UD(:,:,0),2*KU,0.0_wp,AA,2*KU)    
    A=Cj(:,:,1)+AA+A
!    A=Cj(:,:,1)+matmul(Sj(:,:,1),ud(:,:,0))+matmul(Pj(:,:,1),UD(:,:,2))   
    call DGESV(2*KU, NRHS , A, 2*KU, IPIV, CCL(:,1:NRHS), 2*KU, INFO ) ! overwrites CCL 
!    call GaussJordan( 2*KU, NRHS, A ,2*KU, CCL(:,1:NRHS), 2*KU, INFO )  ! overwrites CCL       
    if (info /= 0) then        
     CALL XERBLA( 'DGESV/DCBSV ', INFO )
     RETURN
    else              
    do i=1,2*KU
      Bj(i,1,1:NRHS)=CCL(i,1:NRHS)   ! write Bj with RHS
    end do                                         
    endif    
    
!   BACKSUBSTITUTION z(j+1)=UE(j+1)+UD(:,:,j+1)*z(j) 
    do jj=1,(N-p)/(2*KU)-1
     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UD(:,:,jj+1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
     Bj(:,jj+1,1:NRHS)=UE(:,jj+1,1:NRHS)+CC(:,1:NRHS)
!      Bj(:,jj+1,1:NRHS)=UE(:,jj+1,1:NRHS)+matmul(UD(:,:,jj+1),Bj(:,jj,1:NRHS))
    end do
    
  do kk=1,NRHS         
    ii=1      
!    overwrite RHS with solution Bj      
    do j=1,(N-p)/2+1-KU,KU            
     do i=1,2*KU
      B(j+i-1,kk)=Bj(i,ii,kk)   
      B(N-2*KU+i-j+1,kk)=Bj(i,ii,kk)   
     end do 
     ii=ii+1      
    end do               
  end do    

  
END SUBROUTINE dcbsv

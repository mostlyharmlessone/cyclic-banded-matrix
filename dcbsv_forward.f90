   SUBROUTINE DCBSV_F( N, KU, NRHS, AB, LDAB, B, LDB, INFO )
   Use lapackinterface  
   IMPLICIT NONE
!   Copyright (c) 2021-2023   Anthony M de Beus
!   PURPOSE solves the cyclic/periodic general banded system, see LAPACK routine DGBSV by contrast
!   using an O(N*KU) algorithm 

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
   REAL(wp),ALLOCATABLE :: Bj(:,:,:)  
   REAL(wp),ALLOCATABLE :: Cj(:,:,:)
   REAL(wp),ALLOCATABLE :: Pj(:,:,:)
   REAL(wp),ALLOCATABLE :: Sj(:,:,:)       
   REAL(wp),ALLOCATABLE :: UD(:,:,:)      ! ud is my set of matrices Aj
   REAL(wp),ALLOCATABLE :: UE(:,:,:)      ! ue is my vectors vj         
   REAL(wp),ALLOCATABLE :: AA(:,:),CC(:,:),EE(:,:) ! working copies          
   REAL(wp),ALLOCATABLE :: BjL(:,:),CjL(:,:)      
   REAL(wp),ALLOCATABLE :: IDENTS(:,:),BL(:,:)
   REAL(wp),ALLOCATABLE :: AAL(:,:),CCL(:,:),SjL(:,:)
   Real(wp),ALLOCATABLE :: D(:,:)   
   INTEGER ::  i,j,k,kk,p,ii,jj,LL,allocstat
   INTEGER, ALLOCATABLE :: ipiv(:)
   
   INTERFACE
    subroutine forward_loop(p, L, N, KU, LB, Bj,Cj,Pj,Sj, NRHS, INFO, LU ,UD,UE, LL)
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
      Integer, Intent(IN) :: p, L, LB, LU, KU, N, NRHS  ! L is starting place, LB=size(B,2)=size(C/P/J,3) 
                                                 ! LU+1=size(UD,3)=size(UE,2) index 0 arrays
      INTEGER, INTENT(OUT) :: INFO,LL                 ! LL is number of steps   
      REAL(wp),INTENT(IN) :: Bj(2*KU,LB,NRHS)  
      REAL(wp),INTENT(IN) :: Cj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Pj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Sj(2*KU,2*KU,LB)    
      REAL(wp),INTENT(INOUT) :: UD(2*KU,2*KU,0:LU) ! ud is my set of matrices Aj                 
      REAL(wp),INTENT(INOUT) :: UE(2*KU,0:LU,NRHS)      ! ue is my vectors vj   
     end subroutine forward_loop
     
     SUBROUTINE GaussJordan( N, NRHS, A, LDA, B, LDB, INFO )
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
      INTEGER,INTENT(IN)   :: LDA, LDB, N, NRHS
      INTEGER, INTENT(OUT) :: INFO 
      REAL(wp),INTENT(INOUT) ::  A( LDA, * ), B( LDB, * )
     END SUBROUTINE GaussJordan      
   END INTERFACE      
          
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

   p=mod(N,2*KU)
   
   allocate(Bj(2*KU,(N-p)/(2*KU)+2*KU,NRHS),Cj(2*KU,2*KU,(N-p)/(2*KU)+2*KU))
   allocate(Pj(2*KU,2*KU,(N-p)/(2*KU)+2*KU),Sj(2*KU,2*KU,(N-p)/(2*KU)+2*KU)) 
   allocate(UD(2*KU,2*KU,0:(N-p)/(2*KU)+1),UE(2*KU,0:(N-p)/(2*KU)+1,NRHS),STAT=allocstat)     
    if (allocstat /=0) then
     write(*,*) 'Memory allocation failed'
     stop
    endif    

!  Initialize
   Bj=0;   Cj=0;   Pj=0;  Sj=0   
   jj=0  !index of number of arrays       
!   Setup the arrays       
    do j=1,(N-p)/2+1-KU,KU   
       jj=jj+1                  
         do i=1,KU
           Bj(i,jj,1:NRHS)=B(j+i-1,1:NRHS)
          do k=1,KU
           Cj(i,k,jj)=AB(KU+k-i+1,j+i-1)
          end do
          do k=1,i
            Pj(i,k,jj)=AB(2*KU+1+k-i,j+i-1)
          end do
          do k=i,KU
            Sj(i,k,jj)=AB(1+k-i,j+i-1)
          end do 
         end do               
         do i=KU+1,2*KU
          Bj(i,jj,1:NRHS)=B(n-2*KU+i-j+1,1:NRHS)        
          do k=KU+1,2*KU
           Cj(i,k,jj)=AB(KU+k-i+1,n-2*KU+i-j+1)
          end do
          do k=i,2*KU 
            Pj(i,k,jj)=AB(1+k-i,n-2*KU+i-j+1)
          end do
          do k=KU+1,i 
            Sj(i,k,jj)=AB(2*KU+1+k-i,n-2*KU+i-j+1)
          end do 
         end do                                                     
   end do
     
!  FIRST EQUATION 
!  initial values
!  FIRST EQUATION V(0)=0 & A(0) = permuted identity
   UE(:,0,:)=0
   UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1
      endif
     end do 
   end do

!   main loop function
    call forward_loop(p, 1, N ,KU, size(Bj,2),Bj,Cj,Pj,Sj, NRHS, INFO, size(UD,3)-1,UD, UE, LL) 
    
    deallocate(Cj,Pj)
    if (p /= 0) then
        allocate(D(2*KU+mod(N,2*KU),NRHS))
    endif
    allocate(SjL(2*KU,2*KU))  
!   LAST EQUATION  (last j from above+KU)    
    jj=(N-p)/(2*KU)       
!   j=(N-p)/2+1-KU    
    SjL=Sj(:,:,jj)
    deallocate(Sj) 
    allocate(ipiv(2*KU+mod(N,2*KU)))        
    allocate(BjL(2*KU+mod(N,2*KU),NRHS),CjL(2*KU+mod(N,2*KU),2*KU+mod(N,2*KU)))
    allocate(AAL(2*KU+mod(N,2*KU),2*KU+mod(N,2*KU)),CCL(2*KU+mod(N,2*KU),NRHS))
    allocate(AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS))    
    allocate(IDENTS(2*KU+mod(N,2*KU),2*KU),BL(2*KU,2*KU+mod(N,2*KU)),STAT=allocstat)
    if (allocstat /=0) then
     write(*,*) 'Memory allocation failed'
     stop
    endif 
    
    ipiv=0;   CjL=0;   BjL=0;   AAL=0;   CCL=0;   EE=0
    
!  LAST ARRAYS needs extra p rows in middle of Cj & Pj now 2KU+p square
!       CjL=Cj(:,:,jj)+matmul(Pj(:,:,jj),ud(:,:,0)) when p=0 & jj=N/2
         j=(N-p)/2+1-KU ! is found at middle 
         do i=1,2*KU+p
          do k=1,2*KU+p
          if (ABS(k-i) <=  KU) then
           CjL(i,k)=AB(KU+k-i+1,j+i-1)
          endif
          end do 
         end do   
         
!  IDENTITY MATRIX WITH p extra rows of zeroes in the middle
   IDENTS=0
   forall(j = 1:KU) IDENTS(j,j) = 1    
    do i=1+KU+p,2*KU+p
     do k=KU,2*KU
      if ( i == k+p ) then
       IDENTS(i,k)=1
      endif
     end do 
    end do               
  
!   remaining values of B centrally
    do i=1,2*KU+p
      BjL(i,1:NRHS)=B((N-p)/2-KU+i,1:NRHS)   ! write BjL with RHS
    end do
!    call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,SjL,2*KU,UE(:,jj-1,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)
!    call DGEMM('N','N',2*KU+p,NRHS,2*KU,1.0_wp,IDENTS,2*KU+p,CC(:,1:NRHS),2*KU,0.0_wp,CCL(:,1:NRHS),2*KU+p)        
!    CCL(:,1:NRHS)=BjL(:,1:NRHS)-CCL(:,1:NRHS)
    CCL(:,1:NRHS)=BjL(:,1:NRHS)-matmul(IDENTS,matmul(SjL,UE(:,jj-1,1:NRHS)))    
!    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,SjL,2*KU,UD(:,:,jj-1),2*KU,0.0_wp,AA,2*KU)
!    call DGEMM('N','T',2*KU,2*KU+p,2*KU,1.0_wp,AA,2*KU,IDENTS,2*KU+p,0.0_wp,BL,2*KU)
!    call DGEMM('N','N',2*KU+p,2*KU+p,2*KU,1.0_wp,IDENTS,2*KU+p,BL,2*KU,0.0_wp,AAL,2*KU+p)
!    AAL=CjL+AAL
    AAL=CjL+matmul(IDENTS,matmul(matmul(SjL,UD(:,:,jj-1)),Transpose(IDENTS)))   
    call DGESV(2*KU+p, NRHS , AAL, 2*KU+p, IPIV, CCL(:,1:NRHS), 2*KU+p, INFO ) ! overwrites CCL
!    call GaussJordan( 2*KU+p, NRHS, AAL ,2*KU+p, CCL(:,1:NRHS), 2*KU+p, INFO )  ! overwrites CCL        
    if (info /= 0) then        
     CALL XERBLA( 'DGESV/DCBSV ', INFO )
     RETURN
    else
    if (p /=0 ) then                     
     do i=1,2*KU+p
      D(i,1:NRHS)=CCL(i,1:NRHS)          ! store RHS with solution for inner values      
     end do
    endif        
    do i=1,KU
      Bj(i,jj,1:NRHS)=CCL(i,1:NRHS)      ! write Bj with RHS
      end do                             ! but exclude middle p values      
    do i=KU+1,2*KU
      Bj(i,jj,1:NRHS)=CCL(i+p,1:NRHS)
    end do 
    endif
    
!   BACKSUBSTITUTION z(j-1)=UE(j-1)+UD(:,:,j-1)*z(j) 
    do jj=(N-p)/(2*KU),2,-1
!     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UD(:,:,jj-1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
!     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+CC(:,1:NRHS)
     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+matmul(UD(:,:,jj-1),Bj(:,jj,1:NRHS))
    end do
     
  do kk=1,NRHS         
    ii=(N-p)/(2*KU)      
!    overwrite RHS with solution Bj      
    do j=(N-p)/2+1-KU,1,-KU            
     do i=1,KU
      B(j+i-1,kk)=Bj(i,ii,kk)   
     end do   
     do i=KU+1,2*KU
      B(N-2*KU+i-j+1,kk)=Bj(i,ii,kk)   
     end do 
     ii=ii-1      
    end do               
  end do
! write saved central values
  if (p /= 0) then    
   do i=1,2*KU+p
     B((N-p)/2-KU+i,1:NRHS)=D(i,1:NRHS)
   end do
   deallocate(D)
  endif
  deallocate(SjL,BjL,CjL,AA,ipiv,AAL,CCL,CC,EE,BL,IDENTS)
  
END SUBROUTINE dcbsv_f




         
 






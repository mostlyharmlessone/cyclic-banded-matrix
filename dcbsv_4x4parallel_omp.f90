   SUBROUTINE DCBSV_4( N, KU, NRHS, AB, LDAB, B, LDB, INFO )
   USE lapackinterface 
   USE AB_matrix_fct
   USE OMP_LIB 
        
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!   Copyright (c) 2021   Anthony M de Beus
!   PURPOSE solves the cyclic/periodic general banded system, see LAPACK routine DGBSV by contrast
!   using an O(N*KU) algorithm 

!    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision 

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
   Real(wp), Intent(IN) :: AB( ldab, N)      !ldab==2*KU+1
   Real(wp), Intent(INOUT) ::  B( ldb, NRHS )

!  .. Work space ..

   Real(wp),ALLOCATABLE :: ABP(:,:),BP(:,:)        !paddable version
   Real(wp),ALLOCATABLE :: ABB(:,:),BB(:,:)        !rotated version
   REAL(wp),ALLOCATABLE :: Bj(:,:,:)
   REAL(wp),ALLOCATABLE :: BBj(:,:,:)              !rotated version     
   REAL(wp),ALLOCATABLE :: Cj(:,:,:)
   REAL(wp),ALLOCATABLE :: Pj(:,:,:)
   REAL(wp),ALLOCATABLE :: Sj(:,:,:) 
   REAL(wp),ALLOCATABLE :: CjR(:,:,:)     !rotated version  
   REAL(wp),ALLOCATABLE :: PjR(:,:,:)     !rotated version  
   REAL(wp),ALLOCATABLE :: SjR(:,:,:)     !rotated version
        
   REAL(wp),ALLOCATABLE :: UD(:,:,:)      ! ud is my set of matrices Aj
   REAL(wp),ALLOCATABLE :: UE(:,:,:)      ! ue is my vectors vj 
   REAL(wp),ALLOCATABLE :: UDR(:,:,:)     ! udr is my set of matrices Ajbar
   REAL(wp),ALLOCATABLE :: UER(:,:,:)     ! uer is my vectors vjbar
   REAL(wp),ALLOCATABLE :: VD(:,:,:)      ! vd is my set of matrices Bj
   REAL(wp),ALLOCATABLE :: VE(:,:,:)      ! ve is my vectors uj 
   REAL(wp),ALLOCATABLE :: VDR(:,:,:)     ! vdr is my set of matrices Bjbar
   REAL(wp),ALLOCATABLE :: VER(:,:,:)     ! ver is my vectors ujbar    
   REAL(wp),ALLOCATABLE :: ACL(:,:),CCL(:,:),IDENT(:,:),CC(:,:),JEXC(:,:),Z2WK(:,:),Z2WL(:,:)  
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER ::  i,j,k,kk,nn,p,ii,jj,L1,L2,LL,thread,allocstat,m
   
   INTERFACE
    subroutine forward_loop(L, N, KU, LB, Bj,Cj,Pj,Sj, NRHS, INFO, LU ,UD,UE, LL)
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
      Integer, Intent(IN) ::  L, LB, LU, KU, N, NRHS  ! L is starting place, LB=size(B,2)=size(C/P/J,3) 
                                                 ! LU+1=size(UD,3)=size(UE,2) index 0 arrays
      INTEGER, INTENT(OUT) :: INFO,LL                 ! LL is number of steps   
      REAL(wp),INTENT(IN) :: Bj(2*KU,LB,NRHS)  
      REAL(wp),INTENT(IN) :: Cj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Pj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Sj(2*KU,2*KU,LB)    
      REAL(wp),INTENT(INOUT) :: UD(2*KU,2*KU,0:LU) ! ud is my set of matrices Aj                 
      REAL(wp),INTENT(INOUT) :: UE(2*KU,0:LU,NRHS)      ! ue is my vectors vj   
     end subroutine forward_loop
        
     subroutine backward_loop( L, N, KU, LB,Bj,Cj,Pj,Sj, NRHS, INFO, LR,LU,UDR, UER, LL)
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!  .. Scalar Arguments ..
      Integer, Intent(IN) ::  L, LB, LR, LU, KU, N, NRHS  ! L is starting place, LB=size(B,2)=size(C/P/J,3) 
                                                       ! LU+1=size(UD,3)=size(UE,2) index LR arrays ie. LR=0 for dcbsv_reverse
      INTEGER, INTENT(OUT) :: INFO,LL                     ! LL is number of steps
      REAL(wp),INTENT(IN) :: Bj(2*KU,LB,NRHS)  
      REAL(wp),INTENT(IN) :: Cj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Pj(2*KU,2*KU,LB)
      REAL(wp),INTENT(IN) :: Sj(2*KU,2*KU,LB)    
      REAL(wp),INTENT(INOUT) :: UDR(2*KU,2*KU,LR:LU) ! ud is my set of matrices Aj                 
      REAL(wp),INTENT(INOUT) :: UER(2*KU,LR:LU,NRHS)      ! ue is my vectors vj
     end subroutine backward_loop

     SUBROUTINE GaussJordan( N, NRHS, A, LDA, B, LDB, INFO )
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
      INTEGER,INTENT(IN)   :: LDA, LDB, N, NRHS
      INTEGER, INTENT(OUT) :: INFO 
      REAL(wp),INTENT(INOUT) ::  A( LDA, * ), B( LDB, * )
     END SUBROUTINE GaussJordan      
   END INTERFACE   
      
   p=mod(N,8*KU)
   if (p > 0 ) then
    M=N-p+8*KU
   else
    M=N
   endif
!  pad the matrix out so p=0 formula works and L1 is an integer    
   L1=M/(8*KU)    
   L2=M/(2*KU)+1 

   allocate(ABP(ldab,M),BP(M,NRHS))
   allocate(ABB(ldab,M))      
   allocate(BB(M,NRHS))
   allocate(Bj(2*KU,L2,NRHS),Cj(2*KU,2*KU,L2))
   allocate(BBj(2*KU,L2,NRHS),CjR(2*KU,2*KU,L2))   
   allocate(Pj(2*KU,2*KU,L2),Sj(2*KU,2*KU,L2))
   allocate(PjR(2*KU,2*KU,L2),SjR(2*KU,2*KU,L2))
      
   allocate(UD(2*KU,2*KU,0:L1+1),UE(2*KU,0:L1+1,NRHS))  
   allocate(UDR(2*KU,2*KU,L2-L1-1:L2),UER(2*KU,L2-L1-1:L2,NRHS))   
   
   allocate(VD(2*KU,2*KU,0:L1+1),VE(2*KU,0:L1+1,NRHS))   
   allocate(VDR(2*KU,2*KU,L2-L1-1:L2),VER(2*KU,L2-L1-1:L2,NRHS)) 
  
   allocate(IPIV(8*KU),STAT=allocstat)  ! need IPIV 8*KU not 2*KU for last equations
   if (allocstat /=0) then
    write(*,*) 'Memory allocation failed'
    stop
   endif 
      
!     INFO handling copied/modified from dgbsv.f *  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
      INFO = 0
      IF( N.LT.(2*KU+2) ) THEN
         INFO = -1
      ELSE IF( KU.LT.1 ) THEN
         INFO = -2
      ELSE IF( L2 .LT. 9 ) THEN  !too few blocks for 4 threads
         INFO = -4         
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DCBSV ', INFO )
!         RETURN
      END IF
!     end info handling
!
!  Initialize
   ipiv=0;   Bj=0;   Cj=0;   Pj=0;   Sj=0 ; PjR = 0 ; CjR = 0 ; SjR = 0 ; BBj = 0 
   UER=0;   UDR=0;   UE=0;   UD=0  ; VER=0;   VDR=0;   VE=0;   VD=0 ; ABP=0 ; BP=0 
   jj=0  !index of number of arrays, replaced by loop invariant (j-1)/KU+1 below to enable parallelism
   ii=0 ; kk=0 ; ll=0 

!  pad the data
   ABP(1:ldab,1:N)=AB(1:ldab,1:N)
   BP(1:N,1:NRHS)=B(1:N,1:NRHS)   
   if ( p > 0) then
    ABP(:,N+1:M)=AB(:,1:M-N) !or ABP(KU+1,N+1:M)=1
    BP(N+1:M,:)=B(1:M-N,:)   !BP(N+1:M,:)=0
   endif
   
!  Rotate solution vector and make a copy in BB (multiply B by Asub{n}sup{k})
    k=M/4
    BB(1:M,1:NRHS)=RotateColumns(BP,k)
    BB(1:M,1:NRHS)=ReverseColumns(BB)
    ABB(1:ldab,1:M)=RotateRows(ABP,k)
    ABB(1:ldab,1:M)=ReverseRows(ABB)
    ABB(1:ldab,1:M)=ReverseColumns(ABB)
    
!   Setup the arrays  
!$OMP PARALLEL DO PRIVATE(i,j,k)  
    do j=1,M/2+1-KU,KU             
         do i=1,KU
           Bj(i,(j-1)/KU+1,1:NRHS)=BP(j+i-1,1:NRHS)
           BBj(i,(j-1)/KU+1,1:NRHS)=BB(j+i-1,1:NRHS)
          do k=1,KU
           Cj(i,k,(j-1)/KU+1)=ABP(KU+k-i+1,j+i-1)
           CjR(i,k,(j-1)/KU+1)=ABB(KU+k-i+1,j+i-1)
          end do
          do k=1,i
            Pj(i,k,(j-1)/KU+1)=ABP(2*KU+1+k-i,j+i-1)
            PjR(i,k,(j-1)/KU+1)=ABB(2*KU+1+k-i,j+i-1)
          end do
          do k=i,KU
            Sj(i,k,(j-1)/KU+1)=ABP(1+k-i,j+i-1)
            SjR(i,k,(j-1)/KU+1)=ABB(1+k-i,j+i-1)
          end do 
         end do                        
         do i=KU+1,2*KU
          Bj(i,(j-1)/KU+1,1:NRHS)=BP(m-2*KU+i-j+1,1:NRHS)
          BBj(i,(j-1)/KU+1,1:NRHS)=BB(m-2*KU+i-j+1,1:NRHS)      
          do k=KU+1,2*KU
           Cj(i,k,(j-1)/KU+1)=ABP(KU+k-i+1,m-2*KU+i-j+1)
           CjR(i,k,(j-1)/KU+1)=ABB(KU+k-i+1,m-2*KU+i-j+1)
          end do
          do k=i,2*KU 
            Pj(i,k,(j-1)/KU+1)=ABP(1+k-i,m-2*KU+i-j+1)
            PjR(i,k,(j-1)/KU+1)=ABB(1+k-i,m-2*KU+i-j+1)
          end do
          do k=KU+1,i 
            Sj(i,k,(j-1)/KU+1)=ABP(2*KU+1+k-i,m-2*KU+i-j+1)
            SjR(i,k,(j-1)/KU+1)=ABB(2*KU+1+k-i,m-2*KU+i-j+1)
          end do 
         end do                                                              
   end do
!$OMP END PARALLEL DO  
   deallocate(ABB,ABP)

   nn=2*KU
!  IDENTITY MATRIX
   allocate(IDENT(nn,nn))
   IDENT=0 
   forall(j = 1:nn) IDENT(j,j) = 1
   nn=KU 
!  EXCHANGE MATRIX
   allocate(JEXC(nn,nn))
   JEXC=0
   forall(j = 1:nn) JEXC(nn-j+1,j) = 1 
!  ROTATION MATRIX
!   nn=2*KU
!   allocate(ANK(nn,nn))
!   kk=KU
!   ANK=0
!   forall(j = 1:kk) ANK(j,j+kk) = 1
!   forall(j = 1:nn-kk) ANK(j+kk,j) = 1
!  Z to W transform matrices
   nn=4*KU
   allocate(Z2WK(nn,2*nn))
   Z2WK=0
   Z2WK(1:KU,3*KU+1:4*KU)=JEXC
   Z2WK(KU+1:2*KU,5*KU+1:6*KU)=JEXC
   Z2WK(2*KU+1:3*KU,7*KU+1:8*KU)=JEXC
   Z2WK(3*KU+1:4*KU,KU+1:2*KU)=JEXC
   allocate(Z2WL(nn,2*nn))
   Z2WL=0
   Z2WL(1:KU,6*KU+1:7*KU)=JEXC
   Z2WL(KU+1:2*KU,1:KU)=JEXC
   Z2WL(2*KU+1:3*KU,2*KU+1:3*KU)=JEXC
   Z2WL(3*KU+1:4*KU,4*KU+1:5*KU)=JEXC           

!  ONLY use p=0 formula because use padded AB/ABB
!  already done above initializing UE(:,0,:)=0;  UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1 
        UDR(i,k,L2)=1         
        VD(i,k,0)=1         
        VDR(i,k,L2)=1        
      endif
     end do 
   end do        

!$OMP PARALLEL num_threads(4) PRIVATE(thread)
!  separate iterative parts into subroutines so each thread can work in parallel
   thread = omp_get_thread_num()
   if (thread==0) then
    call forward_loop(L2-L1-1, M ,KU, L2,Bj,Cj,Pj,Sj, NRHS, INFO, L1+1, UD, UE, II)             
   endif
   if (thread==1) then  
    call backward_loop(L2-L1-1, M ,KU, L2,Bj,Cj,Pj,Sj, NRHS, INFO, L2-L1-1, L2,UDR, UER, JJ)   
   endif
   if (thread==2) then 
    call forward_loop(L2-L1-1,M,KU,L2,BBj,CjR,PjR,SjR, NRHS, INFO, L1+1, VD, VE, LL)   
   endif
   if (thread==3) then 
    call backward_loop(L2-L1-1,M,KU,L2,BBj,CjR,PjR,SjR, NRHS, INFO, L2-L1-1, L2 ,VDR, VER, KK)    
   endif
!$OMP END PARALLEL

write(*,*) ii,jj,kk,ll,m

   deallocate(Cj,Pj,Sj,CjR,PjR,SjR)  
   Bj = 0 ; BBj = 0 ! done with input, ready for results
!  LAST EQUATIONS
   allocate(ACL(8*KU,8*KU),CCL(8*KU,NRHS),CC(8*KU,NRHS))
   ACL = 0 ; CCL = 0 ; CC =0
     
   ACL(1:2*KU,1:2*KU)=udR(:,:,JJ)
   ACL(1:2*KU,4*KU+1:6*KU)=-IDENT  
   ACL(2*KU+1:4*KU,2*KU+1:4*KU)=vdR(:,:,KK)
   ACL(2*KU+1:4*KU,4*KU+1:6*KU)=-IDENT
   ACL(4*KU+1:6*KU,2*KU+1:4*KU)=IDENT
   ACL(4*KU+1:6*KU,4*KU+1:6*KU)=-vd(:,:,LL-1)
   ACL(2*KU+1:4*KU,1:8*KU)=matmul(ACL(2*KU+1:4*KU,2*KU+1:6*KU),Z2WK)
   ACL(4*KU+1:6*KU,1:8*KU)=matmul(ACL(4*KU+1:6*KU,2*KU+1:6*KU),Z2WL)              
   ACL(6*KU+1:8*KU,2*KU+1:4*KU)=IDENT
   ACL(6*KU+1:8*KU,6*KU+1:8*KU)=-ud(:,:,II-1)
       
!  RHS   
   CCL(1:2*KU,1:NRHS)=-uER(:,JJ,1:NRHS)      !ue, ueR are vectors v
   CCL(2*KU+1:4*KU,1:NRHS)=-vER(:,KK,1:NRHS) !ve, veR are vectors u
   CCL(4*KU+1:6*KU,1:NRHS)=vE(:,LL-1,1:NRHS)
   CCL(6*KU+1:8*KU,1:NRHS)=uE(:,II-1,1:NRHS)    

   call DGESV(8*KU, NRHS , ACL, 8*KU, IPIV, CCL(:,1:NRHS), 8*KU, INFO ) ! overwrites CCL 
!   call GaussJordan(8*KU, NRHS, ACL ,8*KU, CCL(:,1:NRHS), 8*KU, INFO )  ! overwrites CCL

write(*,*) CCL(:,1)

stop


    if (info /= 0) then     
     CALL XERBLA( 'DGESV ', INFO )
!     CALL XERBLA( 'GAUSSJ ', INFO )     
!     RETURN
    else
      Bj(:,JJ-1,1:NRHS)=CCL(1:2*KU,1:NRHS)   ! write Bjj-1 with RHS zJJ-1      
      Bj(:,II-1,1:NRHS)=CCL(2*KU+1:4*KU,1:NRHS)   ! write Bii-1 with RHS zII-1      
      Bj(:,JJ,1:NRHS)=CCL(4*KU+1:6*KU,1:NRHS)   ! write Bjj with RHS zJJ      
      Bj(:,II,1:NRHS)=CCL(6*KU+1:8*KU,1:NRHS)   ! write Bii with RHS zII
!     now generate w
      CC(1:4*KU,1:NRHS)=matmul(Z2WK,CCL(:,1:NRHS)) 
      CCL(1:4*KU,1:NRHS)=matmul(Z2WL,CCL(:,1:NRHS))          
      BBj(:,KK-1,1:NRHS)=CC(1:2*KU,1:NRHS)        
      BBj(:,LL-1,1:NRHS)=CCL(1:2*KU,1:NRHS)         
      BBj(:,KK,1:NRHS)=CC(2*KU+1:4*KU,1:NRHS)        
      BBj(:,LL,1:NRHS)=CCL(2*KU+1:4*KU,1:NRHS)                 
    endif

    deallocate(ACL,CCL,CC)
    deallocate(IDENT,IPIV,Z2WK,Z2WL,JEXC)    
       
    allocate(CC(2*KU,NRHS)) 
    CC = 0  
       
!   BACKSUBSTITUTION z(j-1)=UE(j-1)+UD(:,:,j-1)*z(j) 
    do jj=LL-1,2,-1
!     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UD(:,:,jj-1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
!     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+CC(:,1:NRHS)
     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+matmul(UD(:,:,jj-1),Bj(:,jj,1:NRHS))
    end do
          
!   BACKSUBSTITUTION z(j+1)=UER(j+1)+UDR(:,:,j+1)*z(j) 
    do jj=KK,M/(2*KU)-1
!     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UDR(:,:,jj+1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
!     Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+CC(:,1:NRHS)
     Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+matmul(UDR(:,:,jj+1),Bj(:,jj,1:NRHS))
    end do    

!   BACKSUBSTITUTION w(j-1)=VE(j-1)+VD(:,:,j-1)*w(j) 
    do jj=LL-1,2,-1
!     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,VD(:,:,jj-1),2*KU,BBj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
!     BBj(:,jj-1,1:NRHS)=VE(:,jj-1,1:NRHS)+CC(:,1:NRHS)
     BBj(:,jj-1,1:NRHS)=VE(:,jj-1,1:NRHS)+matmul(VD(:,:,jj-1),BBj(:,jj,1:NRHS))
    end do

!   BACKSUBSTITUTION w(j+1)=VER(j+1)+VDR(:,:,j+1)*w(j) 
    do jj=KK,M/(2*KU)-1
!     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,VDR(:,:,jj+1),2*KU,BBj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
!     BBj(:,jj+1,1:NRHS)=VER(:,jj+1,1:NRHS)+CC(:,1:NRHS)
     BBj(:,jj+1,1:NRHS)=VER(:,jj+1,1:NRHS)+matmul(VDR(:,:,jj+1),BBj(:,jj,1:NRHS))
    end do 
    
    deallocate(CC,UD,UE,UDR,UER,VD,VE,VDR,VER) 

!   overwrite RHS with solution Bj and BBj, direction irrelevant  
!   Zero out B,BB
    BP(1:M,1:NRHS)=0 ; BB(1:M,1:NRHS)=0
    do jj=1,NRHS                       
     ii=M/(2*KU)+1          
     do j=M/2+1,1,-KU            
      do i=1,KU
        BP(j+i-1,jj)=Bj(i,ii,jj)
        if (ii .ne. LL .and. ii .ne. ll-1 .and. ii .ne. kk .and. ii .ne. kk-1) then ! these are the solved duplicates 
         BB(j+i-1,jj)=BBj(i,ii,jj)                                                  ! from central equation above
        endif
      end do   
      do i=KU+1,2*KU
       BP(M-2*KU+i-j+1,jj)=Bj(i,ii,jj)
       if (ii .ne. LL .and. ii .ne. ll-1 .and. ii .ne. kk .and. ii .ne. kk-1) then
        BB(M-2*KU+i-j+1,jj)=BBj(i,ii,jj)
       endif
      end do
     ii=ii-1      
     end do        
    end do  ! end jj

!  Rotate solution vector back (multiply BB by Asub{n}sup{n-k})
    k=M/4 !(n-k)
    BB(1:M,1:NRHS)=RotateColumns(BB,k)
    BB(1:M,1:NRHS)=ReverseColumns(BB)
!   Replace the missing values into B
    BP(1:M,1:NRHS)=B(1:M,1:NRHS)+BB(1:M,1:NRHS)  
    B(1:N,1:NRHS)=BP(1:N,1:NRHS)
   deallocate (BB,BP,Bj,BBj)
        
  END SUBROUTINE dcbsv_4
  








   SUBROUTINE DCBSV( N, KU, NRHS, AB, LDAB, B, LDB, INFO )
   Use lapackinterface 
   USE OMP_LIB      
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
   REAL(wp),ALLOCATABLE :: Bj(:,:,:)  
   REAL(wp),ALLOCATABLE :: Cj(:,:,:)
   REAL(wp),ALLOCATABLE :: Pj(:,:,:)
   REAL(wp),ALLOCATABLE :: Sj(:,:,:)       
   REAL(wp),ALLOCATABLE :: UD(:,:,:)      ! ud is my set of matrices Aj
   REAL(wp),ALLOCATABLE :: UE(:,:,:)      ! ue is my vectors vj 
   REAL(wp),ALLOCATABLE :: UDR(:,:,:)     ! udr is my set of matrices Ajbar
   REAL(wp),ALLOCATABLE :: UER(:,:,:)     ! uer is my vectors vjbar 
   REAL(wp),ALLOCATABLE :: A(:,:),AA(:,:),CC(:,:),EE(:,:) ! working copies
   REAL(wp),ALLOCATABLE :: CCL(:,:),IDENT(:,:)
   REAL(wp), ALLOCATABLE :: Cp(:,:),Sp(:,:),Bp(:,:),EEK(:,:) ! pxp and px2*KU, and sometimes p=0
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER ::  i,j,k,kk,hh,p,ii,jj,L,LL,thread
    
   p=mod(N,2*KU)
   L=(N-p)/4

   allocate(Bj(2*KU,N/(2*KU)+1,NRHS),Cj(2*KU,2*KU,N/(2*KU)+1))
   allocate(Pj(2*KU,2*KU,N/(2*KU)+1),Sj(2*KU,2*KU,N/(2*KU)+1))
   allocate(UD(2*KU,2*KU,0:N/(4*KU)+1),UE(2*KU,0:N/(4*KU)+1,NRHS))
   allocate(UDR(2*KU,2*KU,N/(4*KU):N/(2*KU)+1),UER(2*KU,N/(4*KU):N/(2*KU)+1,NRHS))
   allocate(A(2*KU,2*KU),AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS))
   allocate(CCL(2*KU,NRHS),IDENT(2*KU,2*KU),IPIV(2*KU))
   if (p /= 0) then
    allocate (Cp(p,p),Sp(p,2*KU),Bp(p,NRHS),EEK(p,2*KU+NRHS))
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
!         RETURN
      END IF
!     end info handling
!
!  Initialize
   ipiv=0;   Bj=0;   Cj=0;   Pj=0;   CCL=0;   EE=0;    Sj=0 
   UER=0;   UDR=0;   UE=0;   UD=0  
   jj=0  !index of number of arrays

!   Setup the arrays             
    do j=1,(N-p)/2+1-KU,KU   
       jj=jj+1
         do i=1,KU
           Bj(i,jj,1:NRHS)=B(j+i-1,1:NRHS)
          do k=1,KU
           Cj(i,k,jj)=AB(KU+k-i+1,j+i-1)
           if ( k <= i ) then
            Pj(i,k,jj)=AB(2*KU+1+k-i,j+i-1)
           endif 
           if ( k >= i ) then
            Sj(i,k,jj)=AB(1+k-i,j+i-1)
           endif                    
          end do 
         end do       
         do i=KU+1,2*KU
          Bj(i,jj,1:NRHS)=B(n-2*KU+i-j+1,1:NRHS)         
          do k=KU+1,2*KU
           Cj(i,k,jj)=AB(KU+k-i+1,n-2*KU+i-j+1)
           if ( k >= i ) then
            Pj(i,k,jj)=AB(1+k-i,n-2*KU+i-j+1)
           endif 
           if ( k <= i ) then
            Sj(i,k,jj)=AB(2*KU+1+k-i,n-2*KU+i-j+1)
           endif                     
          end do 
         end do                              
   end do

!  IDENTITY MATRIX 
   IDENT=0
    do i=1,2*KU
     do k=1,2*KU
      if ( i == k ) then
       IDENT(i,k)=1
      endif
     end do
    end do 

!  FIRST ARRAY forward uses p=0 formula
!  already done above UE(:,0,:)=0;  UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1        
      endif
     end do 
   end do        
               
!  FIRST EQUATION in reverse if p=0 V(N/(2*KU)=0 & A(N/(2*KU) = permuted identity
   if (p == 0) then
!  already done above UER(:,N/(2*KU)+1,:)=0;    UDR(:,:,N/(2*KU)+1)=0
    do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UDR(i,k,N/(2*KU)+1)=1
      endif
     end do 
    end do
   else
   
!  if p /= 0 have to compute UD(:,:,N/(2*KU)+1), UE(:,N/(2*KU)+1,:)
!  already done above UER(:,(N-p)/(2*KU)+1,:)=0;    UDR(:,:,(N-p)/(2*KU)+1)=0
     if (p < KU) then
!    needs A0 in center 2(KU-p) in size
!    result is 2KUx2KU
!     sqsize=2*KU-2*p
      do i=p+1,2*KU-p
       do k=p+1,2*KU-p
        if ( ABS(k - i) == (KU-p) ) then
         UDR(i,k,(N-p)/(2*KU)+1)=1
        endif
       end do 
      end do 
     endif
!  FIRST ARRAYS now p square and px2*KU at the middle values
   Cp=0
   Sp=0
   EEK=0
   j=(N-p)/2+1   ! is the start of the middle p values                        
   do i=1,p
    Bp(i,1:NRHS)=B(j+i-1,1:NRHS)
    do k=1,p
     if (KU+k-i+1 >= 1 .AND. KU+k-i+1 <= 2*KU+1 ) then
      Cp(i,k)=AB(KU+k-i+1,j+i-1)
     endif 
    end do
    do k=KU+1,2*KU-p+i 
     Sp(i,k)=AB(k+p-i+1,j+i-1)
    end do 
    do k=1,KU
     if (k >= i) then
      if (1+k-i >= 1 .AND. 1+k-i <= 2*KU+1 ) then
       Sp(i,k)=AB(1+k-i,j+i-1)
      endif
     endif
    end do
   end do       
!   solve for UE and precursor of UD
!   Cp zpj = -Sp zj + Bp
!   concatenate UD precursor and UD solutions onto EEK        
    EEK(:,1:2*KU)=-Sp(:,:)
    do hh=1,NRHS
      EEK(:,2*KU+hh)=Bp(:,hh)     
    end do     
    call DGESV(p,2*KU+NRHS,Cp,p,IPIV,EEK,p,INFO) ! overwrites EEK into solution, overwrites Cp into factors 
!    call GaussJordan( p,2*KU+NRHS,Cp,p,EEK,p,INFO )   ! overwrites EEK into solution, overwrites Cp into inverse     
    if (info /= 0) then
     CALL XERBLA( 'DGETRS/DCBSV ', INFO )
!     RETURN
    else
!    Two cases, p < KU or p >= KU (p <=2*KU always since p=mod(N,2*KU))    
!    p rows on top and bottom    
     if ( p <= KU ) then     
     do i=1,p    
      UDR(i,:,(N-p)/(2*KU)+1)=EEK(i,1:2*KU)
      UDR(2*KU-i+1,:,(N-p)/(2*KU)+1)=EEK(p-i+1,1:2*KU)                  
      do hh=1,NRHS     
       UER(i,(N-p)/(2*KU)+1,hh)=EEK(i,2*KU+hh)
       UER(2*KU-i+1,(N-p)/(2*KU)+1,hh)=EEK(p-i+1,2*KU+hh)                 
      end do
     end do
     else     ! p > KU
     do i=1,KU    
      UDR(i,:,(N-p)/(2*KU)+1)=EEK(i,1:2*KU)
      UDR(2*KU-i+1,:,(N-p)/(2*KU)+1)=EEK(p-i+1,1:2*KU)                  
      do hh=1,NRHS     
       UER(i,(N-p)/(2*KU)+1,hh)=EEK(i,2*KU+hh)
       UER(2*KU-i+1,(N-p)/(2*KU)+1,hh)=EEK(p-i+1,2*KU+hh)                 
      end do
     end do
     endif          
    endif  ! info =0                        
   endif ! p /=0

!$OMP PARALLEL num_threads(2) PRIVATE(thread)
!  separate iterative parts into subroutines so each thread can work in parallel
   thread = omp_get_thread_num()
   if (thread==0) then
    call forward_dcbsv( N, KU, Bj,Cj,Pj,Sj, NRHS, INFO, UD, UE, JJ)
   else
    call backward_dcbsv( N, KU, Bj,Cj,Pj,Sj, NRHS, INFO, UDR, UER, LL)
   endif
!$OMP END PARALLEL

   jj=LL  ! count forward is count backward + 1
                     
!  LAST EQUATIONS for both is at jj determined above when j=L    
!  B(j)=Inverse[IDENT-matmul(udR(:,:,j),ud(:,:,j-1))][udR(:,:,j)*ue(j-1)+ueR(j)]
!  [IDENT-matmul(udR(:,:,j),ud(:,:,j-1))]B(j)   = udR(:,:,j)*ue(j-1)+ueR(j)
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,udR(:,:,jj),2*KU,ud(:,:,jj-1),2*KU,0.0_wp,AA,2*KU)
    A=IDENT-AA
    do hh=1,NRHS
     call DGEMV('N',2*KU,2*KU,1.0_wp,udR(:,:,jj),2*KU,UE(:,jj-1,hh),1,0.0_wp,CC(:,hh),1)
     CCL(:,hh)=CC(:,hh)+ueR(:,jj,hh)
    end do
!    A=IDENT-matmul(udR(:,:,jj),ud(:,:,jj-1))
!    CCL(:,1:NRHS)=matmul(udR(:,:,jj),ue(:,jj-1,1:NRHS))+ueR(:,jj,1:NRHS)
    call DGESV(2*KU, NRHS , A, 2*KU, IPIV, CCL(:,1:NRHS), 2*KU, INFO ) ! overwrites CCL 
!    call GaussJordan( 2*KU, NRHS, A ,2*KU, CCL(:,1:NRHS), 2*KU, INFO )  ! overwrites CCL       
    if (info /= 0) then        
     CALL XERBLA( 'DGESV/DCBSV ', INFO )
!     RETURN
    else              
      Bj(:,jj,1:NRHS)=CCL(:,1:NRHS)   ! write Bj with RHS                                        
    endif

!  B(j-1)=Inverse[IDENT-matmul(ud(:,:,j-1),udR(:,:,j))][ud(:,:,j-1)*ueR(j)+ue(j-1)]
!  [IDENT-matmul(ud(:,:,j-1),udR(:,:,j))]B(j-1) = ud(:,:,j-1)*ueR(j)+ue(j-1)
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,ud(:,:,jj-1),2*KU,udR(:,:,jj),2*KU,0.0_wp,AA,2*KU)
    A=IDENT-AA
    do hh=1,NRHS
     call DGEMV('N',2*KU,2*KU,1.0_wp,ud(:,:,jj-1),2*KU,UER(:,jj,hh),1,0.0_wp,CC(:,hh),1)
     CCL(:,hh)=CC(:,hh)+ue(:,jj-1,hh)
    end do
!    A=IDENT-matmul(ud(:,:,jj-1),udR(:,:,jj))
!    CCL(:,1:NRHS)=matmul(ud(:,:,jj-1),ueR(:,jj,1:NRHS))+ue(:,jj-1,1:NRHS)
    call DGESV(2*KU, NRHS , A, 2*KU, IPIV, CCL(:,1:NRHS), 2*KU, INFO ) ! overwrites CCL 
!    call GaussJordan( 2*KU, NRHS, A ,2*KU, CCL(:,1:NRHS), 2*KU, INFO )  ! overwrites CCL       
    if (info /= 0) then        
     CALL XERBLA( 'DGESV/DCBSV ', INFO )
!     RETURN
    else              
      Bj(:,jj-1,1:NRHS)=CCL(:,1:NRHS)   ! write Bj-1 with RHS                                         
    endif
  
!   BACKSUBSTITUTION z(j-1)=UE(j-1)+UD(:,:,j-1)*z(j) 
    do jj=LL-1,2,-1
     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UD(:,:,jj-1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+CC(:,1:NRHS)
!      Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+matmul(UD(:,:,jj-1),Bj(:,jj,1:NRHS))
    end do
!   BACKSUBSTITUTION z(j+1)=UER(j+1)+UDR(:,:,j+1)*z(j) 
    do jj=LL,(N-p)/(2*KU)
     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UDR(:,:,jj+1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
     Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+CC(:,1:NRHS)
!      Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+matmul(UDR(:,:,jj+1),Bj(:,jj,1:NRHS))
    end do  

!   overwrite RHS with solution Bj, direction irrelevant    
    do kk=1,NRHS                       
     ii=(N-p)/(2*KU)+1          
     do j=(N-p)/2+1,1,-KU            
      do i=1,KU
       B(j+i-1,kk)=Bj(i,ii,kk) 
      end do   
      do i=KU+1,2*KU
       B(N-2*KU+i-j+1,kk)=Bj(i,ii,kk) 
      end do 
     ii=ii-1      
     end do        
    end do  ! end kk

   deallocate (Bj,Cj,Pj,Sj,UD,UE,UDR,UER)
   deallocate(A,AA,CC,EE,CCL,IDENT,IPIV)
   if (p /= 0) then
    deallocate (Cp,Sp,Bp,EEK)
   endif
        
  END SUBROUTINE dcbsv


  subroutine forward_dcbsv( N, KU, Bj,Cj,Pj,Sj, NRHS, INFO, UD, UE, LL)
  Use lapackinterface
  IMPLICIT NONE   
!  PURPOSE forward iterative loop
!  Copyright (c) 2021   Anthony M de Beus
   INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision

!  .. Scalar Arguments ..
   Integer, Intent(IN) ::  KU, N, NRHS
   INTEGER, INTENT(OUT) :: INFO
!  .. Array Arguments .. 

!  .. Work space ..
   REAL(wp), Intent(IN) :: Bj(2*KU,N/(2*KU)+1,NRHS)  
   REAL(wp), Intent(IN) :: Cj(2*KU,2*KU,N/(2*KU)+1)
   REAL(wp), Intent(IN) :: Pj(2*KU,2*KU,N/(2*KU)+1)
   REAL(wp), Intent(IN) :: Sj(2*KU,2*KU,N/(2*KU)+1)       
   REAL(wp) :: UD(2*KU,2*KU,0:N/(4*KU)+1)      ! ud is my set of matrices Aj
   REAL(wp) :: UE(2*KU,0:N/(4*KU)+1,NRHS)      ! ue is my vectors vj  
   REAL(wp) :: A(2*KU,2*KU),AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS) ! working copies
   INTEGER ::  i,j,k,kk,hh,p,ii,jj,ipiv(2*KU),L,LL

   p=mod(N,2*KU)
   L=(N-p)/4

!! FORWARD   
   !  ALL BUT THE LAST EQUATION, GENERATE UD & UE  ! (Cj+Sj*A(:,:,j-1))*A(:,:,j)=-Pj
   jj=0  !index of number of arrays
   do j=L,(N-p)/2,KU      
    jj=jj+1     
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,Sj(:,:,jj),2*KU,UD(:,:,jj-1),2*KU,0.0_wp,AA,2*KU)
    A=Cj(:,:,jj)+AA
!    A=Cj(:,:,jj)+matmul(Sj(:,:,jj),UD(:,:,jj-1))        
!   concatenate Aj and vj solutions onto EE        
     EE(:,1:2*KU)=-Pj(:,:,jj)
     do hh=1,NRHS
      call DGEMV('N',2*KU,2*KU,1.0_wp,Sj(:,:,jj),2*KU,UE(:,jj-1,hh),1,0.0_wp,CC(:,hh),1) 
      EE(:,2*KU+hh)=Bj(:,jj,hh)-CC(:,hh)     
!      EE(:,2*KU+hh)=Bj(:,jj,hh)-matmul(Sj(:,:,jj),UE(:,jj-1,hh))
     end do
!    compute next UD,UE using factored A
    call DGESV(2*KU,2*KU+NRHS,A,2*KU,IPIV,EE,2*KU,INFO) ! overwrites EE into solution 
!    call GaussJordan( 2*KU, 2*KU+NRHS ,A ,2*KU , EE, 2*KU, INFO )   ! overwrites EE into solution                   
    if (info /= 0) then
     CALL XERBLA( 'DGETRS/DCBSV ', INFO )
!     RETURN
    else
     UD(:,:,jj)=EE(:,1:2*KU)
    do hh=1,NRHS     
     UE(:,jj,hh)=EE(:,2*KU+hh)                
    end do    
    endif              
   end do
   LL=jj
  end subroutine forward_dcbsv

  subroutine backward_dcbsv(  N, KU, Bj,Cj,Pj,Sj, NRHS, INFO, UDR, UER, LL)
  Use lapackinterface
  IMPLICIT NONE   
!  PURPOSE backward iterative loop
!  Copyright (c) 2021   Anthony M de Beus
   INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision

!  .. Scalar Arguments ..
   Integer, Intent(IN) ::  KU, N, NRHS
   INTEGER, INTENT(OUT) :: INFO
!  .. Array Arguments ..

!  .. Work space ..
   REAL(wp), Intent(IN) :: Bj(2*KU,N/(2*KU)+1,NRHS)  
   REAL(wp), Intent(IN) :: Cj(2*KU,2*KU,N/(2*KU)+1)
   REAL(wp), Intent(IN) :: Pj(2*KU,2*KU,N/(2*KU)+1)
   REAL(wp), Intent(IN) :: Sj(2*KU,2*KU,N/(2*KU)+1)       
   REAL(wp) :: UDR(2*KU,2*KU,N/(4*KU):N/(2*KU)+1)     ! udr is my set of matrices Ajbar
   REAL(wp) :: UER(2*KU,N/(4*KU):N/(2*KU)+1,NRHS)     ! uer is my vectors vjbar 
   REAL(wp) :: A(2*KU,2*KU),AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS) ! working copies
   INTEGER ::  i,j,k,kk,hh,p,ii,jj,ipiv(2*KU),L,LL

   p=mod(N,2*KU)
   L=(N-p)/4

!! BACKWARD       
!  ALL BUT THE LAST EQUATION, GENERATE UDR & UER  ! (Cj+Pj*A(:,:,j+1))*A(:,:,j)=-Sj
   jj=(N-p)/(2*KU)+1  !index of number of arrays   
    do j=(N-p)/2,L+KU,-KU ! j not used; just a counter      
    jj=jj-1    
    call DGEMM('N','N',2*KU,2*KU,2*KU,1.0_wp,Pj(:,:,jj),2*KU,UDR(:,:,jj+1),2*KU,0.0_wp,AA,2*KU)
    A=Cj(:,:,jj)+AA
!    A=Cj(:,:,jj)+matmul(Pj(:,:,jj),UDR(:,:,jj+1))        
!   concatenate Aj and vj solutions onto EE        
     EE(:,1:2*KU)=-Sj(:,:,jj)
     do hh=1,NRHS
      call DGEMV('N',2*KU,2*KU,1.0_wp,Pj(:,:,jj),2*KU,UER(:,jj+1,hh),1,0.0_wp,CC(:,hh),1) 
      EE(:,2*KU+hh)=Bj(:,jj,hh)-CC(:,hh)     
!      EE(:,2*KU+hh)=Bj(:,jj,hh)-matmul(Pj(:,:,jj),UER(:,jj+1,hh))
     end do
!    compute next UD,UE using factored A
    call DGESV(2*KU,2*KU+NRHS,A,2*KU,IPIV,EE,2*KU,INFO) ! overwrites EE into solution 
!    call GaussJordan( 2*KU, 2*KU+NRHS ,A ,2*KU , EE, 2*KU, INFO )   ! overwrites EE into solution                    
    if (info /= 0) then
     CALL XERBLA( 'DGETRS/DCBSV ', INFO )
!     RETURN
    else
     UDR(:,:,jj)=EE(:,1:2*KU)
    do hh=1,NRHS     
     UER(:,jj,hh)=EE(:,2*KU+hh)                
    end do    
    endif  
   end do
   LL=jj
  end subroutine backward_dcbsv





   SUBROUTINE DCBSV( N, KU, NRHS, AB, LDAB, B, LDB, INFO )
   Use lapackinterface 
   USE OMP_LIB      
   IMPLICIT NONE
!   Copyright (c) 2021   Anthony M de Beus
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
   REAL(wp),ALLOCATABLE :: UDR(:,:,:)     ! udr is my set of matrices Ajbar
   REAL(wp),ALLOCATABLE :: UER(:,:,:)     ! uer is my vectors vjbar 
   REAL(wp),ALLOCATABLE :: ACL(:,:),CCL(:,:)   !last equation working copies
   REAL(wp),ALLOCATABLE :: CC(:,:),IDENT(:,:)   
   REAL(wp), ALLOCATABLE :: Cp(:,:),Sp(:,:),Bp(:,:),EEK(:,:) ! pxp and px2*KU, and sometimes p=0
   INTEGER, ALLOCATABLE :: ipiv(:)
   INTEGER ::  i,j,k,kk,hh,p,ii,jj,L1,LL,thread,allocstat
    
   p=mod(N,2*KU)
   L1=(N-p)/(4*KU)

   allocate(Bj(2*KU,(N-p)/(2*KU)+1,NRHS),Cj(2*KU,2*KU,(N-p)/(2*KU)+1))
   allocate(Pj(2*KU,2*KU,(N-p)/(2*KU)+1),Sj(2*KU,2*KU,(N-p)/(2*KU)+1))
   allocate(UD(2*KU,2*KU,0:L1+1),UE(2*KU,0:L1+1,NRHS))
   allocate(UDR(2*KU,2*KU,L1:(N-p)/(2*KU)+1),UER(2*KU,L1:(N-p)/(2*KU)+1,NRHS))
   allocate(IDENT(2*KU,2*KU),IPIV(4*KU),STAT=allocstat)   ! need IPIV(4*KU) not 2*KU for alternate last equation solver
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
   ipiv=0;   Bj=0;   Cj=0;   Pj=0;  Sj=0 
   UER=0;   UDR=0;   UE=0;   UD=0  
   jj=0  !index of number of arrays, replaced by loop invariant (j-1)/KU+1 below to enable parallelism

!   Setup the arrays  
!$OMP PARALLEL DO PRIVATE(i,j,k)     
    do j=1,(N-p)/2+1-KU,KU                     
         do i=1,KU
           Bj(i,(j-1)/KU+1,1:NRHS)=B(j+i-1,1:NRHS)
          do k=1,KU
           Cj(i,k,(j-1)/KU+1)=AB(KU+k-i+1,j+i-1)
          end do
          do k=1,i
            Pj(i,k,(j-1)/KU+1)=AB(2*KU+1+k-i,j+i-1)
          end do
          do k=i,KU
            Sj(i,k,(j-1)/KU+1)=AB(1+k-i,j+i-1)
          end do 
         end do                        
         do i=KU+1,2*KU
          Bj(i,(j-1)/KU+1,1:NRHS)=B(n-2*KU+i-j+1,1:NRHS)        
          do k=KU+1,2*KU
           Cj(i,k,(j-1)/KU+1)=AB(KU+k-i+1,n-2*KU+i-j+1)
          end do
          do k=i,2*KU 
            Pj(i,k,(j-1)/KU+1)=AB(1+k-i,n-2*KU+i-j+1)
          end do
          do k=KU+1,i 
            Sj(i,k,(j-1)/KU+1)=AB(2*KU+1+k-i,n-2*KU+i-j+1)
          end do 
         end do                                                              
   end do
!$OMP END PARALLEL DO  

!  IDENTITY MATRIX 
   IDENT=0 
   forall(j = 1:2*KU) IDENT(j,j) = 1

!  FIRST ARRAY forward uses p=0 formula
!  already done above UE(:,0,:)=0;  UD(:,:,0)=0
   do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UD(i,k,0)=1        
      endif
     end do 
   end do        
               
!  FIRST EQUATION in reverse if p=0 V((N-p)/(2*KU)=0 & A((N-p)/(2*KU) = permuted identity
   if (p == 0) then
!  already done above UER(:,(N-p)/(2*KU)+1,:)=0;    UDR(:,:,(N-p)/(2*KU)+1)=0
    do i=1,2*KU
     do k=1,2*KU
      if ( (ABS(k - i) - KU) == 0 ) then
        UDR(i,k,(N-p)/(2*KU)+1)=1
      endif
     end do 
    end do
   else
   allocate (Cp(p,p),Sp(p,2*KU),Bp(p,NRHS),EEK(p,2*KU+NRHS),STAT=allocstat)
   if (allocstat /=0) then
    write(*,*) 'Memory allocation failed'
    stop
   endif       
!  if p /= 0 have to compute UD(:,:,(N-p)/(2*KU)+1), UE(:,(N-p)/(2*KU)+1,:)
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
    deallocate (Cp,Sp,Bp,EEK)              
   endif ! p /=0

!$OMP PARALLEL num_threads(2) PRIVATE(thread)
!  separate iterative parts into subroutines so each thread can work in parallel
   thread = omp_get_thread_num()
   if (thread==0) then
    call forward_loop(L1, N ,KU, size(Bj,2),Bj,Cj,Pj,Sj, NRHS, INFO, L1+1, UD, UE, JJ) !size(UD,3)-1 == L1+1 =(N-p)/(4*KU)+1     
write(*,*) 'thread,jj',thread,JJ
   else
    call backward_loop(L1, N ,KU, size(Bj,2),Bj,Cj,Pj,Sj, NRHS, INFO, L1, size(Bj,2),UDR, UER, LL) !dimensions UDR L1,size(Bj,2)=(N-p)/(2*KU)+1   
write(*,*) 'thread,ll',thread,LL
   endif
!$OMP END PARALLEL


write(*,*) 'L1,jj,ll',L1,jj,ll
jj=ll
write(*,*) ud(:,:,JJ-1)
write(*,*) ' '
write(*,*) udr(:,:,JJ)
write(*,*) ' '

   deallocate(Cj,Pj,Sj)
!  LAST EQUATIONS for both is at jj determined above when j=L  
   jj=LL  ! count forward is count backward + 1
   
   allocate(CC(2*KU,NRHS))

   allocate(ACL(4*KU,4*KU),CCL(4*KU,NRHS))
   ACL=0; CCL=0
   ACL(1:2*KU,1:2*KU)=udR(:,:,JJ)
   ACL(1:2*KU,2*KU+1:4*KU)=-IDENT
   ACL(2*KU+1:4*KU,1:2*KU)=IDENT
   ACL(2*KU+1:4*KU,2*KU+1:4*KU)=-ud(:,:,JJ-1)

   CCL(1:2*KU,1:NRHS)=-uER(:,JJ,1:NRHS)        !ue, ueR are vectors v 
   CCL(2*KU+1:4*KU,1:NRHS)=uE(:,JJ-1,1:NRHS)

    call DGESV(4*KU, NRHS , ACL, 4*KU, IPIV, CCL(:,1:NRHS), 4*KU, INFO ) ! overwrites CCL 
!    call GaussJordan( 4*KU, NRHS, ACL ,4*KU, CCL(:,1:NRHS), 4*KU, INFO )  ! overwrites CCL

    if (info /= 0) then        
     CALL XERBLA( 'DGESV/DCBSV ', INFO )
!     RETURN
    else              
      Bj(:,JJ-1,1:NRHS)=CCL(1:2*KU,1:NRHS)   ! write Bj-1 with RHS z jj-1
      Bj(:,JJ,1:NRHS)=CCL(2*KU+1:4*KU,1:NRHS)   ! write Bj with RHS z hat j                                     
    endif     
    
   deallocate(ACL,CCL)
  
!   BACKSUBSTITUTION z(j-1)=UE(j-1)+UD(:,:,j-1)*z(j) 
    do jj=LL-1,2,-1
     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UD(:,:,jj-1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+CC(:,1:NRHS)
!     Bj(:,jj-1,1:NRHS)=UE(:,jj-1,1:NRHS)+matmul(UD(:,:,jj-1),Bj(:,jj,1:NRHS))
    end do
!   BACKSUBSTITUTION z(j+1)=UER(j+1)+UDR(:,:,j+1)*z(j) 
    do jj=LL,(N-p)/(2*KU)
     call DGEMM('N','N',2*KU,NRHS,2*KU,1.0_wp,UDR(:,:,jj+1),2*KU,Bj(:,jj,1:NRHS),2*KU,0.0_wp,CC(:,1:NRHS),2*KU)    
     Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+CC(:,1:NRHS)
!     Bj(:,jj+1,1:NRHS)=UER(:,jj+1,1:NRHS)+matmul(UDR(:,:,jj+1),Bj(:,jj,1:NRHS))
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

   deallocate (Bj,UD,UE,UDR,UER)
   deallocate(CC,IDENT,IPIV)
        
  END SUBROUTINE dcbsv
  








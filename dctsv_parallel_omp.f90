     SUBROUTINE DCTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
     USE OMP_LIB     
     IMPLICIT NONE   
!      PURPOSE solves the cyclic/periodic tridiagonal system, see LAPACK routine DGTSV for comparison
!      Copyright (c) 2021   Anthony M de Beus
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!      Arguments copied/modified from dgtsv.f *  -- LAPACK routine (version 3.1) --
!      Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!      November 2006
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 4.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input) working precision array, dimension (N)
!          On entry, DL must contain the sub-diagonal elements of A. The first term
!          wraps around to the last column first row
!
!  D       (input) working precision array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!
!  DU      (input) working precision array, dimension (N)
!          On entry, DU must contain the super-diagonal elements
!          of A. The last term wraps around to the first column last row
!
!  B       (input/output) working precision array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N by NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, a determinant is exactly zero, and the solution
!               has not been computed.  
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 9
!              D1  DU1 0   0   0   0   0   0   DL1
!              DL2 D2  DU2 0   0   0   0   0   0
!              0   DL3 D3  DU3 0   0   0   0   0
!              0   0   DL4 D4  DU4 0   0   0   0
!              0   0   0   DL5 D5  DU5 0   0   0
!              0   0   0   0   DL6 D6  DU6 0   0
!              0   0   0   0   0   DL7 D7  DU7 0
!              0   0   0   0   0   0   DL8 D8  DU8
!              DU9 0   0   0   0   0   0   DL9 D9
!
!      .. Scalar Arguments ..
       INTEGER, INTENT(IN) :: LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO 
!      ..
!      .. Array Arguments ..
       REAL(wp), INTENT(IN) :: D( N ), DL( N ), DU( N )  ! no output no LU factors
       REAL(wp), INTENT(INOUT) :: B( LDB, NRHS )         ! on entry RHS, on exit, solution 

!      make a copy to avoid data blocking 
       REAL(wp) :: DR( N ), DLR( N ), DUR( N )  ! no output no LU factors
       REAL(wp) :: BR( N, NRHS )                ! on entry RHS, on exit, solution        
!      ..
       REAL(wp) :: udR(2,2,N/4:N/2+1),ueR(2,N/4:N/2+1,NRHS),DET2  ! udR is my set of matrices Aj(bar) ueR is my vectors vj(bar)
       REAL(wp) :: ud(2,2,0:N/4+1),ue(2,0:N/4+1,NRHS),DET         !  ud is my set of matrices Aj ue is my vectors vj

       INTEGER :: i,j,k,p,L,thread

!      needed for f90+ calling of f77 routines
       INTERFACE
          SUBROUTINE XERBLA( SRNAME, INFO )
!         .. Scalar Arguments ..
          CHARACTER*6        SRNAME          !CHARACTER(6)
          INTEGER            INFO
          END SUBROUTINE XERBLA      
       END INTERFACE

!      INFO handling copied/modified from dgtsv.f *  -- LAPACK routine (version 3.1) --
!      Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!      November 2006
       INFO = 0
       IF( N.LT.4 ) THEN
         INFO = -1
       ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
       ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
       END IF
       IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DCTSV ', -INFO )
!         RETURN
       END IF
!       IF( N.EQ.0 ) RETURN
!      End INFO handling

       DR=D
       DLR=DL
       DUR=DU
       BR=B
       
       p=mod(N,2)
       L=(N-p)/4
      
!      FIRST EQUATION ud(0)=((0,1),(1,0) & ue(0) = (0,0)
       ud(:,:,0)=0
       ud(1,2,0)=1
       ud(2,1,0)=1 
       ue(:,0,:)=0 

!      FIRST EQUATION 
!       udR(0)=((0,1),(1,0) & ueR(0) = (0,0) when p==0
       if (p == 0) then 
        udR(:,:,N/2+1)=0
        udR(1,2,N/2+1)=1
        udR(2,1,N/2+1)=1 
        ueR(:,N/2+1,:)=0  
       else
        if ( D((N+1)/2) /= 0) then
        udR(1,1,(N+1)/2)=-DL((N+1)/2)/D((N+1)/2)
        udR(1,2,(N+1)/2)=-DU((N+1)/2)/D((N+1)/2)
        udR(2,1,(N+1)/2)=-DL((N+1)/2)/D((N+1)/2)
        udR(2,2,(N+1)/2)=-DU((N+1)/2)/D((N+1)/2)
        ueR(1,(N+1)/2,1:NRHS)=B((N+1)/2,1:NRHS)/D((N+1)/2)
        ueR(2,(N+1)/2,1:NRHS)=B((N+1)/2,1:NRHS)/D((N+1)/2) 
        else
         INFO=(N+1)/2
         CALL XERBLA( 'DCTSV ', -INFO )
!         RETURN
        endif                           
       endif        

!$OMP PARALLEL num_threads(2) PRIVATE(thread)
!      separate iterative parts into subroutines so each thread can work in parallel
       thread = omp_get_thread_num()
       if (thread==0) then
        call forward_dctsv(N, NRHS, DL, D, DU, B, LDB, INFO, UD, UE)
       else
        call backward_dctsv(N, NRHS, DLR, DR, DUR, BR, LDB, INFO, UDR, UER)
       endif
!$OMP END PARALLEL

!      LAST EQUATIONS for both is at L
       j=L    
!      B(j)=Inverse[IDENT-matmul(udR(:,:,j),ud(:,:,j-1))][udR(:,:,j)*ue(j-1)+ueR(j)]
              
       DET=(1-udR(1,1,j)*ud(1,1,j-1)-udR(1,2,j)*ud(2,1,j-1))*&
            (1-udR(2,1,j)*ud(1,2,j-1)-udR(2,2,j)*ud(2,2,j-1))-&
              (udR(1,1,j)*ud(1,2,j-1)+udR(1,2,j)*ud(2,2,j-1))*&
              (udR(2,1,j)*ud(1,1,j-1)+udR(2,2,j)*ud(2,1,j-1))
                            

       DET2=(1-ud(1,1,j-1)*udR(1,1,j)-ud(1,2,j-1)*udR(2,1,j))*&
            (1-ud(2,1,j-1)*udR(1,2,j)-ud(2,2,j-1)*udR(2,2,j))-&
              (ud(1,1,j-1)*udR(1,2,j)+ud(1,2,j-1)*udR(2,2,j))*&
              (ud(2,1,j-1)*udR(1,1,j)+ud(2,2,j-1)*udR(2,1,j))


        IF (DET /= 0 .AND. DET2 /= 0) THEN           

         B(j,1:NRHS) = ((1-udR(2,1,j)*ud(1,2,j-1)-udR(2,2,j)*ud(2,2,j-1))*&
               (udR(1,1,j)*ue(1,j-1,1:NRHS)+udR(1,2,j)*ue(2,j-1,1:NRHS)+ ueR(1,j,1:NRHS))+&
               (udR(1,1,j)*ud(1,2,j-1)+udR(1,2,j)*ud(2,2,j-1))*(udR(2,1,j)*ue(1,j-1,1:NRHS)&
               +udR(2,2,j)*ue(2,j-1,1:NRHS)+ ueR(2,j,1:NRHS)))/DET                       
         
         B(N-j+1,1:NRHS) = ((1-udR(1,1,j)*ud(1,1,j-1)-udR(1,2,j)*ud(2,1,j-1))*&
               (udR(2,1,j)*ue(1,j-1,1:NRHS)+udR(2,2,j)*ue(2,j-1,1:NRHS)+ ueR(2,j,1:NRHS))+&
               (udR(2,1,j)*ud(1,1,j-1)+udR(2,2,j)*ud(2,1,j-1))*(udR(1,1,j)*ue(1,j-1,1:NRHS)&
               +udR(1,2,j)*ue(2,j-1,1:NRHS)+ ueR(1,j,1:NRHS)))/DET


        if (L > 1 ) then 
        B(j-1,1:NRHS) = ((1-ud(2,1,j-1)*udR(1,2,j)-ud(2,2,j-1)*udR(2,2,j))*&
               (ud(1,1,j-1)*ueR(1,j,1:NRHS)+ud(1,2,j-1)*ueR(2,j,1:NRHS)+ ue(1,j-1,1:NRHS))+&
               (ud(1,1,j-1)*udR(1,2,j)+ud(1,2,j-1)*udR(2,2,j))*(ud(2,1,j-1)*ueR(1,j,1:NRHS)&
               +ud(2,2,j-1)*ueR(2,j,1:NRHS)+ ue(2,j-1,1:NRHS)))/DET2        
         
         B(N-j+2,1:NRHS) = ((1-ud(1,1,j-1)*udR(1,1,j)-ud(1,2,j-1)*udR(2,1,j))*&
               (ud(2,1,j-1)*ueR(1,j,1:NRHS)+ud(2,2,j-1)*ueR(2,j,1:NRHS)+ ue(2,j-1,1:NRHS))+&
               (ud(2,1,j-1)*udR(1,1,j)+ud(2,2,j-1)*udR(2,1,j))*(ud(1,1,j-1)*ueR(1,j,1:NRHS)&
               +ud(1,2,j-1)*ueR(2,j,1:NRHS)+ ue(1,j-1,1:NRHS)))/DET2
         endif
                                                              
        ELSE

         INFO=j
         CALL XERBLA( 'DCTSV ', -INFO )
!         RETURN
        ENDIF       

!     These are iterative and should have a problem with OMP, but apparently not with num_threads=2.
               
!      BACKSUBSTITUTION B(j-1)=UE(j-1)+UD(:,:,j-1)*B(j)
       do i=L-1,1,-1
        B(i,1:NRHS)=    ue(1,i,1:NRHS)+ud(1,1,i)*B(i+1,1:NRHS)+ud(1,2,i)*B(N-i,1:NRHS)
        B(N-i+1,1:NRHS)=ue(2,i,1:NRHS)+ud(2,1,i)*B(i+1,1:NRHS)+ud(2,2,i)*B(N-i,1:NRHS) 
       end do

!      BACKSUBSTITUTION B(j+1)=UER(j+1)+UDR(:,:,j+1)*B(j)
       do i=L,(N-p)/2-1
        B(i+1,1:NRHS)=ueR(1,i+1,1:NRHS)+udR(1,1,i+1)*B(i,1:NRHS)+udR(1,2,i+1)*B(N-i+1,1:NRHS)
        B(N-i,1:NRHS)=ueR(2,i+1,1:NRHS)+udR(2,1,i+1)*B(i,1:NRHS)+udR(2,2,i+1)*B(N-i+1,1:NRHS) 
       end do

       if (p /=0) then  ! take the average for the middle value
        i=(N-p)/2
        B(i+1,1:NRHS)=(ueR(1,i+1,1:NRHS)+udR(1,1,i+1)*B(i,1:NRHS)+udR(1,2,i+1)*B(N-i+1,1:NRHS)+&
                       ueR(2,i+1,1:NRHS)+udR(2,1,i+1)*B(i,1:NRHS)+udR(2,2,i+1)*B(N-i+1,1:NRHS))/2 
       endif
      
     end subroutine DCTSV
     
!     The two iterative loop routines follow     

       subroutine forward_dctsv(N, NRHS, DL, D, DU, B, LDB, INFO, UD, UE)
       IMPLICIT NONE   
!      PURPOSE forward iterative loop
!      Copyright (c) 2021   Anthony M de Beus
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
       INTEGER, INTENT(IN) :: LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO 
!      ..
!      .. Array Arguments ..
       REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * ), B( LDB, * )  
       REAL(wp) :: ud(2,2,0:N/4+1),ue(2,0:N/4+1,NRHS)
!      ..
       REAL(wp) :: DET  !  ud is my set of matrices Aj ue is my vectors vj 
       INTEGER :: i,j,p,L

       p=mod(N,2)
       L=(N-p)/4

!      ALL BUT THE LAST EQUATION
       do j=1,L                                                                                   
        DET=D(j)*D(1-j+N)+DL(j)*D(1-j+N)*ud(1,1,-1+j)-& 
            DL(j)*DU(1-j+N)*ud(1,2,-1+j)*ud(2,1,-1+j)+&       
            D(j)*DU(1-j+N)*ud(2,2,-1+j)+&
            DL(j)*DU(1-j+N)*ud(1,1,-1+j)*ud(2,2,-1+j)
        IF (DET /= 0) THEN
         ud(1,1,j)=-((DU(j)*(D(1-j+N)+DU(1-j+N)*ud(2,2,-1+j)))/DET)    
         ud(1,2,j)=(DL(j)*DL(1-j+N)*ud(1,2,-1+j))/DET                 
         ud(2,1,j)=(DU(j)*DU(1-j+N)*ud(2,1,-1+j))/DET                
         ud(2,2,j)=-((DL(1-j+N)*(D(j)+DL(j)*ud(1,1,-1+j)))/DET)
         ue(1,j,1:NRHS)=(-DL(j)*(B(1-j+N,1:NRHS)-DU(1-j+N)*ue(2,-1+j,1:NRHS))*ud(1,2,-1+j)+&             
                (B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(D(1-j+N)+DU(1-j+N)*ud(2,2,-1+j)))/DET
         ue(2,j,1:NRHS)=((B(1-j+N,1:NRHS)-DU(1-j+N)*ue(2,-1+j,1:NRHS))*(D(j)+DL(j)*ud(1,1,-1+j))-&
                DU(1-j+N)*(B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*ud(2,1,-1+j))/DET          
        ELSE
         INFO=j
         CALL XERBLA( 'DCTSV ', INFO )
!         RETURN
        ENDIF                                                                                              
       end do
       end subroutine forward_dctsv

       subroutine backward_dctsv(N, NRHS, DL, D, DU, B, LDB, INFO, UDR, UER)
       IMPLICIT NONE   
!      PURPOSE backward iterative loop
!      Copyright (c) 2021   Anthony M de Beus
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
       INTEGER, INTENT(IN) :: LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO 
!      ..
!      .. Array Arguments ..
       REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * ), B( LDB, * ) 
       REAL(wp) :: udR(2,2,N/4:N/2+1),ueR(2,N/4:N/2+1,NRHS) 
!      ..
       REAL(wp) :: DET  !  ud is my set of matrices Aj ue is my vectors vj 
       INTEGER :: i,j,p,L

       p=mod(N,2)
       L=(N-p)/4

!      ALL BUT THE LAST EQUATION  
       do j=(N-p)/2,L,-1                                              
        DET=D(j)*D(1-j+N)+DU(j)*D(1-j+N)*udR(1,1,1+j)-& 
            DU(j)*DL(1-j+N)*udR(1,2,1+j)*udR(2,1,1+j)+&       
            D(j)*DL(1-j+N)*udR(2,2,1+j)+&
            DU(j)*DL(1-j+N)*udR(1,1,1+j)*udR(2,2,1+j)
                                  
        IF (DET /= 0) THEN
         udR(1,1,j)=-((DL(j)*(D(1-j+N)+DL(1-j+N)*udR(2,2,1+j)))/DET)    
         udR(1,2,j)=(DU(j)*DU(1-j+N)*udR(1,2,1+j))/DET                 
         udR(2,1,j)=(DL(j)*DL(1-j+N)*udR(2,1,1+j))/DET                
         udR(2,2,j)=-((DU(1-j+N)*(D(j)+DU(j)*udR(1,1,1+j)))/DET)                                 
         ueR(1,j,1:NRHS)=(-DU(j)*(B(1-j+N,1:NRHS)-DL(1-j+N)*ueR(2,1+j,1:NRHS))*udR(1,2,1+j)+&             
                (B(j,1:NRHS)-DU(j)*ueR(1,1+j,1:NRHS))*(D(1-j+N)+DL(1-j+N)*udR(2,2,1+j)))/DET
         ueR(2,j,1:NRHS)=((B(1-j+N,1:NRHS)-DL(1-j+N)*ueR(2,1+j,1:NRHS))*(D(j)+DU(j)*udR(1,1,1+j))-&
                DL(1-j+N)*(B(j,1:NRHS)-DU(j)*ueR(1,1+j,1:NRHS))*udR(2,1,1+j))/DET                                                  
        ELSE
         INFO=j
         CALL XERBLA( 'DCTSV ', -INFO )
!         RETURN
        ENDIF                                                                                                           
       end do
       end subroutine backward_dctsv
              


       

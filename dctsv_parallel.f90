     SUBROUTINE DCTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
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
       REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * )  ! no output no LU factors
       REAL(wp), INTENT(INOUT) :: B( LDB, * ) ! on entry RHS, on exit, solution 
!      ..

       REAL(wp) :: udR(2,2,0:N/2+1),ueR(2,0:N/2+1,NRHS),DETR
       REAL(wp) :: ud(2,2,0:N/2+1),ue(2,0:N/2+1,NRHS),DET,X(NRHS),Y(NRHS)  !  ud is my set of matrices Aj ue is my vectors vj 
       INTEGER :: i,j,k,p  

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
         RETURN
       END IF
       IF( N.EQ.0 ) RETURN
!      End INFO handling
       p=mod(N,2)
             
!      FIRST EQUATION ud(0)=((0,1),(1,0)
       ud(:,:,0)=0
       ud(1,2,0)=1
       ud(2,1,0)=1 
       ue(:,0,:)=0  
       udR(:,:,(N-p)/2)=0
       udR(1,2,(N-p)/2)=1
       udR(2,1,(N-p)/2)=1 
       ueR(:,(N-p)/2,:)=0  
              
!     FIRST EQUATION ue(0) = (0,0)

!      ALL BUT THE LAST EQUATION  
       do j=1,(N-p)/2-1                                                   
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
         RETURN
        ENDIF 
         
        k=(N-p)/2-j+1
        DETR=udR(1,1,k)*udR(2,2,k)-udR(2,1,k)*udR(1,2,k)
               
        IF (DET /= 0) THEN
         udR(1,1,k-1)=-D(k)/DL(k)-(DU(k)/DL(k))*(udR(2,2,k)/DETR)  
         udR(1,2,k-1)=(DU(k)/DL(k))*(udR(2,1,k)/DETR)  
         udR(2,1,k-1)=(DL(1-k+N)/DU(1-k+N))*(udR(1,2,k)/DETR)                                 
         udR(2,2,k-1)=-D(k)/DU(1-k+N)-(DL(1-k+N)/DU(1-k+N))*(udR(1,1,k)/DETR)
                  
         ueR(1,k-1,1:NRHS)=B(k,1:NRHS)/DL(k)-(ueR(1,k,1:NRHS)*udR(2,2,k)-&
                           ueR(2,k,1:NRHS)*udR(2,1,k))*DU(j)/(DL(k)*DETR)

         ueR(2,k-1,1:NRHS)=B(1-k+N,1:NRHS)/DU(1-k+N)-(ueR(1,k,1:NRHS)*udR(1,2,k)-&
                           ueR(2,k,1:NRHS)*udR(1,1,k))*DL(1-k+N)/(DU(1-k+N)*DETR)
                                                                                                                            
        ELSE
         INFO=j
         CALL XERBLA( 'DCTSV ', -INFO )
         RETURN
        ENDIF
        
        write(*,*) 'j:',j        
        write(*,*) Transpose(ud(:,:,j))
 !       write(*,*) ue(:,j,:)
        write(*,*) ' '
        write(*,*) 'k,j(k):',k-1, (N-p)/2-k+1      
        write(*,*) Transpose(udR(:,:,k-1))
 !       write(*,*) ueR(:,k-1,:)
        write(*,*) ' '        
        
                                                                                            
       end do 

!      LAST EQUATION 
       j=(N-p)/2
       if ( p == 0 ) then
!      EVEN CASE (2x2) 
!      NO ud(,,j); ue(,j) iS THE SOLUTION AT j,j+1  
        DET=D(j)*D(1+j)-DL(1+j)*DU(j)+DL(j)*D(1+j)*ud(1,1,-1+j)-&
            DL(j)*DL(1+j)*ud(1,2,-1+j)-DU(j)*DU(1+j)*ud(2,1,-1+j)-&
            DL(j)*DU(1+j)*ud(1,2,-1+j)*ud(2,1,-1+j)+D(j)*DU(1+j)*ud(2,2,-1+j)+&
            DL(j)*DU(1+j)*ud(1,1,-1+j)*ud(2,2,-1+j) 


        IF (DET /= 0) THEN           
         ue(1,j,1:NRHS)=((B(1+j,1:NRHS)-DU(1+j)*ue(2,-1+j,1:NRHS))*(-DU(j)-DL(j)*ud(1,2,-1+j))+&
           (B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(D(1+j)+DU(1+j)*ud(2,2,-1+j)))/DET                 
         ue(2,j,1:NRHS)=((B(1+j,1:NRHS)-DU(1+j)*ue(2,-1+j,1:NRHS))*(D(j)+DL(j)*ud(1,1,-1+j))+&
           (B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(-DL(1+j)-DU(1+j)*ud(2,1,-1+j)))/DET



         B(j,1:NRHS)=ue(1,j,1:NRHS)         
         B(j+1,1:NRHS)=ue(2,j,1:NRHS)                                                                        
        ELSE
         INFO=j
         CALL XERBLA( 'DCTSV ', -INFO )
         RETURN
        ENDIF       
        else
!       ODD CASE (3x3)
!       NO ud(,,j); B(j ETC.) iS THE SOLUTION AT j,j+1,j+2
        DET=D(j)*D(1+j)*D(2+j)-DL(1+j)*D(2+j)*DU(j)-DL(2+j)*D(j)*DU(1+j)+&           
            DL(j)*D(1+j)*D(2+j)*ud(1,1,-1+j)-DL(j)*DL(2+j)*DU(1+j)*ud(1,1,-1+j)+&
            DL(j)*DL(1+j)*DL(2+j)*ud(1,2,-1+j)+DU(j)*DU(1+j)*DU(2+j)*ud(2,1,-1+j)-&
            DL(j)*D(1+j)*DU(2+j)*ud(1,2,-1+j)*ud(2,1,-1+j)+D(j)*D(1+j)*DU(2+j)*ud(2,2,-1+j)-&
            DL(1+j)*DU(j)*DU(2+j)*ud(2,2,-1+j)+DL(j)*D(1+j)*DU(2+j)*ud(1,1,-1+j)*ud(2,2,-1+j)
        IF (DET /= 0) THEN  
         X(1:NRHS)=((B(2+j,1:NRHS)-DU(2+j)*ue(2,-1+j,1:NRHS))*(DU(j)*DU(1+j)-DL(j)*D(1+j)*ud(1,2,-1+j))+&                  
                (B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(D(1+j)*D(2+j)-DL(2+j)*DU(1+j)+D(1+j)*DU(2+j)*ud(2,2,-1+j))+&
                 B(1+j,1:NRHS)*(-D(2+j)*DU(j)+DL(j)*DL(2+j)*ud(1,2,-1+j))+B(1+j,1:NRHS)*(-DU(j)*DU(2+j)*ud(2,2,-1+j)))/DET 
         Y(1:NRHS)=((B(2+j,1:NRHS)-DU(2+j)*ue(2,-1+j,1:NRHS))*(D(j)*D(1+j)-DL(1+j)*DU(j)+DL(j)*D(1+j)*ud(1,1,-1+j))+& 
                (B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(DL(1+j)*DL(2+j)-D(1+j)*DU(2+j)*ud(2,1,-1+j))+&
                 B(1+j,1:NRHS)*(-DL(2+j)*D(j)-DL(j)*DL(2+j)*ud(1,1,-1+j)+DU(j)*DU(2+j)*ud(2,1,-1+j)))/DET
         B(j+1,1:NRHS)=((B(2+j,1:NRHS)-DU(2+j)*ue(2,-1+j,1:NRHS))*(-D(j)*DU(1+j)-DL(j)*DU(1+j)*ud(1,1,-1+j)+&
                 DL(j)*DL(1+j)*ud(1,2,-1+j))+(B(j,1:NRHS)-DL(j)*ue(1,-1+j,1:NRHS))*(-DL(1+j)*D(2+j)+&
                 DU(1+j)*DU(2+j)*ud(2,1,-1+j)-DL(1+j)*DU(2+j)*ud(2,2,-1+j))+&
                 B(1+j,1:NRHS)*(D(j)*D(2+j)+DL(j)*D(2+j)*ud(1,1,-1+j)-DL(j)*DU(2+j)*ud(1,2,-1+j)*&
                 ud(2,1,-1+j)+D(j)*DU(2+j)*ud(2,2,-1+j)+DL(j)*DU(2+j)*ud(1,1,-1+j)*ud(2,2,-1+j)))/DET                
         B(j,1:NRHS)=X(1:NRHS)
         B(j+2,1:NRHS)=Y(1:NRHS)                                                                       
        ELSE
         INFO=j
         CALL XERBLA( 'DCTSV ', -INFO )
         RETURN
        ENDIF                        
        endif               
!      BACKSUBSTITUTION B(j-1)=UE(j-1)+UD(:,:,j-1)*B(j)
       do i=(N-p)/2-1,1,-1
        B(i,1:NRHS)=    ue(1,i,1:NRHS)+ud(1,1,i)*B(i+1,1:NRHS)+ud(1,2,i)*B(N-i,1:NRHS)
        B(N-i+1,1:NRHS)=ue(2,i,1:NRHS)+ud(2,1,i)*B(i+1,1:NRHS)+ud(2,2,i)*B(N-i,1:NRHS)                
       end do
      
     end subroutine DCTSV
       

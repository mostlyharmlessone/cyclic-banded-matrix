       subroutine forward(N, NRHS, DL, D, DU, B, LDB, INFO, UD, UE)
       IMPLICIT NONE   
!      PURPOSE 
!      Copyright (c) 2021   Anthony M de Beus
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
       INTEGER, INTENT(IN) :: LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO 
!      ..
!      .. Array Arguments ..
       REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * ), B( LDB, * ) 
       REAL(wp) ::  ud(2,2,0:N/2+1),ue(2,0:N/2+1,NRHS) 
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
       end subroutine forward

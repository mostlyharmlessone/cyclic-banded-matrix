       subroutine backward(N, NRHS, DL, D, DU, B, LDB, INFO, UDR, UER)
       IMPLICIT NONE   
!      PURPOSE 
!      Copyright (c) 2021   Anthony M de Beus
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
       INTEGER, INTENT(IN) :: LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO 
!      ..
!      .. Array Arguments ..
       REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * ), B( LDB, * ) 
       REAL(wp) ::  udR(2,2,0:N/2+1),ueR(2,0:N/2+1,NRHS) 
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
       end subroutine backward

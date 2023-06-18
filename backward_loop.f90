  subroutine backward_loop( L, N, KU, LB,Bj,Cj,Pj,Sj, NRHS, INFO, LR,LU,UDR, UER, LL)
  Use lapackinterface
  IMPLICIT NONE   
!  PURPOSE backward iterative loop
!  Copyright (c) 2021   Anthony M de Beus
   INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision

!  .. Scalar Arguments ..
   Integer, Intent(IN) ::  L, LB, LR, LU, KU, N, NRHS  ! L is starting place, LB=size(B,2)=size(C/P/J,3) 
                                                       ! LU+1=size(UD,3)=size(UE,2) index LR arrays ie. LR=0 for dcbsv_reverse
   INTEGER, INTENT(OUT) :: INFO,LL                     ! LL is number of steps
!  .. Array Arguments ..

   REAL(wp),INTENT(IN) :: Bj(2*KU,LB,NRHS)  
   REAL(wp),INTENT(IN) :: Cj(2*KU,2*KU,LB)
   REAL(wp),INTENT(IN) :: Pj(2*KU,2*KU,LB)
   REAL(wp),INTENT(IN) :: Sj(2*KU,2*KU,LB)    
   REAL(wp),INTENT(INOUT) :: UDR(2*KU,2*KU,LR:LU) ! ud is my set of matrices Aj                 
   REAL(wp),INTENT(INOUT) :: UER(2*KU,LR:LU,NRHS)      ! ue is my vectors vj
   
!  .. Work space ..   
   REAL(wp) :: A(2*KU,2*KU),AA(2*KU,2*KU),CC(2*KU,NRHS),EE(2*KU,2*KU+NRHS) ! working copies
   INTEGER ::  hh,p,jj,ipiv(2*KU)

   p=mod(N,2*KU)

!! BACKWARD       
!  ALL BUT THE LAST EQUATION, GENERATE UDR & UER  ! (Cj+Pj*A(:,:,j+1))*A(:,:,j)=-Sj
    do jj=(N-p)/(2*KU),L+1,-1     
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
   LL=jj+1
   
  end subroutine backward_loop
  
  

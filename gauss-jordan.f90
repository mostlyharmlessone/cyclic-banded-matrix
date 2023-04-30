  SUBROUTINE GaussJordan( N, NRHS, A, LDA, B, LDB, INFO )
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
! Copyright (c) 2021   Anthony M de Beus              
!          Arguments copied and modified from -- LAPACK routine (version 3.1) --
!          Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!          November 2006
!
!  A       (input/output) array, dimension (LDA,N)
!          On entry, the N-by-N coefficient matrix A.
!          On exit, the N-by-N inverse of A
!          
!  LDA    (input) INTEGER
!          The leading dimension of the array AB. LDA >= max(1,N) 
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N > 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  B       (input/output) array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N). 
!
!  IPIV    IPIV is INTEGER array, dimension (N)
!          The pivot indices that define the permutation matrix P;
!          row i of the matrix was interchanged with row IPIV(i).
!
!  INFO    INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero so the solution could not be computed. 

!     .. Scalar Arguments ..
  INTEGER,INTENT(IN)   :: LDA, LDB, N, NRHS
  INTEGER, INTENT(OUT) :: INFO 
!     ..
!     .. Array Arguments ..
  REAL(wp),INTENT(INOUT) ::  A( LDA, * ), B( LDB, * )
  
!  .. Work space ..  
  INTEGER ::  ipiv(N),icol(N),irow(N),i,j,k,ii,jj,index_row,index_col
  REAL(wp) :: largest,swap(N),temp,pivot,swap2(NRHS)
     
  info = 0
    IF( N.LT.0 ) THEN
         info = -1
    ELSE IF( NRHS.LT.0 ) THEN
         info = -2
    ELSE IF( LDA.LT.max( 1, N ) ) THEN
         info = -4
    ELSE IF( LDB.LT.max( 1, N ) ) THEN
         info = -6
    END IF
    IF( info.NE.0 ) THEN
       CALL xerbla( 'GAUSSJ ', -info )
       RETURN
    END IF
    
   ipiv=0    
   do i=1,N
     largest=0     
     do j=1,N
       if(ipiv(j) == 1) then 
        cycle
       endif     
      do k=1,N
       if (ipiv(k) > 1) then
        info=ipiv(k)
        write(*,*) 'Singular matrix in GaussJordan'   
        RETURN
       endif
       if(ipiv(k) == 1) then 
        cycle
       endif 
       if(largest >= ABS(A(j,k))) then 
        exit
       endif                        
       index_row=j
       index_col=k
       largest=ABS(A(j,k))
      end do
    end do
    ipiv(index_col)=ipiv(index_col)+1
    irow(i)=index_row     
    icol(i)=index_col
    if (index_row /= index_col) then    
      swap(1:N)=A(index_row,1:N)
      A(index_row,1:N)=A(index_col,1:N)
      A(index_col,1:N)=swap(1:N)
      swap2(1:NRHS)=B(index_row,1:NRHS)
      B(index_row,1:NRHS)=B(index_col,1:NRHS)
      B(index_col,1:NRHS)=swap2(1:NRHS)    
    endif
    pivot=A(index_col,index_col) 
    A(index_col,index_col)=1
    A(index_col,1:N)=A(index_col,1:N)/pivot       
    B(index_col,1:NRHS)=B(index_col,1:NRHS)/pivot
    do jj=1,N
     if (jj == index_col) then
      cycle
     endif
     temp=A(jj,index_col)
     A(jj,index_col)=0 
     A(jj,1:N)=A(jj,1:N)-A(index_col,1:N)*temp
     B(jj,1:NRHS)=B(jj,1:NRHS)-B(index_col,1:NRHS)*temp
    end do        
   end do 
   ii=1
   do i=1,N
    ii=N-i+1
    if (irow(ii) == icol(ii)) then
     cycle
    endif  
    index_row=irow(ii)     
    index_col=icol(ii)   
    swap(1:N)=A(1:N,index_row)
    A(1:N,index_row)=A(1:N,index_col)
    A(1:N,index_col)=swap(1:N)        
   end do  
            
  END SUBROUTINE GaussJordan     

  SUBROUTINE GaussJordan( N, NRHS, A, LDA, B, LDB, INFO )
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
              
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
  INTEGER ::  ipiv(N+1),icol(N),irow(N),i,j,k,ii,jj,index_row,index_col
  REAL(wp) :: largest,swap,temp,det,pivot
     
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
   det=1
   do i=1,N
     largest=0     
     do j=1,N
       if(ipiv(j) == 1) then 
        exit
       endif     
      do k=1,N
       if (ipiv(k) > 1) then
        info=ipiv(k)
        write(*,*) 'Singular matrix in GaussJordan', info,i,k,det       
        RETURN
       endif
       if(ipiv(k) == 1) then 
        exit
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
     det=-det
     do ii=1,N
      swap=A(index_row,ii)
      A(index_row,ii)=A(index_col,ii)
      A(index_col,ii)=swap
     end do    
     if (NRHS >= 1) then 
      do ii=1,NRHS
       swap=B(index_row,ii)
       B(index_row,ii)=B(index_col,ii)
       B(index_col,ii)=swap    
      end do
     endif        
    endif
    pivot=A(index_col,index_col)  
    det=det*pivot
    A(index_col,index_col)=1
    do ii=1,N
      A(index_col,ii)=A(index_col,ii)/pivot
    end do        
    if (NRHS >= 1) then     
     do ii=1,NRHS
       B(index_col,ii)=B(index_col,ii)/pivot
     end do         
    endif
!    jj=1
!    do while (jj /= index_col .AND. jj <= N)
     do jj=1,N
     if (jj == index_col) then
      exit
     endif
     temp=A(jj,index_col)
     A(jj,index_col)=0 
      do ii=1,N
        A(jj,ii)=A(jj,ii)-A(index_col,ii)*temp
      end do        
      if (NRHS >= 1) then     
        do ii=1,NRHS
         B(jj,ii)=B(jj,ii)-B(index_col,ii)*temp
        end do         
      endif        
!     jj=jj+1
    end do        
   end do 
   write(*,*) 'never got here'
   i=1
   ii=1
   do while (irow(ii) /= icol(ii) .AND. i <= N)
    ii=N-i+1  
    index_row=irow(ii)     
    index_col=icol(ii)   
    do k=1,N
     swap=A(k,index_row)
     A(k,index_row)=A(k,index_col)
     A(k,index_col)=swap
    end do 
    i=i+1         
   end do
   do i=1,N   
     if (ipiv(i) /= 1) then
      info=ipiv(i)
      write(*,*) 'Singular matrix in GaussJordan, det =',det      
      RETURN
     endif 
   end do    
            
  END SUBROUTINE GaussJordan     

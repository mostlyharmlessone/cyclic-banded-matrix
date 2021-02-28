    program testme
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=1001           ! size of problem
    INTEGER, PARAMETER :: KU=3             ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    REAL(wp) :: d(n,2),a(n),b(n),c(n),aij(n,n),bij(n,n),cij(n,n),s(n,2)
    REAL(wp) :: AB(2*KU+1,n),time_end,time_start
    INTEGER :: i,j,k,INFO,ipiv(n)
!   Copyright (c) 2021   Anthony M de Beus
    AB=0
   
!   This is an example of a band matrix for testing
    aij=0
    do i=1,n
     do j=1,n
      if (ABS(i - j) <= KU ) then  ! regular band
        aij(i,j)=500+3*i-2*j
        if (i > j) then
        aij(i,j)=aij(i,j)+1        ! assymetry
       endif         
      endif
      if (ABS(i - j) >= N-KU) then ! cyclic elements
       aij(i,j) = 20+i+2*j
       if (i > j) then
        aij(i,j)=aij(i,j)*2        ! assymetry
       endif
      endif
     end do
      s(i,1)=i                     ! solution vectors
      s(i,2)=i**2
    end do
    d=matmul(aij,s )               ! RHS vectors
    cij=aij                        ! store a copy

    call CPU_TIME(time_start)
    call DGESV(N, 2, AIJ, N, IPIV, D, N, INFO ) ! overwrites  aij & d
    call CPU_TIME(time_end)
    write(*,*) 'Using full matrix solver from LAPACK, dgesv.f, O(n^2)'
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '
    
    aij=cij
    d=matmul(aij,s)                           ! tests using two solution vectors
    call CPU_TIME(time_start)
!    call GaussJordan( N, 2, AIJ, N, D, N, INFO )
    call CPU_TIME(time_end)
    write(*,*) 'Using full matrix solver Gauss-Jordan with inverse, O(n^3)'
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    bij=aij
       
    aij=cij

    aij=matmul(aij,bij)
    bij=0
    do i=1,n
      bij(i,i) = 1.0
    enddo
    write(*,*) 'inverse error ',sum((aij-bij)*(aij-bij))   ! tests the inverse
    write(*,*) ' '

    aij=cij
    d=matmul(aij,s)

    do i=1,n                                           ! store aij in band-cyclic format
     do j=1,n
      if (aij(i,j) /= 0) then
      k=2*KU+2-mod(N+KU+1+i-j,N)
      AB(k,i)=aij(i,j)
      endif 
     end do
    end do

    IF (KU == 1) then                                  ! store tridiagonal matrices in vector format
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)

    call CPU_TIME(time_start)
    call DCTSV( n, 2, a, b, c, d, n, INFO )            ! overwrites d into solution
    call CPU_TIME(time_end)
    write(*,*) 'Using dctsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '

    ENDIF

    d=matmul(aij,s)
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using dcbsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
   

   END PROGRAM testme

    program testme
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=1000           ! size of problem
    INTEGER, PARAMETER :: KU=1             ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    REAL(wp) :: d(n),a(n),b(n),c(n),aij(n,n),s(n)
    REAL(wp) :: AB(2*KU+1,n),time_end,time_start
    INTEGER :: i,j,k,INFO,ipiv(n)
!   Copyright (c) 2021   Anthony M de Beus
    AB=0
   
!   This is a made up example of a band matrix for testing
    aij=0
    do i=1,n
     do j=1,n
      if (ABS(i - j) <= KU ) then  ! regular band
        aij(i,j)=500+3*i-2*j
        if (i > j) then
        aij(i,j)=aij(i,j)+1
       endif         
      endif
      if (ABS(i - j) >= N-KU) then
       aij(i,j) = 20+i+2*j
       if (i > j) then
        aij(i,j)=aij(i,j)*2
       endif
      endif
     end do
      s(i)=i
    end do
    d=matmul(aij,s)
    call CPU_TIME(time_start)
    call DGESV(N, 1, AIJ, N, IPIV, D, N, INFO ) ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using full matrix solver from LAPACK, DGESV.f'
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'error',dot_product((s-d),(s-d))
    write(*,*) ' '
    
!   This is a made up example of a band matrix for testing    
    aij=0
    do i=1,n
     do j=1,n
      if (ABS(i - j) <= KU ) then  ! regular band
        aij(i,j)=500+3*i-2*j
       if (i > j) then
        aij(i,j)=aij(i,j)+1
       endif             
      endif
      if (ABS(i - j) >= (N-KU)) then
       aij(i,j) = 20+i+2*j
       if (i > j) then
        aij(i,j)=aij(i,j)*2
       endif
      endif
     end do
      s(i)=i
    end do
    d=matmul(aij,s)

    do i=1,n
     do j=1,n
      if (aij(i,j) /= 0) then
      k=2*KU+2-mod(N+KU+1+i-j,N)
      AB(k,i)=aij(i,j)
      endif 
     end do
    end do

    IF (KU == 1) then
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)


    call CPU_TIME(time_start)
    call DCTSV( n, 1, a, b, c, d, n, INFO ) ! overwrites d into solution
    call CPU_TIME(time_end)
    write(*,*) 'Using dgtsv.f90'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'error',dot_product((s-d),(s-d))
    write(*,*) ' '

    ENDIF

    d=matmul(aij,s)
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 1, AB, 2*KU+1, d, N, INFO )  ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using dgtsv.f90'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'error',dot_product((s-d),(s-d))
   

   END PROGRAM testme

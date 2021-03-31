    program testme_large
!   doesn't compare with full matrix routine to allow for larger n    
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=901            ! size of problem
    INTEGER, PARAMETER :: KU=5             ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    INTEGER, PARAMETER :: KL=5             ! for testing vs lapack version only
                                           ! KL=KU to run non-periodic version of matrix KL=0 runs periodic version
    REAL(wp) :: d(n,2),a(n),b(n),c(n),s(n,2)
    REAL(wp) :: AB(2*KU+1,n),time_end,time_start ! AB(2*KU+1,n) for dcbsv; ! AB(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: CD(2*KL+KU+1,n)                  ! CD(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: a_short(n-1)                     ! truncated a() for dgtsv
    INTEGER :: i,j,k,INFO,ipiv(n)
!   Copyright (c) 2021   Anthony M de Beus
    AB=0
    CD=0
!   only runs dctsv.f90 if KU=1
!   only runs dgbsv and/or dgtsv if KL > 0 (and in fact KL=KU); non-periodic lapack routines, KU=KL=1 for dgtsv
   
    do i=1,n                                           ! store AB in band-cyclic format
     do j=1,n
      k=2*KU+2-mod(N+KU+1+i-j,N)
      AB(k,i)=500+3*i-2*j
     end do
    end do

    if (KL > 0) then
    do i=1,n                                           ! store CD in lapack band format if we're testing non-cyclic
     do j=1,n
      CD(KL+KU+1+i-j,j)=aij(i,j) 
     end do
    end do
    endif

    do i=1,n
      s(i,1)=i                     ! solution vectors
      s(i,2)=i**2
    end do

!    d=mat_mul(AB,s)                ! needs matrix multiplication for cyclic stored matrices
    do k=1,2 
     do j=1,n               
      do i=-KU,KU
       d(j,k)=AB(KU+1+i,j)*s(i+j,k)
      end do
     end do
    end do

    IF (KU == 1) then                                  ! store tridiagonal matrices in vector format
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)
    do i=1,n-1
     a_short(i)=a(i+1)                                 !truncated "a" for dgtsv, don't have to truncate "c"
    end do

    call CPU_TIME(time_start)
    call DCTSV( n, 2, a, b, c, d, n, INFO )            ! overwrites d into solution
    call CPU_TIME(time_end)
    write(*,*) 'Using dctsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '
    
     if (KL > 0) then   ! or KL == 1
      d=matmul(aij,s)       
      call CPU_TIME(time_start)
      call DGTSV( n, 2, a_short, b, c, d, n, INFO )     ! overwrites d into solution
      call CPU_TIME(time_end)
      write(*,*) 'Using dgtsv, O(n)'    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
      write(*,*) ' '
     endif

    ENDIF

    d=matmul(aij,s)
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using dcbsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '

     if (KL > 0) then   
      d=matmul(aij,s)       
      call CPU_TIME(time_start)          
      call dgbsv( N, KL, KU, 2, CD, 2*KL+KU+1, IPIV, d, N, INFO ) ! overwrites d into solution
      call CPU_TIME(time_end)
      write(*,*) 'Using dgbsv, O(n)'    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
     endif
 

   END PROGRAM testme

FUNCTION Area_Circle(r)

IMPLICIT NONE
REAL :: Area_Circle
REAL, INTENT(IN) :: r

! Declare local constant Pi
REAL, PARAMETER :: Pi = 3.1415927

Area_Circle = Pi * r * r

END FUNCTION Area_Circle




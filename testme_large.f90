    program testme_large
!   doesn't compare with full matrix routine to allow for larger n    
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=50000 ! size of problem
    INTEGER, PARAMETER :: KU=1    ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    INTEGER, PARAMETER :: KL=1  ! for testing vs lapack version only
                                  ! KL=KU to run non-periodic version of matrix KL=0 runs periodic version
    REAL(wp) :: d(n,2),a(n),b(n),c(n),s(n,2),dd(n,2)
    REAL(wp) :: AB(2*KU+1,n),time_end,time_start ! AB(2*KU+1,n) for dcbsv; ! AB(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: CD(2*KL+KU+1,n)                  ! CD(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: a_short(n-1)                     ! truncated a() for dgtsv
    INTEGER :: i,j,k,INFO,ipiv(n)
!   Copyright (c) 2021   Anthony M de Beus
    AB=0
    CD=0
!   only runs dctsv.f90 if KU=1
!   only runs dgbsv and/or dgtsv if KL > 0 (and in fact KL=KU)
!   non-periodic lapack routines, KU=KL=1 for dgtsv
   
    do i=-KU,KU                                ! generate AB in band-cyclic format
     do j=1,n
      AB(i+KU+1,j)=20*(KU**2-(i)**2) !+ 5.*j/(1.0*n)
      if (i == 0) then
       AB(i+KU+1,j)=AB(i+KU+1,j)+100           ! emphasize diagonal dominance
      end if
      if (KL == KU) then
      if ((i+j-1) /= mod(N+i+j-1,N)) then
        AB(i+KU+1,j)=0                         ! off diagonal elements zero for non-cyclic tests
      endif 
      endif
     end do
    end do

    if (KL > 0) then
    do i=1,2*KU+1                       ! store CD in lapack band format if we're testing non-cyclic
     do j=1,n
      CD(i+KL,j)=AB(2*KU+2-i,j)          !CD(KL+KU+1+i-j,j) vs  AB(2*KU+2-mod(N+KU+1+i-j,N),i)
     end do                                 !needs transpose of i,j too?
    end do
    endif

    do i=1,n
      s(i,1)=i                     ! solution vectors
      s(i,2)=10*sin(5.0*i)         ! i**2 gets too ill-conditioned with large n
    end do

! needs matrix multiplication for cyclic stored matrices
    dd=0
    do k=1,2 
     do j=1,n               
      do i=-KU,KU
       dd(j,k)=dd(j,k)+AB(KU+1+i,j)*s(1+mod(N+i+j-1,N),k)
      end do
     end do
    end do
    
    d=dd
    
    IF (KU == 1) then                           ! store tridiagonal matrices in vector format
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)
    do i=1,n-1
     a_short(i)=a(i+1)                          !truncated "a" for dgtsv, don't have to truncate "c"
    end do

    call CPU_TIME(time_start)
    call DCTSV( n, 2, a, b, c, d, n, INFO )     ! overwrites d into solution
    call CPU_TIME(time_end)
    write(*,*) 'Using dctsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '
    
     if (KL > 0) then   ! or KL == 1 
      d=dd      
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

    d=dd
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using dcbsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '

     if (KL > 0) then    
      d=dd      
      call CPU_TIME(time_start)          
      call dgbsv( N, KL, KU, 2, CD, 2*KL+KU+1, IPIV, d, N, INFO ) ! overwrites d into solution
      call CPU_TIME(time_end)
      write(*,*) 'Using dgbsv, O(n)'    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
     endif
 
   END PROGRAM testme_large


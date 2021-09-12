
    program testme_large
    USE OMP_LIB
    use AB_matrix_fct, only : multiply, RowConvert, ColumnConvert
!   only runs with banded matrix routines to allow for larger n    
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=90000 ! size of problem
    INTEGER, PARAMETER :: KU=378 ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    INTEGER, PARAMETER :: KL=378  ! for testing vs lapack version only
                                ! KL=KU to run non-periodic version of matrix KL=0 runs periodic version
    REAL(wp), ALLOCATABLE :: d(:,:),a(:),b(:),c(:),s(:,:),dd(:,:),z(:,:),zz(:,:),a_short(:) ! a_short truncated a() for dgtsv    
    REAL(wp), ALLOCATABLE :: AB(:,:),CD(:,:)                ! AB(2*KU+1,n) for dcbsv; ! AB(KL+KU+1+i-j,j) for dgbsv 
                                                            ! CD(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: time_end,time_start                     
    INTEGER :: i,j,k,INFO
    INTEGER, ALLOCATABLE :: ipiv(:)
!   Copyright (c) 2021   Anthony M de Beus

!   only runs tridiagonal routines dctsv.f90, dgtsv.f90 or thomas.f90 if KU=1
!   only runs non-periodic lapack routines dgbsv and/or dgtsv if KL > 0 (and in fact KL=KU)   
   
    allocate(AB(2*KU+1,n),CD(2*KL+KU+1,n),d(n,2),a(n),b(n),c(n),s(n,2),dd(n,2),z(n,2),zz(n,2),a_short(n-1),ipiv(n))

    AB=0
    CD=0

    do i=-KU,KU                                ! generate AB in band-cyclic format
     do j=1,n
      AB(i+KU+1,j)=20.0*(i)**2 + 5.*j/(1.0*n)  ! without the second term, can be ill-conditioned, eg N=9, KU=2 
      if (i == 0) then
       AB(i+KU+1,j)=AB(i+KU+1,j)+300           ! emphasize diagonal dominance
      end if
      if (i >= KU) then
       AB(i+KU+1,j)=AB(i+KU+1,j)+1.5           ! asymmetry
      end if   
      if (KL == KU) then
      if ((i+j-1) /= mod(N+i+j-1,N)) then
        AB(i+KU+1,j)=0                         ! off diagonal elements zero for non-cyclic tests
      endif 
      endif
     end do
    end do
       
    if (KL == KU .AND. KL > 0 ) then
     CD=ColumnConvert(AB)                ! for LAPACK use make a dgbsv compatible matrix
    endif

    IF (N > 100) then    
     do i=1,n
      s(i,1)=57.3*cos(40.0*i)        ! solution vectors
      s(i,2)=10*sin(5.0*i)           ! i and i**2 get too ill-conditioned with large n
     end do
    else
     do i=1,n
      s(i,1)=i                       ! solution vectors
      s(i,2)=i**2 
     end do
    endif  
    
    dd=multiply(AB,s)                           ! RHS vectors
    d=dd
            
    IF (KU == 1) then                           ! store tridiagonal matrices in vector format
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)
    do i=1,n-1
     a_short(i)=a(i+1)                          !truncated "a" for dgtsv, don't have to truncate "c"
    end do

    write(*,*) 'Using dctsv.f90, O(n)'    
    time_start=omp_get_wtime()                  ! have to use this to get actual time with omp
    call DCTSV( n, 2, a, b, c, d, n, INFO )     ! overwrites d into solution
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '


!    LAPACK routine for non-cyclic system    
     if (KL > 0) then   ! and KL == KU == 1 
      d=dd
      write(*,*) 'Using dgtsv, O(n)'           
      time_start=omp_get_wtime()
      call DGTSV( n, 2, a_short, b, c, d, n, INFO )     ! overwrites d into solution
      time_end=omp_get_wtime()   
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
      write(*,*) ' '
      
!     simple tridiagonal algorithm: should be fastest with -O3 compilation
      a=AB(1,:)
      b=AB(2,:)
      c=AB(3,:)
      d=dd 
      write(*,*) 'Using thomas, O(n)'        
      time_start=omp_get_wtime()
      call thomas(a,b,c,d,z,n,2)                        ! overwrites b and d, output is z           
      time_end=omp_get_wtime()    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-z(:,1)),(s(:,1)-z(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-z(:,2)),(s(:,2)-z(:,2)))
      write(*,*) ' '
          
     endif

    ENDIF

    d=dd
    write(*,*) 'Using dcbsv.f90, O(n)'
    time_start=omp_get_wtime()
    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    z=multiply(AB,s(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1)) 
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))
    write(*,*) ' '
    
!    zz=d    ! save the solution
!    d=z-dd  !RHS error
!    write(*,*) 'iterative step'
!    time_start=omp_get_wtime()
!    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites dd
!    time_end=omp_get_wtime()    
!    write(*,*) 'time: ',time_end-time_start 
!    z=zz-d
!    write(*,*) 'solution error',dot_product((s(:,1)-z(:,1)),(s(:,1)-z(:,1))) 
!    write(*,*) 'solution error',dot_product((s(:,2)-z(:,2)),(s(:,2)-z(:,2))) 
    
           
!    LAPACK routine for non-cyclic system
     if (KL > 0) then    
      d=dd  
      write(*,*) 'Using dgbsv, O(n)'    
      time_start=omp_get_wtime()          
      call dgbsv( N, KL, KU, 2, CD, 2*KL+KU+1, IPIV, d, N, INFO ) ! overwrites d into solution
      time_end=omp_get_wtime()    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
     endif
          
     deallocate(AB,CD,d,a,b,c,s,dd,z,a_short,ipiv)


   END PROGRAM testme_large


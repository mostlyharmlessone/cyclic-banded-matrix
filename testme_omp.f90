    program testme_omp
    USE OMP_LIB
    USE lapackinterface
    use AB_matrix_fct, only : multiply, RowConvert, ColumnConvert, ABConvert
!   only runs with banded matrix routines to allow for larger n    
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=24!8592! 8950 ! size of problem
    INTEGER, PARAMETER :: KU=2!358 ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    INTEGER, PARAMETER :: KL=0   ! for testing vs lapack version only
                                ! KL=KU to run non-periodic version of matrix KL=0 runs periodic version
    REAL(wp), ALLOCATABLE :: d(:,:),a(:),b(:),c(:),s(:,:),dd(:,:),z(:,:),zz(:,:),a_short(:) ! a_short truncated a() for dgtsv    
    REAL(wp), ALLOCATABLE :: AB(:,:),CD(:,:)                ! AB(2*KU+1,n) for dcbsv; ! AB(KL+KU+1+i-j,j) for dgbsv 
                                                            ! CD(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: time_end,time_start                     
    INTEGER :: i,j,p,INFO
    INTEGER, ALLOCATABLE :: ipiv(:)
!   Copyright (c) 2021   Anthony M de Beus

!   only runs tridiagonal routines dctsv.f90, dgtsv.f90 or thomas.f90 if KU=1
!   only runs non-periodic lapack routines dgbsv and/or dgtsv if KL > 0 (and in fact KL=KU) 

    INTERFACE  
     SUBROUTINE thomas(a,b,c,d,z,n,k) 
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
      INTEGER, INTENT(IN) :: n,k
      REAL(wp), INTENT(INOUT) :: a(n),b(n),c(n),d(n,k)
      REAL (wp), INTENT(OUT) :: z(n,k)    
     END SUBROUTINE thomas
     
     SUBROUTINE GaussJordan( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, N, NRHS
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
     END SUBROUTINE GaussJordan          
   
     SUBROUTINE DCTSV( N, NRHS, DL, D, DU, B, LDB, INFO )   
      INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision    
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: LDB, N, NRHS
      INTEGER, INTENT(OUT) :: INFO 
!     .. Array Arguments ..
      REAL(wp), INTENT(IN) :: D( * ), DL( * ), DU( * )  ! no output no LU factors
      REAL(wp), INTENT(INOUT) :: B( LDB, * ) ! on entry RHS, on exit, solution
     END SUBROUTINE DCTSV
      
     SUBROUTINE DCBSV( N, KU, NRHS, AB, LDAB, B, LDB, INFO ) 
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!       .. Scalar Arguments ..
       Integer, Intent(IN) ::  KU, LDAB, LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO
!        .. Array Arguments ..
       Real(wp), Intent(IN) :: AB( ldab, * )
       Real(wp), Intent(INOUT) ::  B( ldb, * )
     END SUBROUTINE DCBSV

     SUBROUTINE DCBSV_F( N, KU, NRHS, AB, LDAB, B, LDB, INFO ) 
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!       .. Scalar Arguments ..
       Integer, Intent(IN) ::  KU, LDAB, LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO
!        .. Array Arguments ..
       Real(wp), Intent(IN) :: AB( ldab, * )
       Real(wp), Intent(INOUT) ::  B( ldb, * )
     END SUBROUTINE DCBSV_F
     
     SUBROUTINE DCBSV_R( N, KU, NRHS, AB, LDAB, B, LDB, INFO ) 
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!       .. Scalar Arguments ..
       Integer, Intent(IN) ::  KU, LDAB, LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO
!        .. Array Arguments ..
       Real(wp), Intent(IN) :: AB( ldab, * )
       Real(wp), Intent(INOUT) ::  B( ldb, * )
     END SUBROUTINE DCBSV_R     
      
     SUBROUTINE DCBSV_4( N, KU, NRHS, AB, LDAB, B, LDB, INFO ) 
       INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
!       .. Scalar Arguments ..
       Integer, Intent(IN) ::  KU, LDAB, LDB, N, NRHS
       INTEGER, INTENT(OUT) :: INFO
!        .. Array Arguments ..
       Real(wp), Intent(IN) :: AB( ldab, * )
       Real(wp), Intent(INOUT) ::  B( ldb, * )
     END SUBROUTINE DCBSV_4                   
    END INTERFACE  
   
    allocate(AB(2*KU+1,n),CD(2*KL+KU+1,n),d(n,4),a(n),b(n),c(n),s(n,4),dd(n,4),z(n,4),zz(n,4),a_short(n-1),ipiv(n))

    AB=0
    CD=0
    p=mod(N,2*KU)
    write(*,*) 'p,KU,(N-p)/2KU: ',p,KU,(n-p)/(2*KU)

    do i=-KU,KU                                ! generate AB in band-cyclic format
     do j=1,n
      AB(i+KU+1,j)=20.0*(i)**2 + 5.*j/(1.0*n)  ! without the second term, can be ill-conditioned, eg N=9, KU=2 
      if (i == 0) then
       AB(i+KU+1,j)=AB(i+KU+1,j)+300           ! emphasize diagonal dominance
      end if

      AB(i+KU+1,j)=KU*KU-I*I+2+j

      if (i >= KU) then
       AB(i+KU+1,j)=AB(i+KU+1,j)+1!.5           ! asymmetry
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

!  Reverse solution vectors for alternate solutions
    s(1:n,3:4)=s(1:n,1:2)
    s(1:n,3:4) = s( size(s,1):1:-1,3:4)    
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
    call DCTSV( n, 4, a, b, c, d, n, INFO )     ! overwrites d into solution
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
      call DGTSV( n, 4, a_short, b, c, d, n, INFO )     ! overwrites d into solution
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
      call thomas(a,b,c,d,z,n,4)                        ! overwrites b and d, output is z           
      time_end=omp_get_wtime()    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-z(:,1)),(s(:,1)-z(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-z(:,2)),(s(:,2)-z(:,2)))
      write(*,*) ' '
          
     endif

    ENDIF

    d=dd
    write(*,*) 'Using dcbsv_f.f90, O(n)'
    time_start=omp_get_wtime()
    call DCBSV_F( N, KU, 4, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1)))/dot_product(s(:,1),s(:,1)) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))/dot_product(s(:,2),s(:,2))
    write(*,*) 'solution error',dot_product((s(:,3)-d(:,3)),(s(:,3)-d(:,3)))/dot_product(s(:,3),s(:,3)) 
    write(*,*) 'solution error',dot_product((s(:,4)-d(:,4)),(s(:,4)-d(:,4)))/dot_product(s(:,4),s(:,4))
    z=multiply(AB,d(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1))/dot_product(dd(:,1),dd(:,1))
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))/dot_product(dd(:,2),dd(:,2))
    write(*,*) 'RHS error',dot_product(z(:,3)-dd(:,3),z(:,3)-dd(:,3))/dot_product(dd(:,3),dd(:,3)) 
    write(*,*) 'RHS error',dot_product(z(:,4)-dd(:,4),z(:,4)-dd(:,4))/dot_product(dd(:,4),dd(:,4))
    write(*,*) ' '
   
    d=dd
    write(*,*) 'Using dcbsv_r.f90, O(n)'
    time_start=omp_get_wtime()
    call DCBSV_R( N, KU, 4, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1)))/dot_product(s(:,1),s(:,1)) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))/dot_product(s(:,2),s(:,2))
    write(*,*) 'solution error',dot_product((s(:,3)-d(:,3)),(s(:,3)-d(:,3)))/dot_product(s(:,3),s(:,3)) 
    write(*,*) 'solution error',dot_product((s(:,4)-d(:,4)),(s(:,4)-d(:,4)))/dot_product(s(:,4),s(:,4))
    z=multiply(AB,d(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1))/dot_product(dd(:,1),dd(:,1))
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))/dot_product(dd(:,2),dd(:,2))
    write(*,*) 'RHS error',dot_product(z(:,3)-dd(:,3),z(:,3)-dd(:,3))/dot_product(dd(:,3),dd(:,3)) 
    write(*,*) 'RHS error',dot_product(z(:,4)-dd(:,4),z(:,4)-dd(:,4))/dot_product(dd(:,4),dd(:,4))
    write(*,*) ' '
    
     d=dd
    write(*,*) 'Using dcbsv.f90, O(n)'
    time_start=omp_get_wtime()
    call DCBSV( N, KU, 4, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1)))/dot_product(s(:,1),s(:,1)) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))/dot_product(s(:,2),s(:,2))
    write(*,*) 'solution error',dot_product((s(:,3)-d(:,3)),(s(:,3)-d(:,3)))/dot_product(s(:,3),s(:,3)) 
    write(*,*) 'solution error',dot_product((s(:,4)-d(:,4)),(s(:,4)-d(:,4)))/dot_product(s(:,4),s(:,4))
    z=multiply(AB,d(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1))/dot_product(dd(:,1),dd(:,1))
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))/dot_product(dd(:,2),dd(:,2))
    write(*,*) 'RHS error',dot_product(z(:,3)-dd(:,3),z(:,3)-dd(:,3))/dot_product(dd(:,3),dd(:,3)) 
    write(*,*) 'RHS error',dot_product(z(:,4)-dd(:,4),z(:,4)-dd(:,4))/dot_product(dd(:,4),dd(:,4))
    write(*,*) ' '   

     d=dd
    write(*,*) 'Using dcbsv4.f90, O(n)'
    time_start=omp_get_wtime()
    call DCBSV_4( N, KU, 4, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    time_end=omp_get_wtime()    
    write(*,*) 'time: ',time_end-time_start

    write(*,*) 'solution ',s(:,1)
    write(*,*) ' '
    write(*,*) 'solution ',d(:,1)

    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1)))/dot_product(s(:,1),s(:,1)) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))/dot_product(s(:,2),s(:,2))
    write(*,*) 'solution error',dot_product((s(:,3)-d(:,3)),(s(:,3)-d(:,3)))/dot_product(s(:,3),s(:,3)) 
    write(*,*) 'solution error',dot_product((s(:,4)-d(:,4)),(s(:,4)-d(:,4)))/dot_product(s(:,4),s(:,4))
    z=multiply(AB,d(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1))/dot_product(dd(:,1),dd(:,1))
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))/dot_product(dd(:,2),dd(:,2))
    write(*,*) 'RHS error',dot_product(z(:,3)-dd(:,3),z(:,3)-dd(:,3))/dot_product(dd(:,3),dd(:,3)) 
    write(*,*) 'RHS error',dot_product(z(:,4)-dd(:,4),z(:,4)-dd(:,4))/dot_product(dd(:,4),dd(:,4))
    write(*,*) ' ' 
      
           
!    LAPACK routine for non-cyclic system
     if (KL > 0) then    
      d=dd  
      write(*,*) 'Using dgbsv, O(n)'    
      time_start=omp_get_wtime()          
      call dgbsv( N, KL, KU, 4, CD, 2*KL+KU+1, IPIV, d, N, INFO ) ! overwrites d into solution
      time_end=omp_get_wtime()    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1)))/dot_product(s(:,1),s(:,1)) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))/dot_product(s(:,2),s(:,2))
      write(*,*) 'solution error',dot_product((s(:,3)-d(:,3)),(s(:,3)-d(:,3)))/dot_product(s(:,3),s(:,3)) 
      write(*,*) 'solution error',dot_product((s(:,4)-d(:,4)),(s(:,4)-d(:,4)))/dot_product(s(:,4),s(:,4))
      z=multiply(AB,d(:,:))
      write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1))/dot_product(dd(:,1),dd(:,1))
      write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))/dot_product(dd(:,2),dd(:,2))
      write(*,*) 'RHS error',dot_product(z(:,3)-dd(:,3),z(:,3)-dd(:,3))/dot_product(dd(:,3),dd(:,3)) 
      write(*,*) 'RHS error',dot_product(z(:,4)-dd(:,4),z(:,4)-dd(:,4))/dot_product(dd(:,4),dd(:,4))
      write(*,*) ' '
     endif
          
     deallocate(AB,CD,d,a,b,c,s,dd,z,a_short,ipiv)


   END PROGRAM testme_omp


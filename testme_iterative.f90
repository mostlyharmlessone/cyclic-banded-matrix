module AB_matrix_fct

contains 

!z=matmul(aij,s)  z=multiply(AB,s(:,:))  ! matmul for AB matrix

function multiply(AB,s) result (z)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL (wp) :: AB(:,:),s(:,:) 
 REAL(wp) :: z(size(s,1),size(s,2))
 INTEGER :: i,j,k
 INTEGER :: n, KU
   KU=(size(AB,1)-1)/2
   n=size(AB,2)
   z=0
    do j=1,n 
      do i=1,2*KU+1
       z(j,:)=z(j,:)+AB(i,j)*s(1+mod(N+i+j-2-KU,N),:)                       
      end do
   end do
 end function
 
 function convert(AB) result (CD)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL (wp) :: AB(:,:)
 REAL(wp) :: CD((size(AB,1)+1)/2,size(AB,2))
 INTEGER :: i,j,m
 INTEGER :: n, KU
   m=size(AB,1)
   KU=(m-1)/2
   n=size(AB,2)
   CD=0
    do j=1,n 
      do i=1,2*KU+1  
       CD(i+KU,j)=AB(m-i+1,1+mod(N+i+j+KU-m-1,n))              
      end do
   end do
 end function

end module


    program testme
    use AB_matrix_fct, only : multiply, convert  
    IMPLICIT NONE
    INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
    INTEGER, PARAMETER :: n=1701          ! size of problem
    INTEGER, PARAMETER :: KU=23          ! bandwidth of matrix, KU=1 for dctsv.f90  KU>1 needs dcbsv.f90
    INTEGER, PARAMETER :: KL=23          ! for testing vs lapack version only
                                          ! KL=KU to run non-periodic version of matrix KL=0 runs periodic version
    REAL(wp) :: d(n,2),a(n),b(n),c(n),aij(n,n),bij(n,n),cij(n,n),s(n,2),z(n,2),dd(n,2)
    REAL(wp) :: AB(2*KU+1,n),time_end,time_start ! AB(2*KU+1,n) for dcbsv; ! AB(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: CD(2*KL+KU+1,n)                  ! CD(KL+KU+1+i-j,j) for dgbsv
    REAL(wp) :: a_short(n-1)                     ! truncated a() for dgtsv
    INTEGER :: i,j,k,INFO,ipiv(n)
!   Copyright (c) 2021   Anthony M de Beus

!   only runs dctsv.f90 if KU=1
!   only runs gauss-jordan if N < 1000 (gets too slow)
!   only runs dgbsv and/or dgtsv if KL > 0 (and in fact KL=KU); non-periodic lapack routines, KU=KL=1 for dgtsv
   
!   This is an example of a band matrix for testing
    aij=0
    do i=1,n
     do j=1,n
      if (ABS(i - j) <= KU ) then  ! regular band
        aij(i,j)=500+3.1*i-2*j
        if (i > j) then
        aij(i,j)=aij(i,j)+1        ! asymmetry
       endif         
      endif
      if (KL == 0) then             ! add the off diagonal terms in KL == 0
       if (ABS(i - j) >= N-KU) then ! cyclic elements
        aij(i,j) = 20+i+2*j
        if (i > j) then
         aij(i,j)=aij(i,j)*2        ! asymmetry
        endif
       endif
      endif
     end do
     end do  
      
    IF (N > 500) then    
     do i=1,n
      s(i,1)=57.3*cos(40.0*i)        ! solution vectors
      s(i,2)=10*sin(5.0*i)           ! i and i**2 get too ill-conditioned with large n
     end do
    else
     do i=1,n
      s(i,1)=i                     ! solution vectors
      s(i,2)=i**2 
     end do      
    endif
                      
    d=matmul(aij,s )               ! RHS vectors
    cij=aij                        ! store a copy
    dd=d

    call CPU_TIME(time_start)
    call DGESV(N, 2, AIJ, N, IPIV, D, N, INFO ) ! overwrites  aij & d
    call CPU_TIME(time_end)
    write(*,*) 'Using full matrix solver from LAPACK, dgesv.f, O(n^2)'
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '
 
    if (N < 1000) then   ! don't do this unless you want to wait a long time
    aij=cij
    d=dd                          
    call CPU_TIME(time_start)
    call GaussJordan( N, 2, AIJ, N, D, N, INFO )
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
    endif

    aij=cij
    d=dd
    AB=0

!     do i=1,n                                           ! store aij in band-cyclic format
!     do j=1,n
!      if (aij(i,j) /= 0) then
!      k=2*KU+2-mod(N+KU+1+i-j,N)
!      AB(k,i)=aij(i,j)
!      endif 
!     end do 
!    end do
    
    do i=1,1+2*KU                                          ! store aij in band-cyclic format
     do j=1,n        
      k=1+mod(N+i+j-2-KU,N)
      AB(i,j)=aij(j,k) 
     end do 
    end do  
     
    CD=0
    
!    if (KL > 0) then
!      do i=1,n                             ! store aij in lapack band format if we're testing non-cyclic
!      do j=1,n
!       if (aij(i,j) /= 0) then
!       CD(KL+KU+1+i-j,j)=aij(i,j)
!       endif 
!      end do
!     end do
!    endif
           
    if (KL > 0 .AND. KL==KU) then    
     do i=1,1+KL+KU                         ! store aij in lapack band format if we're testing non-cyclic
      do j=1,n        
       k=1+mod(N+i+j-2-KU,N)
       CD(i+KL,j)=aij(k,j) 
      end do 
     end do 
    endif 
              
    IF (KU == 1) then                                  ! store tridiagonal matrices in vector format
    a=AB(1,:)
    b=AB(2,:)
    c=AB(3,:)
    do i=1,n-1
     a_short(i)=a(i+1)                                 ! truncated "a" for dgtsv, don't have to truncate "c"
    end do

    call CPU_TIME(time_start)
    call DCTSV( n, 2, a, b, c, d, n, INFO )            ! overwrites d into solution
    call CPU_TIME(time_end)
    write(*,*) 'Using dctsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' '  

!    LAPACK routine for non-cyclic system    
     if (KL > 0) then   ! or KL == 1
      d=dd      
      call CPU_TIME(time_start)
      call DGTSV( n, 2, a_short, b, c, d, n, INFO )     ! overwrites b and d into solution
      call CPU_TIME(time_end)
      write(*,*) 'Using dgtsv, O(n)'    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
      write(*,*) ' '
     
!     simple tridiagonal algorithm: should be fastest with -O3 compilation
      a=AB(1,:)
      b=AB(2,:)
      c=AB(3,:)
      d=dd       
      call CPU_TIME(time_start)
      call thomas(a,b,c,d,z,n,2)                        ! overwrites b and d, output is z           
      call CPU_TIME(time_end)
      write(*,*) 'Using thomas, O(n)'    
      write(*,*) 'time: ',time_end-time_start
      write(*,*) 'solution error',dot_product((s(:,1)-z(:,1)),(s(:,1)-z(:,1))) 
      write(*,*) 'solution error',dot_product((s(:,2)-z(:,2)),(s(:,2)-z(:,2)))
      write(*,*) ' '

     endif
              
    ENDIF

!    LAPACK routine for non-cyclic system
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

!   put this last to use dd     
    d=dd
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 2, AB, 2*KU+1, d, N, INFO )      ! overwrites d
    call CPU_TIME(time_end)
    write(*,*) 'Using dcbsv.f90, O(n)'    
    write(*,*) 'time: ',time_end-time_start
    write(*,*) 'solution error',dot_product((s(:,1)-d(:,1)),(s(:,1)-d(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-d(:,2)),(s(:,2)-d(:,2)))
    write(*,*) ' ' 
    z=multiply(AB,s(:,:))
    write(*,*) 'RHS error',dot_product(z(:,1)-dd(:,1),z(:,1)-dd(:,1)) 
    write(*,*) 'RHS error',dot_product(z(:,2)-dd(:,2),z(:,2)-dd(:,2))
    
    d=z-dd
    write(*,*) 'iterative step'
    call CPU_TIME(time_start)
    call DCBSV( N, KU, 2, AB, 2*KU+1, dd, N, INFO )      ! overwrites dd
    call CPU_TIME(time_end)    
    write(*,*) 'time: ',time_end-time_start 
    z=d+dd
    write(*,*) dd
    write(*,*) 'solution error',dot_product((s(:,1)-z(:,1)),(s(:,1)-z(:,1))) 
    write(*,*) 'solution error',dot_product((s(:,2)-z(:,2)),(s(:,2)-z(:,2)))    
     
     
 

   END PROGRAM testme

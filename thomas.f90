subroutine thomas(a,b,c,d,z,n,k) ! "Llewellyn Thomas" algorithm for tridiagonal banded matrices"
! Adapted from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
integer, INTENT(IN) :: n,k
real(wp), INTENT(INOUT) :: a(n),b(n),c(n),d(n,k)
real (wp), INTENT(OUT) :: z(n,k)
real(wp) :: g
integer i
do i=2,n
 g=a(i)/b(i-1)
 b(i)= b(i)-g*c(i-1)
 d(i,:)=d(i,:)-g*d(i-1,:)
end do
 z(n,:)= d(n,:)/b(n)
do i=n-1,1,-1
 z(i,:)=(d(i,:)-c(i)*z(i+1,:))/b(i)
end do 
return
end subroutine thomas

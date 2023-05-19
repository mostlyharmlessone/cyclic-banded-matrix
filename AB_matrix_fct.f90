module AB_matrix_fct

contains 

!z=matmul(aij,s)  z=multiply(AB,s(:,:))  matmul for AB matrix

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
 
! makes an AB matrix with bandwidth KU from a square matrix Aij
 function ABConvert(KU,A) result (AB) 
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL (wp) :: A(:,:)
 REAL(wp) :: AB(2*KU+1,size(A,2))
 INTEGER :: KU,i,j,m,n
  n=size(A,2)
  m=2*KU+1
   do j=1,n   
    do i=1,m
     AB(i,j)=a(mod(N+i+j+KU-m-1,n),j)
    end do
   end do
 end function
 
! CD=ColumnConvert(AB)  makes a LAPACK band matrix compatible with dgbsv and KU=KL out of periodic AB matrix
  
 function ColumnConvert(AB) result (CD)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL (wp) :: AB(:,:)
 REAL(wp) :: CD((3*size(AB,1)-1)/2,size(AB,2))
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
 
 ! AB=RowConvert(CD)  makes a periodic AB matrix out of a LAPACK band matrix compatible with dgbsv and KU=KL 
  
 function RowConvert(CD) result (AB)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL (wp) :: CD(:,:)
 REAL(wp) :: AB((2*size(CD,1)+1)/3,size(CD,2))
 INTEGER :: i,j,m
 INTEGER :: n, KU
   m=(2*size(CD,1)+1)/3
   KU=(m-1)/2
   n=size(CD,2)
   AB=0
    do j=1,n 
      do i=1,2*KU+1  
       AB(m-i+1,1+mod(N+i+j+KU-m-1,n))=CD(i+KU,j)             
      end do
   end do
 end function

 ! Rotates columns of a matrix

 function RotateColumns(A,k) result(B)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL(wp) :: A(:,:)
 REAL(wp) :: B(size(A,1),size(A,2))
 INTEGER :: k,kk
  kk = mod(k,size(A,1))
  B = A(size(A,1):1:-1,1:size(A,2))
  B(1:size(A,1)-kk,1:size(A,2)) = B(size(A,1)-kk:1:-1,1:size(A,2))
  B(size(A,1)-kk+1:size(A,1),1:size(A,2)) = B(size(A,1):size(A,1)-kk+1:-1,1:size(A,2))
 end function 

 function ReverseColumns(A) result(B)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL(wp) :: A(:,:)
 REAL(wp) :: B(size(A,1),size(A,2))
  B = A(size(A,1):1:-1,1:size(A,2))
 end function 

 function ReverseRotation(A) result(B)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL(wp) :: A(:,:)
 REAL(wp) :: B(size(A,1),size(A,2)),C(size(A,1)/4,size(A,2))
  B = A(size(A,1):1:-1,1:size(A,2)) 
  C(1:(size(A,1)/4),size(A,2)) = B(1:(size(A,1)/4),size(A,2))
  B(1:(size(A,1)/4),size(A,2)) = B((size(A,1)/4+1):size(A,1)/2,size(A,2))
  B((size(A,1)/4+1):size(A,1)/2,size(A,2)) = C(1:(size(A,1)/4),size(A,2))
  C(1:(size(A,1)/4),size(A,2)) = B(size(A,1)/2+1:(3*size(A,1)/4),size(A,2))
  B(size(A,1)/2+1:(3*size(A,1)/4),size(A,2)) = B((3*size(A,1)/4+1):size(A,1),size(A,2))
  B((3*size(A,1)/4+1):size(A,1),size(A,2)) = C(1:(size(A,1)/4),size(A,2))
 end function

 ! Rotates rows of a matrix

 function RotateRows(A,k) result(B)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL(wp) :: A(:,:)
 REAL(wp) :: B(size(A,1),size(A,2))
 INTEGER :: k,kk
  kk = mod(k,size(A,2))
  B = A(1:size(A,1),size(A,2):1:-1)
  B(1:size(A,1),1:size(A,2)-kk) = B(1:size(A,1),size(A,2)-kk:1:-1)
  B(1:size(A,1),size(A,2)-kk+1:size(A,2)) = B(1:size(A,1),size(A,2):size(A,2)-kk+1:-1)
 end function

 function ReverseRows(A) result(B)
 implicit none
 INTEGER, PARAMETER :: wp = KIND(0.0D0) ! working precision
 REAL(wp) :: A(:,:)
 REAL(wp) :: B(size(A,1),size(A,2))
  B = A(1:size(A,1),size(A,2):1:-1)
 end function

end module

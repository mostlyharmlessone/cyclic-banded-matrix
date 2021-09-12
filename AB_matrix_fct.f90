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

end module

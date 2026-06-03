      MODULE LapackInterface
      IMPLICIT NONE
!  -- interfaces for FORTRAN77 LAPACK to fortran 90+      
!  -- suggested from http://www.siam.org/books/ot134 Numerical Computing with Modern Fortran Richard J.Hanson and Tim Hopkins SIAM 
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016

        INTERFACE

         SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

!     .. Scalar Arguments ..
         CHARACTER          JOBZ, UPLO
         INTEGER            INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
         DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
         END SUBROUTINE DSYEV


         SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!     .. Scalar Arguments ..
         INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
         INTEGER            IPIV( * )
         DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
         END SUBROUTINE DGESV
         
         SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
!     .. Scalar Arguments ..
         CHARACTER          TRANS
         INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     .. Array Arguments ..
         DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
         END SUBROUTINE DGELS

          SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!     .. Scalar Arguments ..
          INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * )
          END SUBROUTINE DGETRF


          SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!     .. Scalar Arguments ..
          CHARACTER          TRANS
          INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
          END SUBROUTINE DGETRS


          SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!      .. Scalar Arguments ..
          INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!      ..
!      .. Array Arguments ..
          INTEGER            IPIV( * )
          DOUBLE PRECISION   AB( ldab, * ), B( ldb, * )
          
          END SUBROUTINE DGBSV
          
          SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )

!      .. Scalar Arguments ..
          INTEGER            INFO, LDB, N, NRHS
!      ..
!      .. Array Arguments ..
          DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )

          END SUBROUTINE DGTSV
                  

          SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!      .. Scalar Arguments ..
          INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
          INTEGER            IPIV( * )
          DOUBLE PRECISION   A( LDA, * ), WORK( * )
          END SUBROUTINE DGETRI

          SUBROUTINE XERBLA( SRNAME, INFO )
!         .. Scalar Arguments ..
          CHARACTER(6)        SRNAME
          INTEGER            INFO
          END SUBROUTINE XERBLA  

          SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!         .. Scalar Arguments ..
          DOUBLE PRECISION ALPHA,BETA
          INTEGER K,LDA,LDB,LDC,M,N
          CHARACTER TRANSA,TRANSB
!         .. Array Arguments ..
          DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*) 
          END SUBROUTINE DGEMM  

          SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!         .. Scalar Arguments ..
          DOUBLE PRECISION ALPHA,BETA
          INTEGER INCX,INCY,LDA,M,N
          CHARACTER TRANS
!         .. Array Arguments ..
          DOUBLE PRECISION A(LDA,*),X(*),Y(*) 
          END SUBROUTINE DGEMV   

          FUNCTION DNRM2(N,DX,INCX) RESULT(RES)
          DOUBLE PRECISION RES
!          .. Scalar Arguments ..
          INTEGER  INCX,N
!         .. Array Arguments ..
          DOUBLE PRECISION DX(*)
          END FUNCTION DNRM2

        END INTERFACE    

      END MODULE LapackInterface

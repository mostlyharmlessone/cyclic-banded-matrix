# Cyclic-banded-matrix

Cyclic-banded-matrix is fortran 90+ source code for solving cyclic/periodic square matrices.  The algorithm is detailed in [band_matrix.pdf](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/band_matrix.pdf)  

I have not found any previous version of this algorithm and it is entirely my derivation. If I have reinvented the wheel, there being nothing new under the sun, it would be interesting to see the previous attribution.

The error rate increases more rapidly with N than an LU decomposition, but is O(N*KU) in speed.  The matrix needs to be diagonally dominant.

## Usage

Two routines are included, dctsv.f90 which solves periodic tridiagonal systems, and dcbsv.f90, which solves general banded periodic matrices.
They are written to be reasonably easily incorporated into code using [LAPACK (Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.November 2006)](http://www.netlib.org/lapack/) and the interfaces have been deliberatedly made as compatible as possible with their non-cyclic LAPACK counterparts dgbsv.f and dgtsv.f

Both use XERBLA error reporting. dcbsv.f90 requires linking with lapack and blas. dctsv.f90, with minimal editing to remove XERBLA, can stand alone.  More extensive editing of dgbsv.f90 with substitution of intrinsic Fortran operators instead of BLAS calls and substitution of LAPACK's dgesv.f by another suitable general matrix solver (e.g. the included gauss-jordan.f90) is also reasonably easy.

A fortran 90 module lapackinterface.f90 is included for compatibility with FORTRAN 77. Note that the module only needs to be compiled once.

A simple test program is included, 

Using gfortran
```
gfortran -c LapackInterface.f90
gfortran testme.f90 dctsv.f90 dcbsv.f90 gauss-jordan.f90 -llapack -lblas

or

gfortran testme.f90 LapackInterface.f90 dctsv.f90 dcbsv.f90 gauss-jordan.f90 -llapack -lblas

./a.out

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[LICENSE](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/LICENSE)


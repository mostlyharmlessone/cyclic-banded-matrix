# Cyclic-banded-matrix

Cyclic-banded-matrix is fortran 90+ source code for solving a system Ax=b with A represented by a cyclic band matrix, aka periodic band matrix, or cyclic/periodic banded matrix. These are all square matrices, and arise from problems with periodic boundary conditions.

The algorithm is detailed in [band_matrix.pdf](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/band_matrix.pdf).  I have not found any previous version of this algorithm and it is entirely my derivation. If I have reinvented the wheel, there being nothing new under the sun, it would be interesting to see the previous attribution.

The error rate increases more rapidly with N than an LU decomposition, but is O(N*KU) in speed.  The matrix needs to be diagonally dominant.

## Usage

Two routines are included, dctsv.f90 which solves periodic tridiagonal systems, and dcbsv.f90, which solves general banded periodic matrices.
They are written to be reasonably easily incorporated into code using [LAPACK (Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.November 2006)](http://www.netlib.org/lapack/) and the interfaces have been deliberately made as compatible as possible with their non-cyclic LAPACK counterparts dgbsv.f and dgtsv.f
Both use XERBLA error reporting and require linking with LAPACK & BLAS. 

Removing XERBLA from dctsv.f90 removes the LAPACK/BLAS dependency.

Removing the LAPACK/BLAS dependency in dcbsv requires substitution by intrinsic Fortran array operators instead of LAPACK & BLAS calls and substitution of LAPACK's dgesv.f  & dgetr(fsi).f by another suitable general matrix solver (e.g. the included gauss-jordan.f90).  The intrinsic operators and calls to the alternate matrix solver are commented out immediately after their respective LAPACK calls, making it easy to switch.

The algorithm can be run forwards or backwards, with dcbsv.f90, dcbsv_reverse.f90, or both simultaneously using OpenMP with dcbsv_parallel_omp.f90.   The parallel OMP versions are predictably nearly twice as fast as the non parallel versions.

A fortran 90 module LapackInterface.f90 is included for compatibility with FORTRAN 77. Note that the module only needs to be compiled once and is only needed for LAPACK/BLAS compatibility.

Two simple test programs are included, testme.f90 using full matrix routines as well as the LAPACK banded noncyclic routines dgtsv and dgbsv for comparison and testme_omp.f90 for larger N with a simple tridiagonal solver thomas.f90 for comparison purposes.

Using gfortran
```
gfortran -c LapackInterface.f90
gfortran testme.f90 dctsv.f90 dcbsv.f90 gauss-jordan.f90 -llapack -lblas

or

gfortran testme.f90 LapackInterface.f90 dctsv.f90 dcbsv.f90 gauss-jordan.f90 -llapack -lblas

or (after commenting out dgesv in testme.f90)

gfortran testme.f90 dctsv.f90 dcbsv.f90 gauss-jordan.f90


parallel versions with OpenMP

gfortran -fopenmp -O3 testme_omp.f90 dcbsv_parallel_omp.f90 dctsv_parallel_omp.f90 thomas.f90 -llapack -lblas

./a.out

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[LICENSE](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/LICENSE)


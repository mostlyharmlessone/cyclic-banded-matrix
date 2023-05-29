# Cyclic-banded-matrix

Cyclic-banded-matrix is fortran 90+ source code for solving a system Ax=b with A represented by a cyclic band matrix, aka periodic band matrix, or cyclic/periodic banded matrix. These are all square matrices, and arise from problems with periodic boundary conditions.

The algorithm is detailed in [band_matrix.pdf](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/band_matrix.pdf).  I have not found any previous version of this algorithm and it is entirely my derivation. If I have reinvented the wheel, there being nothing new under the sun, it would be interesting to see the previous attribution.

The error rate increases more rapidly with N than an LU decomposition, but is O(N*KU) in speed.  The matrix needs to be diagonally dominant.

## Usage

Two routines are included, dctsv.f90 which solves periodic tridiagonal systems, and dcbsv.f90, which solves general banded periodic matrices.
They are written to be reasonably easily incorporated into code using [LAPACK (Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.November 2006)](http://www.netlib.org/lapack/) and the interfaces have been deliberately made as compatible as possible with their non-cyclic LAPACK counterparts dgbsv.f and dgtsv.f

Both use XERBLA error reporting and require linking with LAPACK & BLAS. 
Removing XERBLA from dctsv.f90 removes the LAPACK/BLAS dependency.

Removing the LAPACK/BLAS dependency in dcbsv requires substitution by intrinsic Fortran array operators instead of LAPACK & BLAS calls and substitution of LAPACK's dgesv.f  & dgetr(fsi).f by another suitable general matrix solver (e.g. the included admittedly slower O(N^3) gauss-jordan.f90).  The intrinsic operators and calls to the alternate matrix solver are commented out immediately after their respective LAPACK calls, making it easy to switch.

The algorithm can be run forwards or backwards, with dcbsv_forward.f90, dcbsv_reverse.f90, or both simultaneously using OpenMP with dcbsv_parallel_omp.f90.   The parallel OMP versions can be nearly twice as fast as the non parallel versions depending on the parameters.

A fortran 90 module LapackInterface.f90 is included for compatibility with FORTRAN 77. Note that the module only needs to be compiled once and is only needed for LAPACK/BLAS compatibility.

Band matrices are defined per LAPACK convention https://www.netlib.org/lapack/lug/node124.html for non-periodic matrices where the **columns** of the full matrix are the columns of the banded form while the rows of the banded form are the diagonals of the full matrix.  For periodic banded matrices, the definition of the banded form is similar with the **rows** of the full matrix as the columns of the banded form while the rows of the banded form are the diagonals of the full matrix.  See [band_matrix.pdf](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/band_matrix.pdf)

  The band diagonal storage scheme is illustrated by the following example, when
  N = 9, KU = 2

     AB(2*KU+2-mod(N+KU+1+i-j,N),i)                  A(i,j)=aij  
     a18  a29  a31  a42  a53 a64 a75 a86 a97         a11 a12 a13  0   0   0   0 a18  a19
     a19  a21  a32  a43  a54 a65 a76 a87 a98         a21 a22 a23 a24  0   0   0  0   a29
     a11  a22  a33  a44  a55 a66 a77 a88 a99         a31 a32 a33 a34 a35  0   0  0    0
     a12  a23  a34  a45  a56 a67 a78 a89 a91         0   a42 a43 a44 a45 a46  0  0    0
     a13  a24  a35  a46  a57 a68 a79 a81 a92         0   0   a53 a54 a55 a56 a57 0    0
                                                     0   0   0   a64 a65 a66 a67 a68  0
                                                     0   0   0   0   a75 a76 a77 a78 a79
                                                     a81 0   0   0    0  a86 a87 a88 a89
                                                     a91 a92 0   0    0   0  a97 a98 a99
													 
													 
A fortran 90 module AB_matrix_fct.f90 to work with the defined periodic banded matrix and their non-cyclic LAPACK counterpart is included.  Note that the modules only need to be compiled once.  The main loops of the algorithms are included as forward_loop.f90 and backward_loop.f90.
 
Two simple test programs are included, testme.f90 using full matrix routines as well as the LAPACK banded noncyclic routines dgtsv and dgbsv for comparison and testme_omp.f90 for larger N with a simple tridiagonal solver thomas.f90 for comparison purposes. CityPlots https://math.nist.gov/MatrixMarket/ of the matrices can be generated with the included Mathematica(R) notebook CityPlot.nb

Using gfortran/ifort/nvfortran
```
gfortran -c LapackInterface.f90 AB_matrix_fct.f90
gfortran testme.f90 dctsv.f90 dcbsv_forward.f90 dcbsv_reverse.f90 forward_loop.f90 backward_loop.f90 gauss-jordan.f90 -llapack -lblas

parallel versions with OpenMP, different compilers

gfortran -fopenmp -O3 testme_omp.f90 dcbsv_parallel_omp.f90 dcbsv_forward.f90 dcbsv_reverse.f90 dctsv_parallel_omp.f90 forward_loop.f90 backward_loop.f90 thomas.f90 LapackInterface.f90 AB_matrix_fct.f90 -llapack -lblas

ifort -qopenmp -shared-intel testme_omp.f90 dcbsv_parallel_omp.f90 dcbsv_4x4parallel_omp.f90 dcbsv_forward.f90 dcbsv_reverse.f90 dctsv_parallel_omp.f90 forward_loop.f90 backward_loop.f90 thomas.f90 LapackInterface.f90 AB_matrix_fct.f90 gauss-jordan.f90 -llapack -lblas

nvfortran -mp testme_omp.f90 dcbsv_parallel_omp.f90 dcbsv_4x4parallel_omp.f90 dcbsv_forward.f90 dcbsv_reverse.f90 dctsv_parallel_omp.f90 forward_loop.f90 backward_loop.f90 thomas.f90 LapackInterface.f90 AB_matrix_fct.f90 gauss-jordan.f90 -llapack -lblas

./a.out

Alternatively, 
  cmake ./    
  make   
  ./testme_cyclic 
  ./testme_cyclic_omp
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[LICENSE](https://github.com/mostlyharmlessone/cyclic-banded-matrix/blob/main/LICENSE)


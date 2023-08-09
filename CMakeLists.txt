cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.16)

project(test_cyclic VERSION 0.1.0 LANGUAGES Fortran)
enable_language(Fortran)

set (CMAKE_Fortran_COMPILER "ifort")
set (CMAKE_Fortran_COMPILER "nvfortran")
set (CMAKE_Fortran_COMPILER "gfortran")

# should read again how to manage compiler flags with Cmake
if(CMAKE_Fortran_COMPILER MATCHES "gfortran")
    set(FC "gfortran")
    set(dialect "-ffree-form -std=f2008 -fexternal-blas -fno-optimize-sibling-calls") #gfortran/lapack/blas bug?
     set(bounds "-Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -fbounds-check -g")
     set(omp "-fopenmp")
     set(NaN "-ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999")
#  profiling with gprof 
 #     set(-pg)    
#  src profiling with gcov -a -b xxx.f90
#     set(-g -coverage) 
#    set(optimize "-O3 -march=native -ffast-math -funroll-loops")    
endif()
if(CMAKE_Fortran_COMPILER MATCHES "ifort")
    set(FC "ifort")
    set(dialect "-stand f08 -free -implicitnone")
    set(type "-shared-intel")
    set(bounds "-check bounds")
    set(omp "-qopenmp")   
endif()
#  nvfortran==pgi
if(CMAKE_Fortran_COMPILER MATCHES "nvfortran")
    set(FC "nvfortran")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-Mbounds -g")
    set(omp "-mp")    
endif()

#set(CMAKE_BUILD_TYPE RELEASE) # no debug or bounds check
set(CMAKE_BUILD_TYPE DEBUG) #this include debugging with addition of bounds flags as specified below
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_RELEASE} ${optimize}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${omp}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${NaN}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${type}")

set_property(GLOBAL PROPERTY GLOBAL_DEPENDS_DEBUG_MODE 0) #set to 1 to show dependencies

set (CMAKE_PREFIX_PATH "/usr/local")

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package(BLAS)
 if ( BLAS_FOUND )
  message(STATUS "BLAS found at ${BLAS_LIBRARIES}")
 else()
  set(BLAS_LIBRARIES "/usr/local/lib/libblas.so")
 endif()
find_package(LAPACK)
 if ( LAPACK_FOUND )
  message(STATUS "LAPACK found at ${LAPACK_LIBRARIES}")
 else()
  set(LAPACK_LIBRARIES "/usr/lib/i386-linux-gnu/lapack/liblapack.so.3")
 endif()
find_package(OpenMP)

# Compile.

#   Fortran Source
    file(GLOB mod_source  LapackInterface.f90 AB_matrix_fct.f90)

    file(GLOB source      dctsv.f90
                          dcbsv_forward.f90
                          dcbsv_reverse.f90 
                          forward_loop.f90
                          backward_loop.f90
                          gauss-jordan.f90)
                          
    file(GLOB source_omp  dcbsv_p4x4parallel_omp.f90 
                          dcbsv_4x4parallel_omp.f90
                          dcbsv_parallel_omp.f90  
                          dcbsv_forward.f90 
                          dcbsv_reverse.f90 
                          dctsv_parallel_omp.f90 
                          forward_loop.f90 
                          backward_loop.f90                      
                          gauss-jordan.f90 
                          thomas.f90)
                           
    add_library( modules ${mod_source})   
    add_library( source_omp ${source_omp})
    add_library( source ${source}) 
    add_dependencies(source_omp modules)
    add_dependencies(source modules)     

    set(LIBS modules ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} )

    add_executable(test_cyclic_omp testme_omp.f90)                                                         
    target_link_libraries(test_cyclic_omp ${LIBS} source_omp)
 
    add_executable(test_cyclic testme.f90)                                                        
    target_link_libraries(test_cyclic ${LIBS} source)      





#Source files grouped by a directory
GET_SOURCE_FILES(SOURCES_MATHLIB)
set(SOURCES ${SOURCES_MATHLIB})

GET_SOURCE_FILES(SOURCES_INTERPOLATIONALGORITHMS InterpolationAlgorithms)
set(SOURCES ${SOURCES} ${SOURCES_INTERPOLATIONALGORITHMS})

GET_SOURCE_FILES(SOURCES_INTEGRATION Integration)
set(SOURCES ${SOURCES} ${SOURCES_INTEGRATION})

GET_SOURCE_FILES(SOURCES_LINALG LinAlg)
set(SOURCES ${SOURCES} ${SOURCES_LINALG})

GET_SOURCE_FILES(SOURCES_LINALG_DENSE LinAlg/Dense)
set(SOURCES ${SOURCES} ${SOURCES_LINALG_DENSE})

GET_SOURCE_FILES(SOURCES_LINALG_SPARSE LinAlg/Sparse)
set(SOURCES ${SOURCES} ${SOURCES_LINALG_SPARSE})

GET_SOURCE_FILES(SOURCES_LINALG_SOLVERS LinAlg/Solvers)
set(SOURCES ${SOURCES} ${SOURCES_LINALG_SOLVERS})

GET_SOURCE_FILES(SOURCES_LINALG_PRECOND LinAlg/Preconditioner)
set(SOURCES ${SOURCES} ${SOURCES_LINALG_PRECOND})

GET_SOURCE_FILES(SOURCES_ODE ODE)
set(SOURCES ${SOURCES} ${SOURCES_ODE})

if(OGS_USE_EIGEN)
    GET_SOURCE_FILES(SOURCES_LINALG_EIGEN LinAlg/Eigen)
    set(SOURCES ${SOURCES} ${SOURCES_LINALG_EIGEN})
endif()

if(OGS_USE_LIS)
    GET_SOURCE_FILES(SOURCES_LINALG_LIS LinAlg/Lis)
    set(SOURCES ${SOURCES} ${SOURCES_LINALG_LIS})
endif()

if(OGS_USE_EIGEN AND OGS_USE_LIS)
    GET_SOURCE_FILES(SOURCES_LINALG_EIGENLIS LinAlg/EigenLis)
    set(SOURCES ${SOURCES} ${SOURCES_LINALG_EIGENLIS})
endif()

if(OGS_USE_PETSC)
    GET_SOURCE_FILES(SOURCES_LINALG_PETSC LinAlg/PETSc)
    set(SOURCES ${SOURCES} ${SOURCES_LINALG_PETSC})
endif()

if(OGS_USE_MKL)
    GET_SOURCE_FILES(SOURCES_LINALG_PARDISO LinAlg/Pardiso)
    set(SOURCES ${SOURCES} ${SOURCES_LINALG_PARDISO})
endif()

if(METIS_FOUND)
	GET_SOURCE_FILES(SOURCES_LINALG_SPARSE_NESTEDDISSECTION LinAlg/Sparse/NestedDissectionPermutation)
	set(SOURCES ${SOURCES} ${SOURCES_LINALG_SPARSE_NESTEDDISSECTION})
endif ()

GET_SOURCE_FILES(SOURCES_NONLINEAR Nonlinear)
set(SOURCES ${SOURCES} ${SOURCES_NONLINEAR})

if(METIS_FOUND)
	include_directories(${METIS_INCLUDE_DIR})
endif()

# Create the library
add_library(MathLib STATIC ${SOURCES})

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(MathLib INTERFACE
	GeoLib
	logog
    sundials_cvode
    sundials_nvecserial
)

if(LAPACK_FOUND)
    target_link_libraries(MathLib INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

if (OGS_USE_LIS)
    target_link_libraries(MathLib INTERFACE ${LIS_LIBRARIES})
endif()

if (OGS_USE_PETSC)
    target_link_libraries(MathLib INTERFACE ${PETSC_LIBRARIES})
endif()


#if(OGS_USE_EIGEN)
#	if(OGS_USE_MKL)
#		find_library(INTEL_MKL_CORE_LIBRARIES NAMES mkl_core PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#		find_library(INTEL_MKL_LIBRARIES NAMES mkl_intel_lp64 PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#		message(STATUS "mkl libs: ${INTEL_MKL_CORE_LIBRARIES} ${INTEL_MKL_LIBRARIES}")
#		# target_link_libraries(ApplicationsLib ${INTEL_MKL_LIBRARIES} PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#	endif()
#endif()

if (OGS_USE_MKL)
	LINK_DIRECTORIES(/home/lehmannc/lib/MKL/mkl/lib/intel64)

	set(MKL_LIBRARIES
		/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_core.so
		/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_gnu_thread.so
		/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_intel_lp64.so
	)

	include_directories(SYSTEM /home/lehmannc/lib/MKL/mkl/include)
    target_link_libraries(MathLib INTERFACE ${MKL_LIBRARIES})
endif()


if(TARGET Boost)
    add_dependencies(MathLib Boost)
endif()
if(TARGET Eigen)
	add_dependencies(MathLib Eigen)
endif()

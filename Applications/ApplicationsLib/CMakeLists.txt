# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib STATIC ${LIB_SOURCES})

#if(OGS_USE_EIGEN)
#	if(OGS_USE_MKL)
#		find_library(INTEL_MKL_CORE_LIBRARIES NAMES mkl_core PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#		find_library(INTEL_MKL_LIBRARIES NAMES mkl_intel_lp64 PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#		message(STATUS "mkl libs: ${INTEL_MKL_CORE_LIBRARIES} ${INTEL_MKL_LIBRARIES}")
#		# target_link_libraries(ApplicationsLib ${INTEL_MKL_LIBRARIES} PATHS /home/lehmannc/lib/MKL/mkl/lib/intel64)
#	endif()
#endif()

LINK_DIRECTORIES(/home/lehmannc/lib/MKL/mkl/lib/intel64)

target_link_libraries(ApplicationsLib INTERFACE
	GeoLib
	MeshGeoToolsLib
	NumLib
	ProcessLib
	${VTK_LIBRARIES}
	logog
	/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_core.so
	/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_gnu_thread.so
	/home/lehmannc/lib/MKL/mkl/lib/intel64/libmkl_intel_lp64.so
)
ADD_VTK_DEPENDENCY(ApplicationsLib)

if(TARGET Eigen)
	add_dependencies(ApplicationsLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(ApplicationsLib Boost)
endif()

project(mergeSTP)

add_subdirectory( utility )

# ------------------ Boilerplate

cmake_minimum_required(VERSION 2.4.5)
if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# ----- OpenCASCADE
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(OpenCASCADE REQUIRED TKPrim TKSTL TKSTEP TKernel TKMath
    TKG2d TKShHealing TKBRep TKTopAlgo)
if (OpenCASCADE_FOUND)
    message("Found OpenCASCADE")
    message("Include dirs: ${OpenCASCADE_INCLUDE_DIR}")
else(OpenCASCADE_FOUND)
    message("Did not find OpenCASCADE")
endif(OpenCASCADE_FOUND)


# ---- done with the sillies

find_package(CGAL QUIET COMPONENTS Core )
if (NOT CGAL_FOUND)
    message(FATAL_ERROR "Cannot find CGAL.")
endif (NOT CGAL_FOUND)

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

include_directories(
    ${OpenCASCADE_INCLUDE_DIR}
#    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_executable( mergeSTP
    main.cpp
    InteriorVolumes-inl.h
    PointFacetDistance-inl.h
    ShellMapper.h
    CGALUtilities.h
)

target_link_libraries( mergeSTP
    ${CGAL_LIBRARIES}
    TKPrim
    TKSTL
    TKSTEP
    TKernel
    TKMath
    TKG2d
    TKShHealing
    TKBRep
    TKTopAlgo
    dl # what the heck is this anyway???  didn't need it on other machine.
)


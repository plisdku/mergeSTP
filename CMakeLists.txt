project(mergeSTP)

add_subdirectory( utility )

# ------------------ Boilerplate

cmake_minimum_required(VERSION 2.4.5)
if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

# ----- OCE

find_package(OCE)
if(OCE_FOUND)
    message(STATUS "Found OCE version ${OCE_VERSION}")
    # Include files reside in ${OCE_INCLUDE_DIRS};
    include_directories(${OCE_INCLUDE_DIRS})
    # We do not need library path, they will be automatically imported.
else(OCE_FOUND)
    # OCE not found; either it is not found and user
    # has to set OCE_DIR to the directory containing
    # OCEConfig.cmake, or OCE is not installed and we
    # try to find OpenCascade files.
    if(DEFINED ENV{CASROOT})
        if(NOT DEFINED OCC_INCLUDE_PATH)
            set(OCC_INCLUDE_PATH "$ENV{CASROOT}/inc")
        endif(NOT DEFINED OCC_INCLUDE_PATH)
        if(NOT DEFINED OCC_LIB_PATH)
            if(WIN32)
                set(OCC_LIB_PATH "$ENV{CASROOT}/win32/lib")
            else(WIN32)
                if(APPLE)
                    set(OCC_LIB_PATH "/Library/OpenCASCADE/6.3.0/lib")
                else(APPLE)
                    set(OCC_LIB_PATH "$ENV{CASROOT}/lin/lib")
                endif(APPLE)
            endif(WIN32)
        endif(NOT DEFINED OCC_LIB_PATH)
    else(DEFINED ENV{CASROOT})
        if(NOT DEFINED OCC_INCLUDE_PATH OR NOT DEFINED OCC_LIB_PATH)
            message(WARNING "To specify paths of OpenCascade files, you may either define the CASROOT environment variable, or set both OCC_INCLUDE_PATH and OCC_LIB_PATH variables.")
        endif(NOT DEFINED OCC_INCLUDE_PATH OR NOT DEFINED OCC_LIB_PATH)
    endif(DEFINED ENV{CASROOT})
    if(DEFINED OCC_INCLUDE_PATH)
        message(STATUS "OCC search path for include files: OCC_INCLUDE_PATH=${OCC_INCLUDE_PATH}")
        include_directories(${OCC_INCLUDE_PATH})
    endif(DEFINED OCC_INCLUDE_PATH)
    if(DEFINED OCC_LIB_PATH)
        message(STATUS "OCC search path for libraries: OCC_LIB_PATH=${OCC_LIB_PATH}")
        link_directories(${OCC_LIB_PATH})
    endif(DEFINED OCC_LIB_PATH)
endif(OCE_FOUND)

# ---- done with the sillies

find_package(CGAL QUIET COMPONENTS Core )
if (NOT CGAL_FOUND)
    message(FATAL_ERROR "Cannot find CGAL.")
endif (NOT CGAL_FOUND)

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

include_directories(
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
)


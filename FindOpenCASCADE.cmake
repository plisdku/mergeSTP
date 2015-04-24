# - Try to find OpenCASCADE libraries
### Does NOT test what version has been found,though
### that could be done by parsing Standard_Version.hxx

message("Paul says, CASROOT is $ENV{CASROOT}.")
message("Paul says OpenCASCADE_FOUND = ${OpenCASCADE_FOUND}.")

# Once done, this will define
#  OpenCASCADE_FOUND - true if OCC has been found
#  OpenCASCADE_INCLUDE_DIR - the OCC include dir
#  OpenCASCADE_LIBRARIES - names of OCC libraries
#  OpenCASCADE_LINK_DIRECTORY - location of OCC libraries

# ${OpenCASCADE_FOUND} is cached, so once OCC is found this block shouldn't have to run again
if( NOT OpenCASCADE_FOUND STREQUAL TRUE )
    message(STATUS "Searching for OpenCASCADE.")
    if(UNIX)
        set( _incsearchpath /usr/include/opencascade /opt/occ/inc $ENV{CASROOT}/inc )
        if (APPLE)
            set( _testlibname libTKernel.dylib )
        else (APPLE)
            set( _testlibname libTKernel.so )
        endif (APPLE)
        set( _libsearchpath /usr/lib /opt/occ/lib $ENV{CASROOT}/lib )
    else(UNIX)
        if (WIN32)
            message("************ FindOpenCASCADE.cmake has NOT been tried on windows and may or may NOT work! *************")
            set( _incsearchpath $ENV{CASROOT}\\inc C:\\OpenCASCADE6.3.0\\ros\\inc )
            set( _testlibname TKernel.dll )
            set( _libsearchpath $ENV{CASROOT}\\win32\\bin C:\\OpenCASCADE6.3.0\\ros\\win32\\bin )
        else(WIN32)
            message( FATAL_ERROR "Unknown system! Exiting." )
        endif (WIN32)
    endif (UNIX)

    message("Paul says, finding the include directory by looking for Standard_Real.hxx")
    #find the include dir by looking for Standard_Real.hxx

    find_path( OpenCASCADE_INCLUDE_DIR Standard_Real.hxx PATHS ${_incsearchpath} DOC "Path to OCC includes" )
    if( OpenCASCADE_INCLUDE_DIR STREQUAL Standard_Real.hxx-NOTFOUND )
        set( OpenCASCADE_FOUND FALSE CACHE BOOL FORCE )
        message( FATAL_ERROR "Cannot find OCC include dir. Install opencascade or set CASROOT or create a symlink /opt/occ/inc pointing to the correct directory." )
    else(OpenCASCADE_INCLUDE_DIR STREQUAL Standard_Real.hxx-NOTFOUND)
        message("Found the include directory.")
    endif( OpenCASCADE_INCLUDE_DIR STREQUAL Standard_Real.hxx-NOTFOUND )

    # Find one lib and save its directory to OpenCASCADE_LINK_DIRECTORY. Because
    #  OCC has so many libs, there is increased risk of a name collision.
    #  Requiring that all libs be in the same directory reduces the risk.
    find_path( OpenCASCADE_LINK_DIRECTORY ${_testlibname} PATH ${_libsearchpath} DOC "Path to OCC libs" )
    if( OpenCASCADE_LINK_DIRECTORY STREQUAL ${_testlibname}-NOTFOUND )
        set( OpenCASCADE_FOUND FALSE CACHE BOOL FORCE )
        message( FATAL_ERROR "Cannot find OCC lib dir. Install opencascade or set CASROOT or create a symlink /opt/occ/lib pointing to the dir where the OCC libs are." )
    else( OpenCASCADE_LINK_DIRECTORY STREQUAL ${_testlibname}-NOTFOUND )
        set( OpenCASCADE_FOUND TRUE CACHE BOOL "Has OCC been found?" FORCE )
        set( _firsttime TRUE ) #so that messages are only printed once
        message( STATUS "Found OCC include dir: ${OpenCASCADE_INCLUDE_DIR}" )
        message( STATUS "Found OCC lib dir: ${OpenCASCADE_LINK_DIRECTORY}" )
    endif( OpenCASCADE_LINK_DIRECTORY STREQUAL ${_testlibname}-NOTFOUND )
else( NOT OpenCASCADE_FOUND STREQUAL TRUE )
    message(STATUS "Not searching for OpenCASCADE, because we found it.")
    set( _firsttime FALSE ) #so that messages are only printed once
endif( NOT OpenCASCADE_FOUND STREQUAL TRUE )

message( "Include dir is ${OpenCASCADE_INCLUDE_DIR}")
message( "Link dir is ${OpenCASCADE_LINK_DIRECTORY}")
message( "_libsearchpath is ${_libsearchpath}" )


if( OpenCASCADE_FOUND STREQUAL TRUE )
    if( DEFINED OpenCASCADE_FIND_COMPONENTS )
        foreach( _libname ${OpenCASCADE_FIND_COMPONENTS} )
        
            #look for libs in OpenCASCADE_LINK_DIRECTORY
            find_library( ${_libname}_OCCLIB ${_libname} ${OpenCASCADE_LINK_DIRECTORY} NO_DEFAULT_PATH)
            set( _foundlib ${${_libname}_OCCLIB} )
            
            if( _foundlib STREQUAL ${_libname}_OCCLIB-NOTFOUND )
                message( FATAL_ERROR "Cannot find ${_libname}. Is it spelled correctly? Correct capitalization? Do you have another package with similarly-named libraries, installed at ${OpenCASCADE_LINK_DIRECTORY}? (That is where this script thinks the OCC libs are.)" )
                message( "Link dir is ${OpenCASCADE_LINK_DIRECTORY}" )
            endif( _foundlib STREQUAL ${_libname}_OCCLIB-NOTFOUND )
            set( OpenCASCADE_LIBRARIES ${OpenCASCADE_LIBRARIES} ${_foundlib} )
        endforeach( _libname ${OpenCASCADE_FIND_COMPONENTS} )

        if (UNIX)
            ADD_DEFINITIONS( -DLIN -DLININTEL )
        endif (UNIX)

        # 32 bit or 64 bit?
        if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
            if( _firsttime STREQUAL TRUE )
                message( STATUS "This is a 32-bit system." )
            endif( _firsttime STREQUAL TRUE )
        else( CMAKE_SIZEOF_VOID_P EQUAL 4 )
            if( _firsttime STREQUAL TRUE )
                message( STATUS "This is a 64-bit system. Adding appropriate compiler flags for OCC." )
            endif( _firsttime STREQUAL TRUE )
            add_definitions( -D_OCC64 )
            
            if (UNIX)
                add_definitions( -m64 )
            endif (UNIX)
        endif( CMAKE_SIZEOF_VOID_P EQUAL 4 )

        add_definitions( -DHAVE_CONFIG_H -DHAVE_IOSTREAM -DHAVE_FSTREAM -DHAVE_LIMITS_H -DHAVE_IOMANIP )
    else( DEFINED OpenCASCADE_FIND_COMPONENTS )
        message( AUTHOR_WARNING "Developer must specify required libraries to link against in the cmake file, i.e. find_package( OpenCASCADE REQUIRED COMPONENTS TKernel TKBRep) . Otherwise no libs will be added - linking against ALL OCC libraries is slow!")
    endif( DEFINED OpenCASCADE_FIND_COMPONENTS )
endif( OpenCASCADE_FOUND STREQUAL TRUE )


set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")

find_program(R_EXECUTABLE NAMES R)

# ---searching R installtion unsing R executable
if(R_EXECUTABLE)
    execute_process(
        COMMAND ${R_EXECUTABLE} --slave -e "cat(strsplit(version[['version.string']], ' ')[[1]][3])"
        OUTPUT_VARIABLE R_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    execute_process(COMMAND ${R_EXECUTABLE} RHOME
        OUTPUT_VARIABLE R_ROOT_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    find_path(R_INCLUDE_DIR R.h
        HINTS ${R_ROOT_DIR}
        PATHS /usr/local/lib /usr/local/lib64 /usr/share
        PATH_SUFFIXES include R/include
        DOC "Path to file R.h")

    set(CMAKE_FIND_LIBRARY_PREFIXES lib)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
    find_library(R_LIBRARY
        NAMES R
        HINTS ${R_ROOT_DIR}/lib HINTS ${R_ROOT_DIR}/libs
        DOC "R library (example libR.a, libR.dylib, etc.).")
endif()

# ---setting include dirs and libraries
set(R_LIBRARIES ${R_LIBRARY})
set(R_INCLUDE_DIRS ${R_INCLUDE_DIR})

foreach(_cpt ${R_FIND_COMPONENTS})
    execute_process(COMMAND echo "cat(find.package('${_cpt}'))"
        COMMAND ${R_EXECUTABLE} --vanilla --slave
        RESULT_VARIABLE _rc
        ERROR_QUIET
        OUTPUT_VARIABLE _cpt_path
        OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT _rc)
        set(R_${_cpt}_FOUND 1)
    endif()

    find_library(R_${_cpt}_LIBRARY
        lib${_cpt}.so lib${_cpt}.dylib
        HINTS ${_cpt_path}/libs ${_cpt_path}/lib)

    if(R_${_cpt}_LIBRARY)
        mark_as_advanced(R_${_cpt}_LIBRARY)
        list(APPEND R_LIBRARIES ${R_${_cpt}_LIBRARY})
    endif()

    find_path(R_${_cpt}_INCLUDE_DIR ${_cpt}.h HINTS ${_cpt_path} PATH_SUFFIXES include R/include)

    if(R_${_cpt}_INCLUDE_DIR)
        mark_as_advanced(R_${_cpt}_INCLUDE_DIR)
        list(APPEND R_INCLUDE_DIRS ${R_${_cpt}_INCLUDE_DIR})
    endif()
endforeach()

# Handle the QUIETLY and REQUIRED arguments and set R_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(R HANDLE_COMPONENTS
    REQUIRED_VARS R_EXECUTABLE R_INCLUDE_DIR R_LIBRARY
    VERSION_VAR R_VERSION)
mark_as_advanced(R_FOUND R_EXECUTABLE R_INCLUDE_DIR R_LIBRARY)

if(R_FOUND AND NOT TARGET R::R)
    add_library(R::R INTERFACE IMPORTED)
    set_target_properties(R::R PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${R_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${R_LIBRARIES}")
endif()

# - Try to find BGL-Extra
# Once done this will define
#
#  BGLEXTRA_FOUND - system has BGL-Extra
#  BGLEXTRA_INCLUDE_DIR - the BGL-Extra include directory
#  BGLEXTRA_LIBRARIES - Link these to use BGL-Extra (empty)
#  BGLEXTRA_VERSION_STRING - the version of BGL-Extra found (since CMake 2.8.8)

#=============================================================================
# Copyright 2013 Sven Mikael Persson <mikael.s.persson@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# set(_BGLEXTRA_PATHS 
# PATHS
#   ""
# )

# This is a header-only library:
set(BGLEXTRA_LIBRARIES "")
set(BGLEXTRA_LIBRARY ${BGLEXTRA_LIBRARIES})

# Currently, only a single version (pending proposal to Boost anyways):
set(BGLEXTRA_VERSION_STRING "1.0")

if(WIN32)
  set(_BGLEXTRA_PATH_HINTS "C:/" "C:/Boost/")
else()
  # Add hints for unix-like systems.. not really needed.
  set(_BGLEXTRA_PATH_HINTS "/usr/include" "/usr/local/include")
endif()

find_path(_BGLEXTRA_INCLUDE_DIR adjacency_list_BC.hpp HINTS ${_BGLEXTRA_PATH_HINTS} PATH_SUFFIXES "libbgl-extra/boost/graph" "boost/graph")
string(REPLACE "/boost/graph" "" BGLEXTRA_INCLUDE_DIR ${_BGLEXTRA_INCLUDE_DIR})
set(BGLEXTRA_INCLUDE_DIRS ${BGLEXTRA_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set BGLEXTRA_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BGLEXTRA
                                  REQUIRED_VARS BGLEXTRA_INCLUDE_DIR
                                  VERSION_VAR BGLEXTRA_VERSION_STRING)

# if (BGLEXTRA_FOUND)
#    include(${CMAKE_CURRENT_LIST_DIR}/CheckLibraryExists.cmake)
#    CHECK_LIBRARY_EXISTS("${BGLEXTRA_LIBRARIES}" BZ2_bzCompressInit "" BZIP2_NEED_PREFIX)
# endif ()

mark_as_advanced(BGLEXTRA_INCLUDE_DIR)
mark_as_advanced(BGLEXTRA_INCLUDE_DIRS)
mark_as_advanced(BGLEXTRA_LIBRARIES)
mark_as_advanced(BGLEXTRA_LIBRARY)

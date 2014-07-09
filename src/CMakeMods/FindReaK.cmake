# - Try to find ReaK Libraries
# Once done this will define
#
#  ReaK_FOUND - system has ReaK
#  ReaK_INCLUDE_DIR - the ReaK include directory
#  ReaK_INCLUDE_DIRS - the ReaK include directories (same as ReaK_INCLUDE_DIR)
#  ReaK_LIBRARY - Link these to use ReaK
#  ReaK_LIBRARIES - Link these to use ReaK (same as ReaK_LIBRARY)
#
# And for all the components requested in the find_package command, there will 
# be the following variables:
#  
#  ReaK_component_FOUND - system has ReaK.component
#  ReaK_component_INCLUDE_DIR - the ReaK.component include directory
#  ReaK_component_INCLUDE_DIRS - the ReaK.component include directories (same as ReaK_component_INCLUDE_DIR)
#  ReaK_component_LIBRARY - Link these to use ReaK.component
#  ReaK_component_LIBRARIES - Link these to use ReaK.component (same as ReaK_component_LIBRARY)
#

#=============================================================================
# Copyright 2014 Sven Mikael Persson <mikael.s.persson@gmail.com>
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

# This is a header-only library:
set(ReaK_LIBRARIES "")

# Currently, only a single version (pending proposal to Boost anyways):
set(ReaK_VERSION_STRING "0.28")

find_path(_ReaK_INCLUDE_DIR reak_core_version.hpp PATH_SUFFIXES "ReaK/core")
string(REPLACE "/ReaK/core" "" ReaK_INCLUDE_DIR ${_ReaK_INCLUDE_DIR})

# Extract ReaK_VERSION from reak_core_version.hpp
set(ReaK_VERSION 0)
file(STRINGS "${_ReaK_INCLUDE_DIR}/reak_core_version.hpp" _ReaK_VERSION_HPP_CONTENTS REGEX "#define REAK_VERSION ")
set(_ReaK_VERSION_REGEX "([0-9]+)")
if("${_ReaK_VERSION_HPP_CONTENTS}" MATCHES ".*#define REAK_VERSION ${_ReaK_VERSION_REGEX}.*")
  set(ReaK_VERSION "${CMAKE_MATCH_1}")
endif()
unset(_ReaK_VERSION_HPP_CONTENTS)
unset(_ReaK_INCLUDE_DIR)

math(EXPR ReaK_MAJOR_VERSION "${ReaK_VERSION} / 100000")
math(EXPR ReaK_MINOR_VERSION "${ReaK_VERSION} / 100 % 1000")
math(EXPR ReaK_SUBMINOR_VERSION "${ReaK_VERSION} % 100")

# handle the QUIETLY and REQUIRED arguments and set REAK_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ReaK 
  FOUND_VAR ReaK_FOUND
  REQUIRED_VARS ReaK_INCLUDE_DIR
  VERSION_VAR "${ReaK_MAJOR_VERSION}.${ReaK_MINOR_VERSION}-${ReaK_SUBMINOR_VERSION}")

set(ReaK_INCLUDE_DIRS ${ReaK_INCLUDE_DIR})

foreach(COMPONENT ${ReaK_FIND_COMPONENTS})
  if(${COMPONENT} STREQUAL "core" OR
     ${COMPONENT} STREQUAL "math")
    set(COMPONENT_BASE_REL_PATH "ReaK/core")
  elseif(${COMPONENT} STREQUAL "mbd" OR
         ${COMPONENT} STREQUAL "ctrl_est" OR
         ${COMPONENT} STREQUAL "motion_plan" OR
         ${COMPONENT} STREQUAL "geometry")
    set(COMPONENT_BASE_REL_PATH "ReaK/ctrl")
  endif()
  
  find_path(_ReaK_${COMPONENT}_INCLUDE_DIR reak_${COMPONENT}_version.hpp PATH_SUFFIXES ${COMPONENT_BASE_REL_PATH})
  string(REPLACE ${COMPONENT_BASE_REL_PATH} "" ReaK_${COMPONENT}_INCLUDE_DIR ${_ReaK_${COMPONENT}_INCLUDE_DIR})
  
  find_library(ReaK_${COMPONENT}_LIBRARY reak_${COMPONENT})
  if(${ReaK_${COMPONENT}_LIBRARY} STREQUAL "ReaK_${COMPONENT}_LIBRARY-NOTFOUND")
    unset(ReaK_${COMPONENT}_LIBRARY)
  endif()
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(ReaK_${COMPONENT} 
    FOUND_VAR ReaK_${COMPONENT}_FOUND
    REQUIRED_VARS ReaK_${COMPONENT}_INCLUDE_DIR ReaK_${COMPONENT}_LIBRARY
    VERSION_VAR "${ReaK_MAJOR_VERSION}.${ReaK_MINOR_VERSION}-${ReaK_SUBMINOR_VERSION}")
  
  set(ReaK_LIBRARIES ${ReaK_LIBRARIES} ${ReaK_${COMPONENT}_LIBRARY})
  
  set(ReaK_${COMPONENT}_WITH_UNIQUE_INCL TRUE)
  foreach(INCL_DIR ${ReaK_INCLUDE_DIRS})
    if(NOT ${ReaK_${COMPONENT}_WITH_UNIQUE_INCL} AND ${INCL_DIR} STREQUAL ${ReaK_${COMPONENT}_INCLUDE_DIR})
      set(ReaK_${COMPONENT}_WITH_UNIQUE_INCL FALSE)
    endif()
  endforeach()
  if(${ReaK_${COMPONENT}_WITH_UNIQUE_INCL})
    set(ReaK_INCLUDE_DIRS ${ReaK_INCLUDE_DIRS} ${ReaK_${COMPONENT}_INCLUDE_DIR})
  endif()
  unset(ReaK_${COMPONENT}_WITH_UNIQUE_INCL)
  
  mark_as_advanced(ReaK_${COMPONENT}_INCLUDE_DIR)
  mark_as_advanced(ReaK_${COMPONENT}_LIBRARY)
endforeach()

mark_as_advanced(ReaK_INCLUDE_DIR)
mark_as_advanced(ReaK_INCLUDE_DIRS)
mark_as_advanced(ReaK_LIBRARIES)
set(ReaK_LIBRARY ${ReaK_LIBRARIES})
mark_as_advanced(ReaK_LIBRARY)

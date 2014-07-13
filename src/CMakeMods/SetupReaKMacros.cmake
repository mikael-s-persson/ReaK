# This cmake module sets up compiler flags, doxygen macros, and custom target macros
# for use when building parts of the ReaK library.
# 
# 
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


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Configured for build-type: ${CMAKE_BUILD_TYPE}")


enable_testing()

include(CheckCXXCompilerFlag)

# Check and enable C++11 or C++0x features:
CHECK_CXX_COMPILER_FLAG("-std=c++11" _COMPILER_SUPPORTS_CXX11)
if(_COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  message(STATUS "Configured compiler ${CMAKE_CXX_COMPILER} to use C++11 support.")
else()
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" _COMPILER_SUPPORTS_CXX0X)
  if(_COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    message(STATUS "Configured compiler ${CMAKE_CXX_COMPILER} to use C++0x support. Consider using a newer compiler for full C++11 support.")
  else()
    if(NOT WIN32 AND NOT MSVC)
      message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Consider using a newer compiler for optimal performance.")
    endif()
  endif()
endif()

if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox")
    message(STATUS "Configured compiler options and output directories for MSVC toolset.")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -O3 -Wl,--no-as-needed -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} --enable-stdcall-fixup")
    message(STATUS "Configured compiler options and output directories for MinGW toolset.")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -pthread -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    # TODO This is just a temporary hack because of a version of clang failing to compile with a experimental 4.9 version of libstdc++.
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -nostdinc++ -isystem /usr/include/c++/4.8 -isystem /usr/include/x86_64-linux-gnu/c++/4.8")
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  else()
#     set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  endif()
  message(STATUS "Configured compiler options and output directories for *nix GCC toolset.")
endif()


include(ConfigDoxygenTargets)


# Now set the global macros for setting up targets.
macro(setup_custom_target target_name)
  message(STATUS "Registered target ${target_name}.")
endmacro(setup_custom_target)

macro(setup_custom_test_program target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests")
  add_test(NAME "${target_name}" WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests/" COMMAND "$<TARGET_FILE:${target_name}>")
endmacro(setup_custom_test_program)

macro(setup_headers HEADER_FILES HEADER_PATH)
  foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
    install(FILES "${SRCINCLUDEROOT}${CURRENT_HEADER_FILE}" DESTINATION "include${HEADER_PATH}")
  endforeach(CURRENT_HEADER_FILE)
endmacro(setup_headers)



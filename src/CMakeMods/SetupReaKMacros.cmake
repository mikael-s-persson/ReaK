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
if(MSVC)
  if(MSVC_VERSION GREATER 1700)
    message(STATUS "This MSVC compiler version (>1700) has C++11 support enabled by default.")
  elseif(MSVC_VERSION GREATER 1600)
    message(STATUS "This MSVC compiler version (1700) has limited C++11 support. Consider using a newer version with better C++11 support, such as Visual Studio 2013 and above.")
  else()
    message(STATUS "This MSVC compiler version (<1700) has no C++11 support. Consider using a newer compiler for optimal performance, such as Visual Studio 2013 and above.")
  endif()
else()
  CHECK_CXX_COMPILER_FLAG("-std=c++11" _COMPILER_SUPPORTS_CXX11)
  if(_COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    message(STATUS "The compiler has C++11 support. C++11 mode was enabled.")
  else()
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" _COMPILER_SUPPORTS_CXX0X)
    if(_COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
      message(STATUS "The compiler only has experimental C++11 support. Consider using a newer compiler with full C++11 support.")
    else()
      message(STATUS "The compiler has no detectable C++11 support. Consider using a newer compiler for optimal performance.")
    endif()
    mark_as_advanced(_COMPILER_SUPPORTS_CXX0X)
  endif()
  mark_as_advanced(_COMPILER_SUPPORTS_CXX11)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
#   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -ftemplate-depth=2000 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
#   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  if (WIN32)
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} --enable-stdcall-fixup")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    # This is a hack because the major-minor version variables (${CLANG_VERSION_MAJOR}.${CLANG_VERSION_MINOR}) from cmake seem unreliable 
    execute_process( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE clang_full_version_string )
    string (REGEX REPLACE ".*clang version ([0-9]+\\.[0-9]+).*" "\\1" CLANG_VERSION_STRING ${clang_full_version_string})
    
    # TODO This is just a temporary hack because of a version of clang failing to compile with a experimental 4.9 version of libstdc++.
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -nostdinc++ -isystem /usr/include/c++/4.8 -isystem /usr/include/x86_64-linux-gnu/c++/4.8")
    
    message(STATUS "Detected Clang version: ${CLANG_VERSION_STRING}")
    if( CLANG_VERSION_STRING VERSION_LESS 3.6 )
#     if( "${CLANG_VERSION_MAJOR}.${CLANG_VERSION_MINOR}" STRLESS "3.6") 
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef -Wno-missing-braces")
    endif()
  endif()
endif()
message(STATUS "Configured compiler options for ${CMAKE_SYSTEM_NAME} system with ${CMAKE_CXX_COMPILER_ID} toolset.")
message(STATUS "   using C++ options: ${CMAKE_CXX_FLAGS}")


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



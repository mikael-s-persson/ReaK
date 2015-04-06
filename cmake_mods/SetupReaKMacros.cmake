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
    set(_COMPILER_SUPPORTS_CXX11 TRUE)
    set(_COMPILER_SUPPORTS_CXX0X TRUE)
  elseif(MSVC_VERSION GREATER 1600)
    message(STATUS "This MSVC compiler version (1700) has limited C++11 support. Consider using a newer version with better C++11 support, such as Visual Studio 2013 and above.")
    set(_COMPILER_SUPPORTS_CXX11 FALSE)
    set(_COMPILER_SUPPORTS_CXX0X TRUE)
  else()
    message(FATAL_ERROR "This MSVC compiler version (<1700) has no C++11 support. Use a newer compiler, such as Visual Studio 2013 and above.")
    set(_COMPILER_SUPPORTS_CXX11 FALSE)
    set(_COMPILER_SUPPORTS_CXX0X FALSE)
  endif()
else()
  CHECK_CXX_COMPILER_FLAG("-std=c++11" _COMPILER_SUPPORTS_CXX11)
  if(_COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    message(STATUS "The compiler has C++11 support. C++11 mode was enabled.")
    set(_COMPILER_SUPPORTS_CXX0X TRUE)
  else()
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" _COMPILER_SUPPORTS_CXX0X)
    if(_COMPILER_SUPPORTS_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
      message(STATUS "The compiler only has experimental C++11 support. Consider using a newer compiler with full C++11 support.")
    else()
      message(FATAL_ERROR "The compiler has no detectable C++11 support. Use a newer compiler.")
      set(_COMPILER_SUPPORTS_CXX0X FALSE)
    endif()
  endif()
endif()
mark_as_advanced(_COMPILER_SUPPORTS_CXX0X)
mark_as_advanced(_COMPILER_SUPPORTS_CXX11)


if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /bigobj /wd4996 /D_SCL_SECURE_NO_WARNINGS /DBOOST_NO_CXX11_EXTERN_TEMPLATE")
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


# Look for winsock2:

if( WIN32 )
# #  find_library( WINSOCK_LIB wsock32 PATHS "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib")
#   find_library( WINSOCK_LIB wsock32 )
#   find_library( WS2_32_LIB ws2_32 )7
  set( WINSOCK_LIB wsock32 )
  set( WS2_32_LIB ws2_32 )
  if( NOT WINSOCK_LIB AND NOT WS2_32_LIB )
    message( FATAL_ERROR  "Could not located the winsock library (ws2_32)")
  else()
    if( WS2_32_LIB )
      set(EXTRA_SYSTEM_LIBS ${EXTRA_SYSTEM_LIBS} ${WS2_32_LIB})
    endif()
    if( WINSOCK_LIB )
      set(EXTRA_SYSTEM_LIBS ${EXTRA_SYSTEM_LIBS} ${WINSOCK_LIB})
    endif()
  endif()
endif()

# Look for Boost:

# set(Boost_DEBUG TRUE)

set(Boost_ADDITIONAL_VERSIONS "1.49" "1.49.0" "1.50" "1.50.0" "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0" "1.54" "1.54.0" "1.55" "1.55.0" "1.56" "1.56.0" "1.57" "1.57.0")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

# if a custom path for boost is provided, than use that (and suppress system paths).
if(EXISTS "${CUSTOM_BOOST_PATH}/include/boost")
  set(Boost_INCLUDE_DIR "${CUSTOM_BOOST_PATH}/include")
  set(Boost_LIBRARY_DIR "${CUSTOM_BOOST_PATH}/lib")
  set(Boost_NO_SYSTEM_PATHS TRUE)
endif()

if (NOT WIN32)
  # make sure that the *nix suffixes and prefixes are correct (some cmake installs of findBoost.cmake are wrong with this).
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES})
  if( Boost_USE_STATIC_LIBS )
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
find_package(Boost 1.49)
else()
find_package(Boost 1.49 COMPONENTS thread chrono date_time system program_options unit_test_framework filesystem random REQUIRED)
endif()
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  add_definitions( "-DREAK_HAS_BOOST" )
  if(MSVC)
    # Disable the libraries, since it uses automatic linking:
    set(Boost_LIBRARIES "")
    find_path(Boost_LIBRARY_DIRS DEPENDENCY_VERSIONS.txt HINTS ${Boost_INCLUDE_DIR})
  endif()
  link_directories(${Boost_LIBRARY_DIRS})
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
endif()

if( NOT WIN32 )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_PREFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES})
endif()




macro(get_subdir_list result curdir)
  file(GLOB children RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/${curdir}" "${CMAKE_CURRENT_SOURCE_DIR}/${curdir}/*")
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${curdir}/${child})
      if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${curdir}/${child}/CMakeLists.txt")
        list(APPEND dirlist ${child})
      endif()
    endif()
  endforeach(child)
  set(${result} ${dirlist})
endmacro(get_subdir_list)

macro(add_subdirectory_if_cmake subdir)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/${subdir}")
  endif()
endmacro(add_subdirectory_if_cmake)

include(ConfigDoxygenTargets)


# Now set the global macros for setting up targets.
macro(ReaK_setup_tool_program target_name)
  install(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT ReaK_${REAK_CURRENT_MODULE}_tools)
  message(STATUS "Registered ReaK tool program ${target_name}.")
endmacro(ReaK_setup_tool_program)

macro(ReaK_setup_static_library target_name)
  install(TARGETS ${target_name} ARCHIVE DESTINATION lib COMPONENT ReaK_${REAK_CURRENT_MODULE})
  message(STATUS "Registered ReaK static library ${target_name}.")
endmacro(ReaK_setup_static_library)

macro(ReaK_setup_shared_library target_name)
  install(TARGETS ${target_name} LIBRARY DESTINATION lib COMPONENT ReaK_${REAK_CURRENT_MODULE})
  message(STATUS "Registered ReaK shared library ${target_name}.")
endmacro(ReaK_setup_shared_library)

macro(ReaK_setup_target target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples")
  message(STATUS "Registered ReaK target ${target_name}.")
endmacro(ReaK_setup_target)

macro(ReaK_setup_test_program target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests")
  add_test(NAME "${target_name}" WORKING_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/unit_tests/" COMMAND "$<TARGET_FILE:${target_name}>")
  message(STATUS "Registered ReaK test program ${target_name}.")
endmacro(ReaK_setup_test_program)


macro(ReaK_add_current_module) 
  
  configure_doxyfile(ReaK_${REAK_CURRENT_MODULE} 
                    "\"${REAK_CURRENT_MODULE_TITLE}\"" 
                    "\"${REAK_CURRENT_MODULE_DESCRIPTION}\"" 
                    "include")
  add_doxygen_target(ReaK_${REAK_CURRENT_MODULE})
  
  include_directories(AFTER "include")
  
  add_subdirectory_if_cmake("src")
  add_subdirectory_if_cmake("tools")
  add_subdirectory_if_cmake("test")
  add_subdirectory_if_cmake("perf")
  add_subdirectory_if_cmake("qt_ui")
  add_subdirectory_if_cmake("example")
  
  install(DIRECTORY include/ DESTINATION include COMPONENT ReaK_${REAK_CURRENT_MODULE}_dev FILES_MATCHING PATTERN "*.hpp")
  
  cpack_add_component(ReaK_${REAK_CURRENT_MODULE}
    DISPLAY_NAME "${REAK_CURRENT_MODULE_TITLE}"
    DESCRIPTION "${REAK_CURRENT_MODULE_DESCRIPTION}"
  )
  
  cpack_add_component(ReaK_${REAK_CURRENT_MODULE}_dev
    DISPLAY_NAME "${REAK_CURRENT_MODULE_TITLE} - Development files"
    DESCRIPTION "${REAK_CURRENT_MODULE_DESCRIPTION} Development files, mainly headers."
  )
  
  cpack_add_component(ReaK_${REAK_CURRENT_MODULE}_tools
    DISPLAY_NAME "${REAK_CURRENT_MODULE_TITLE} - Tools"
    DESCRIPTION "${REAK_CURRENT_MODULE_DESCRIPTION} Utility applications associated to this module."
  )
  
  
endmacro()




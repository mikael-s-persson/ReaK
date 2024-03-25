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

set(CMAKE_CXX_STANDARD 20)


enable_testing()

include(CheckCXXCompilerFlag)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -stdlib=libc++ -pthread -ftemplate-depth=2000 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# These flags are required for Boost libraries, which don't work with libc++.
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
# add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)

execute_process( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE clang_full_version_string )
message(STATUS "Detected compiler version: ${clang_full_version_string}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedef -Wno-missing-braces")

# search for clang-tidy
find_program(CLANG_TIDY_EXE NAMES "clang-tidy" REQUIRED)
# setup clang-tidy command from executable + options
set(CLANG_TIDY_COMMAND "${CLANG_TIDY_EXE}" "-checks=-*,cppcoreguidelines-init-variables,google-explicit-constructor,misc-unused-using-decls,performance-noexcept-move-constructor,performance-trivially-destructible,performance-unnecessary-copy-initialization,readability-*,-readability-identifier-length,-readability-magic-numbers,-readability-function-cognitive-complexity,-readability-convert-member-functions-to-static,modernize-*,-modernize-use-trailing-return-type,-modernize-use-nodiscard,-modernize-pass-by-value" "-fix" "-header-filter=include/ReaK")

message(STATUS "Configured compiler options for ${CMAKE_SYSTEM_NAME} system.")
message(STATUS "   using C++ options: ${CMAKE_CXX_FLAGS}")


# Process Abseil's CMake build system
set(ABSL_PROPAGATE_CXX_STD ON)
add_subdirectory(third_party/abseil-cpp)

# Process googletest's CMake build system
add_subdirectory(third_party/googletest)

# Process boost_graph_ext_mp's CMake build system
add_subdirectory(third_party/boost_graph_ext_mp)
include_directories(SYSTEM "third_party/boost_graph_ext_mp/include")

# Look for Python (needed by Boost.Python):

find_package(Python3 COMPONENTS Development)
if(Python3_FOUND)
  include_directories(SYSTEM ${Python3_INCLUDE_DIRS})
  link_directories(${Python3_LIBRARY_DIRS})
  link_libraries(${Python3_LIBRARIES})
else()
  message(WARNING "The ReaK python bindings module requires the presence of the python libraries!")
endif()


# Look for Boost:

set(Boost_DEBUG TRUE)

set(Boost_ADDITIONAL_VERSIONS "1.70" "1.70.0" "1.71" "1.71.0" "1.72" "1.72.0" "1.73" "1.73.0" "1.74" "1.74.0")
set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_MULTITHREADED ON)

find_package(Boost 1.70 COMPONENTS python REQUIRED)

if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
endif()

include(ConfigDoxygenTargets)


# Now set the global macros for setting up targets.

# reak_cc_library()
#
# CMake function to imitate Bazel's cc_library rule.
#
# Parameters:
# NAME: name of target (see Note)
# HDRS: List of public header files for the library
# SRCS: List of source files for the library
# DEPS: List of other libraries to be linked in to the binary targets
# COPTS: List of private compile options
# DEFINES: List of public defines
# LINKOPTS: List of link options
# PUBLIC: Add this so that this library will be exported under ReaK::
#
# Note:
# By default, reak_cc_library will always create a library named reak_${NAME},
# and alias target ReaK::${NAME}.  The ReaK:: form should always be used.
# This is to reduce namespace pollution.
#
# reak_cc_library(
#   NAME
#     awesome
#   HDRS
#     "a.h"
#   SRCS
#     "a.cc"
# )
# reak_cc_library(
#   NAME
#     fantastic_lib
#   SRCS
#     "b.cc"
#   DEPS
#     absl::awesome # not "awesome" !
#   PUBLIC
# )
#
# reak_cc_library(
#   NAME
#     main_lib
#   ...
#   DEPS
#     absl::fantastic_lib
# )
function(reak_cc_library)
  cmake_parse_arguments(REAK_CC_LIB
    "DISABLE_INSTALL;PUBLIC"
    "NAME"
    "HDRS;SRCS;COPTS;DEFINES;LINKOPTS;DEPS"
    ${ARGN}
  )

  set(_NAME "reak_${REAK_CC_LIB_NAME}")

  # Check if this is a header-only library
  # Note that as of February 2019, many popular OS's (for example, Ubuntu
  # 16.04 LTS) only come with cmake 3.5 by default.  For this reason, we can't
  # use list(FILTER...)
  set(REAK_CC_SRCS "${REAK_CC_LIB_SRCS}")
  foreach(src_file IN LISTS REAK_CC_SRCS)
    if(${src_file} MATCHES ".*\\.(h|inc|hpp|tpp)")
      list(REMOVE_ITEM REAK_CC_SRCS "${src_file}")
    endif()
  endforeach()

  if(REAK_CC_SRCS STREQUAL "")
    set(REAK_CC_LIB_IS_INTERFACE 1)
  else()
    set(REAK_CC_LIB_IS_INTERFACE 0)
  endif()
  
  foreach(hdr_file IN LISTS REAK_CC_HDRS)
    if(${hdr_file} MATCHES ".*\\.(h|inc|hpp|tpp)")
      string(FIND "${CMAKE_CURRENT_SOURCE_DIR}/${hdr_file}" "ReaK/" _hdr_file_reak_pos REVERSE)
      if(${_hdr_file_reak_pos} LESS "0")
        string(SUBSTRING "${CMAKE_CURRENT_SOURCE_DIR}/${hdr_file}" ${_hdr_file_reak_pos} -1 REAK_CC_LIB_BUILD_INCLUDE_DIR)
        break()
      endif()
    endif()
  endforeach()

  #  set(_build_type "shared")
  set(_build_type "static")

  # Generate a pkg-config file for every library:
  if(REAK_ENABLE_INSTALL)
    set(PC_VERSION "head")
    set(LNK_LIB "${LNK_LIB} -l${_NAME}")
    foreach(cflag ${REAK_CC_LIB_COPTS})
      set(PC_CFLAGS "${PC_CFLAGS} ${cflag}")
    endforeach()
    string(REPLACE ";" " " PC_LINKOPTS "${REAK_CC_LIB_LINKOPTS}")
    file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/lib/pkgconfig/${_NAME}.pc" CONTENT "\
prefix=${CMAKE_INSTALL_PREFIX}\n\
exec_prefix=\${prefix}\n\
libdir=${CMAKE_INSTALL_FULL_LIBDIR}\n\
includedir=${CMAKE_INSTALL_FULL_INCLUDEDIR}\n\
\n\
Name: ${_NAME}\n\
Description: ReaK ${REAK_CC_LIB_NAME} library\n\
URL: https://github.com/mikael-s-persson/ReaK/\n\
Version: ${PC_VERSION}\n\
Requires:${PC_DEPS}\n\
Libs: -L\${libdir} ${PC_LINKOPTS} $<$<NOT:$<BOOL:${REAK_CC_LIB_IS_INTERFACE}>>:${LNK_LIB}>\n\
Cflags: -I\${includedir}${PC_CFLAGS}\n")
    install(FILES "${CMAKE_BINARY_DIR}/lib/pkgconfig/${_NAME}.pc"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  endif()

  if(NOT REAK_CC_LIB_IS_INTERFACE)
    if(_build_type STREQUAL "static" OR _build_type STREQUAL "shared")
      add_library(${_NAME} "")
      target_sources(${_NAME} PRIVATE ${REAK_CC_LIB_SRCS} ${REAK_CC_LIB_HDRS})
      target_link_libraries(${_NAME}
      PUBLIC ${REAK_CC_LIB_DEPS}
      PRIVATE
        ${REAK_CC_LIB_LINKOPTS}
      )
    else()
      message(FATAL_ERROR "Invalid build type: ${_build_type}")
    endif()

    # Linker language can be inferred from sources, but because "CXX" is always the
    # correct linker language for static or for shared libraries, we set it
    # unconditionally.
    set_property(TARGET ${_NAME} PROPERTY LINKER_LANGUAGE "CXX")

    target_include_directories(${_NAME} ${REAK_INTERNAL_INCLUDE_WARNING_GUARD}
      PUBLIC
        "$<BUILD_INTERFACE:${REAK_CC_LIB_BUILD_INCLUDE_DIR}>"
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
    target_compile_options(${_NAME} PRIVATE ${REAK_CC_LIB_COPTS})
    target_compile_definitions(${_NAME} PUBLIC ${REAK_CC_LIB_DEFINES})

    # When being installed, we lose the absl_ prefix.  We want to put it back
    # to have properly named lib files.  This is a no-op when we are not being
    # installed.
    if(REAK_ENABLE_INSTALL)
      set_target_properties(${_NAME} PROPERTIES
        OUTPUT_NAME "${_NAME}"
        SOVERSION 0
      )
    endif()
  else()
    # Generating header-only library
    add_library(${_NAME} INTERFACE)
    target_include_directories(${_NAME} ${REAK_INTERNAL_INCLUDE_WARNING_GUARD}
      INTERFACE
        "$<BUILD_INTERFACE:${REAK_CC_LIB_BUILD_INCLUDE_DIR}>"
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      )

    target_link_libraries(${_NAME}
      INTERFACE
        ${REAK_CC_LIB_DEPS}
        ${REAK_CC_LIB_LINKOPTS}
    )
    target_compile_definitions(${_NAME} INTERFACE ${REAK_CC_LIB_DEFINES})
  endif()
  
  if(REAK_ENABLE_CLANG_TIDY)
    set_target_properties(${_NAME} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  endif()

  if(REAK_ENABLE_INSTALL)
    install(TARGETS ${_NAME} EXPORT ${PROJECT_NAME}Targets
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  endif()

  add_library(ReaK::${REAK_CC_LIB_NAME} ALIAS ${_NAME})

  message(STATUS "Registered cc_library: ReaK::${REAK_CC_LIB_NAME} (internally: ${_NAME}).")
endfunction()

# reak_cc_test()
#
# CMake function to imitate Bazel's cc_test rule.
#
# Parameters:
# NAME: name of target (see Usage below)
# SRCS: List of source files for the binary
# DEPS: List of other libraries to be linked in to the binary targets
# COPTS: List of private compile options
# DEFINES: List of public defines
# LINKOPTS: List of link options
# TESTDATA: List of data files needed for the test.
#
# Note:
# By default, reak_cc_test will always create a binary named reak_${NAME}.
# This will also add it to ctest list as reak_${NAME}.
#
# Usage:
# reak_cc_test(
#   NAME
#     awesome_test
#   SRCS
#     "awesome_test.cc"
#   DEPS
#     absl::awesome
#     GTest::gmock
#     GTest::gtest_main
# )
function(reak_cc_test)
  cmake_parse_arguments(REAK_CC_TEST
    ""
    "NAME"
    "SRCS;COPTS;DEFINES;LINKOPTS;DEPS;TESTDATA"
    ${ARGN}
  )

  set(_NAME "reak_${REAK_CC_TEST_NAME}")

  add_executable(${_NAME} "")
  target_sources(${_NAME} PRIVATE ${REAK_CC_TEST_SRCS})

  target_compile_definitions(${_NAME} PUBLIC ${REAK_CC_TEST_DEFINES})
  target_compile_options(${_NAME} PRIVATE ${REAK_CC_TEST_COPTS})

  target_link_libraries(${_NAME} PUBLIC ${REAK_CC_TEST_DEPS} PRIVATE ${REAK_CC_TEST_LINKOPTS})

  add_test(NAME ${_NAME} COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${_NAME}")
  
  foreach(_testdata_file ${REAK_CC_TEST_TESTDATA})
    get_filename_component(_testdata_file_dir ${_testdata_file} DIRECTORY)
    add_custom_command(TARGET ${_NAME} PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E 
      make_directory "${CMAKE_CURRENT_BINARY_DIR}/${_testdata_file_dir}")
    add_custom_command(TARGET ${_NAME} PRE_BUILD
      COMMAND ${CMAKE_COMMAND} -E 
      copy "${CMAKE_CURRENT_SOURCE_DIR}/${_testdata_file}" "${CMAKE_CURRENT_BINARY_DIR}/${_testdata_file_dir}")
  endforeach()
  
  if(REAK_ENABLE_CLANG_TIDY)
    set_target_properties(${_NAME} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  endif()
  
  message(STATUS "Registered cc_test: ${_NAME}")
endfunction()

# reak_cc_binary()
#
# CMake function to imitate Bazel's cc_binary rule.
#
# Parameters:
# NAME: name of target (see Usage below)
# SRCS: List of source files for the binary
# DEPS: List of other libraries to be linked in to the binary targets
# COPTS: List of private compile options
# DEFINES: List of public defines
# LINKOPTS: List of link options
#
# Note:
# By default, reak_cc_binary will always create a binary named ${NAME}.
#
# Usage:
# reak_cc_binary(
#   NAME
#     awesome_tool
#   SRCS
#     "awesome_tool.cc"
#   DEPS
#     absl::awesome
#     GTest::gmock
#     GTest::gtest_main
# )
function(reak_cc_binary)
  cmake_parse_arguments(REAK_CC_BIN
    ""
    "NAME"
    "SRCS;COPTS;DEFINES;LINKOPTS;DEPS"
    ${ARGN}
  )

  set(_NAME "${REAK_CC_BIN_NAME}")

  add_executable(${_NAME} "")
  target_sources(${_NAME} PRIVATE ${REAK_CC_BIN_SRCS})

  target_compile_definitions(${_NAME} PUBLIC ${REAK_CC_BIN_DEFINES})
  target_compile_options(${_NAME} PRIVATE ${REAK_CC_BIN_COPTS})

  target_link_libraries(${_NAME} PUBLIC ${REAK_CC_BIN_DEPS} PRIVATE ${REAK_CC_BIN_LINKOPTS})
  
  if(REAK_ENABLE_CLANG_TIDY)
    set_target_properties(${_NAME} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  endif()

  message(STATUS "Registered cc_binary: ${REAK_CC_BIN_NAME}")
endfunction()




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

set(CMAKE_CXX_STANDARD 17)


enable_testing()

include(CheckCXXCompilerFlag)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libc++ -pthread -ftemplate-depth=2000 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
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
  add_definitions( "-DREAK_HAS_BOOST" )
  link_directories(${Boost_LIBRARY_DIRS})
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
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
  set_target_properties(${target_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  install(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT ReaK_${REAK_CURRENT_MODULE}_tools)
  message(STATUS "Registered ReaK tool program ${target_name}.")
endmacro(ReaK_setup_tool_program)

macro(ReaK_setup_static_library target_name)
  set_target_properties(${target_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  install(TARGETS ${target_name} ARCHIVE DESTINATION lib COMPONENT ReaK_${REAK_CURRENT_MODULE})
  message(STATUS "Registered ReaK static library ${target_name}.")
endmacro(ReaK_setup_static_library)

macro(ReaK_setup_shared_library target_name)
  set_target_properties(${target_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  install(TARGETS ${target_name} LIBRARY DESTINATION lib COMPONENT ReaK_${REAK_CURRENT_MODULE})
  message(STATUS "Registered ReaK shared library ${target_name}.")
endmacro(ReaK_setup_shared_library)

macro(ReaK_setup_target target_name)
  set_target_properties(${target_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples")
  message(STATUS "Registered ReaK target ${target_name}.")
endmacro(ReaK_setup_target)

macro(ReaK_setup_test_program target_name)
  set_target_properties(${target_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
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




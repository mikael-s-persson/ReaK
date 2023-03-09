cmake_minimum_required(VERSION 3.14)

project(ReaK)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_mods")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(DOXROOT_PATH "${CMAKE_BINARY_DIR}/dox")
set(DOXYFILE_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in")

if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  option(REAK_ENABLE_INSTALL "Enable install rule" OFF)
  set(REAK_INTERNAL_INCLUDE_WARNING_GUARD SYSTEM)
else()
  option(REAK_ENABLE_INSTALL "Enable install rule" ON)
  set(REAK_INTERNAL_INCLUDE_WARNING_GUARD "")
endif()

option(REAK_ENABLE_CLANG_TIDY "Enable clang-tidy rule" OFF)

include(SetupReaKMacros)

# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fprofile-use" )
# set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fprofile-use" )
# set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -march=native -flto" )
# set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -O3 -march=native -flto" )

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "TGZ;DEB")
set(CPACK_DEB_COMPONENT_INSTALL 1)
set(CPACK_PACKAGE_NAME "reak")
set(CPACK_PACKAGE_VENDOR "Sven Mikael Persson")
set(CPACK_PACKAGE_CONTACT "Sven Mikael Persson <mikael.s.persson@gmail.com>")
# set(CPACK_PACKAGE_DESCRIPTION "ReaK library package description here.")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "28")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
# set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "ReaK library debian package description here.")
include(CPack)


include_directories("./")

add_subdirectory("ReaK/core")
add_subdirectory("ReaK/math")
add_subdirectory("ReaK/geometry")
add_subdirectory("ReaK/mbd")
add_subdirectory("ReaK/control")
add_subdirectory("ReaK/planning")
add_subdirectory("ReaK/topologies")



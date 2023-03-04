cmake_minimum_required(VERSION 3.14)

project(ReaK)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_mods")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(DOXROOT_PATH "${CMAKE_BINARY_DIR}/dox")
set(DOXYFILE_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in")

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


if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
  add_definitions( "-DREAK_HAS_BOOST" )
endif()

# These are the main sub-directories
get_subdir_list(SUBDIRS "libs")
message(STATUS "Obtained the following list of sub-dirs: ${SUBDIRS}")
foreach(subdir ${SUBDIRS})
  add_subdirectory("libs/${subdir}")
endforeach(subdir)



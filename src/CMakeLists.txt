cmake_minimum_required(VERSION 2.8.8)


set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# This defines the basic input (src) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}/src")
set(SRCINCLUDEROOT "${CODEROOT}/src")
# set(SRCINCLUDEROOT "${CODEROOT}/include")

set(DOXROOT_PATH "${CODEROOT}/dox")
set(DOXTOSRCDIR_PATH "../src")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SRCROOT}/CMakeMods)

if(NOT EXISTS "${OUTPUTROOT}")
  set(OUTPUTROOT "${CODEROOT}")
endif()
message(STATUS "Configured to output to: ${OUTPUTROOT}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${OUTPUTROOT}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${OUTPUTROOT}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OUTPUTROOT}/bin")

include(SetupReaKMacros)

# These are the main sub-directories
set(REAKDIR "/ReaK")

set(SUBDIRS ${REAKDIR})

# Add the sub-directories to cmake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)






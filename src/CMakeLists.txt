cmake_minimum_required(VERSION 2.6)

# Use the following root directory on wopa laptop
# set(CODEROOT "/home/wopa/mpersson/mpersson/code")
# Use the following root directory on Mikael Persson's home computer
# set(CODEROOT "/media/S_M_P_0/Ecole/TKKThesis/mpersson/code")
set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# This defines the basic input (src) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}/src")


set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")

set(LIBROOT "${CODEROOT}/lib")
set(INCLUDEROOT "${CODEROOT}/include")
set(BINROOT "${CODEROOT}/bin")

# Now set the global macros for setting up targets.
macro(setup_custom_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}")
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}" LIBRARY DESTINATION "${LIBROOT}" ARCHIVE DESTINATION "${LIBROOT}")
endmacro(setup_custom_target)

macro(setup_headers HEADER_FILES HEADER_PATH)
  foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
    install(FILES "${SRCROOT}${CURRENT_HEADER_FILE}" DESTINATION "${INCLUDEROOT}${HEADER_PATH}")
  endforeach(CURRENT_HEADER_FILE)
endmacro(setup_headers)

# These are the main sub-directories
set(REAKDIR "/ReaK")

set(SUBDIRS ${REAKDIR})

# Add the sub-directories to cmake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)






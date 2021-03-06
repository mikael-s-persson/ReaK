cmake_minimum_required(VERSION 2.8.8)


set(CODEROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

# This defines the basic input (src) and output (lib & include) directories.
set(SRCROOT "${CODEROOT}/src")
set(SRCINCLUDEROOT "${CODEROOT}/src")
# set(SRCINCLUDEROOT "${CODEROOT}/include")
set(DOXROOT "${CODEROOT}/dox")
set(DOXTOSRCDIR "../src")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SRCROOT}/CMakeMods)


if(NOT EXISTS "${DOXROOT}")
  file(MAKE_DIRECTORY "${DOXROOT}")
endif()


if(NOT EXISTS "${OUTPUTROOT}")
  set(OUTPUTROOT "${CODEROOT}")
endif()
message(STATUS "Configured to output to: ${OUTPUTROOT}")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Configured for build-type: ${CMAKE_BUILD_TYPE}")


enable_testing()

if (WIN32)
  if (MSVC)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
    set(LIBROOT "${OUTPUTROOT}/lib_msvc")
    set(INCLUDEROOT "${OUTPUTROOT}/include")
    set(BINROOT "${OUTPUTROOT}/bin_msvc")
    message(STATUS "Configured compiler options and output directories for MSVC toolset.")
  else()
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11 -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pthread -std=c++11 -Wl,--no-as-needed -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pthread -std=c++11 -Wl,--no-as-needed -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} --enable-stdcall-fixup")
    set(LIBROOT "${OUTPUTROOT}/lib_w32")
    set(INCLUDEROOT "${OUTPUTROOT}/include")
    set(BINROOT "${OUTPUTROOT}/bin_w32")
    message(STATUS "Configured compiler options and output directories for MinGW toolset.")
  endif()
else()
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -pthread -std=c++11 -Wl,--no-as-needed -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -pthread -std=c++11 -Wl,--no-as-needed -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -pthread -std=c++11 -nostdinc++ -isystem /usr/include/c++/4.8 -isystem /usr/include/x86_64-linux-gnu/c++/4.8 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -pthread -std=c++11 -nostdinc++ -isystem /usr/include/c++/4.8 -isystem /usr/include/x86_64-linux-gnu/c++/4.8 -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  else()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -pthread -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -pthread -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  endif()
  set(LIBROOT "${OUTPUTROOT}/lib")
  set(INCLUDEROOT "${OUTPUTROOT}/include")
  set(BINROOT "${OUTPUTROOT}/bin")
  message(STATUS "Configured compiler options and output directories for *nix GCC toolset.")
endif()

# Git related macros and commands:
find_package(Git)

if(GIT_FOUND)
  message(STATUS "Git package was found: '${GIT_EXECUTABLE}'")
  set(PRINT_VERSION_COMMAND ${GIT_EXECUTABLE} describe)
  set(PRINT_FILE_DATE_COMMAND_STR "${GIT_EXECUTABLE} log --pretty=format:'Author: %an <%ae>%nDate: %ad' --date=short -n 1")
else()
  set(PRINT_VERSION_COMMAND ${CMAKE_COMMAND} -E echo unknown version)
  set(PRINT_FILE_DATE_COMMAND_STR "${CMAKE_COMMAND} -E echo")
endif()


# Now set the global macros for setting up targets.
macro(setup_custom_target target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}")
  set_property(TARGET ${target_name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${LIBROOT}")
  set_property(TARGET ${target_name} PROPERTY ARCHIVE_OUTPUT_DIRECTORY "${LIBROOT}")
  get_target_property(is_target_excluded_from_all ${target_name} EXCLUDE_FROM_ALL)
  if( NOT is_target_excluded_from_all )
    install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}" LIBRARY DESTINATION "${LIBROOT}" ARCHIVE DESTINATION "${LIBROOT}")
  endif()
  message(STATUS "Registered target ${target_name}.")
endmacro(setup_custom_target)

macro(setup_custom_test_program target_name CUSTOM_OUTPUT_PATH)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}/unit_tests")
  add_test(NAME "${target_name}" WORKING_DIRECTORY "${BINROOT}/unit_tests/" COMMAND "${BINROOT}/unit_tests/${target_name}")
endmacro(setup_custom_test_program)

macro(setup_headers HEADER_FILES HEADER_PATH)
  foreach(CURRENT_HEADER_FILE ${HEADER_FILES})
    install(FILES "${SRCROOT}${CURRENT_HEADER_FILE}" DESTINATION "${INCLUDEROOT}${HEADER_PATH}")
  endforeach(CURRENT_HEADER_FILE)
endmacro(setup_headers)

# Doxygen related macros and commands:
find_package(Doxygen)

macro(configure_doxyfile TARGET_NAME TARGET_TITLE TARGET_BRIEF TARGET_INPUT_PATH)
  if(EXISTS "${SRCROOT}/Doxyfile.in")
    set(DOXY_TARGET_OUTPUT_DIR "${TARGET_NAME}")
    set(DOXY_TARGET_ROOT_DIR "${DOXTOSRCDIR}")
    set(DOXY_TARGET_NAME "${TARGET_NAME}")
    set(DOXY_TARGET_TITLE "${TARGET_TITLE}")
    set(DOXY_TARGET_BRIEF "${TARGET_BRIEF}")
    set(DOXY_TARGET_INPUT_PATH "${TARGET_INPUT_PATH}")
    set(DOXY_TARGET_FILE_VERSION_FILTER "\"${PRINT_FILE_DATE_COMMAND_STR} \"")
    if(NOT EXISTS "${DOXROOT}/${DOXY_TARGET_OUTPUT_DIR}")
      file(MAKE_DIRECTORY "${DOXROOT}/${DOXY_TARGET_OUTPUT_DIR}")
    endif()
    configure_file("${SRCROOT}/Doxyfile.in" "${DOXROOT}/${TARGET_NAME}_Doxyfile" @ONLY)
    message(STATUS "Configured Doxyfile: '${SRCROOT}/Doxyfile.in' --> '${DOXROOT}/${TARGET_NAME}_Doxyfile'.")
  else()
    message(STATUS "WARNING : The '${SRCROOT}/Doxyfile.in' file does not exist!")
  endif()
endmacro(configure_doxyfile)

if (DOXYGEN_FOUND)
  message(STATUS "Doxygen package was found.")
  add_custom_target(dox)
  if(WIN32)
    find_package(HTMLHelp)
    if(HTML_HELP_COMPILER)
      message(STATUS "Compressed help-file generator was found: '${HTML_HELP_COMPILER}'")
      set(DOXY_COMPRESSED_HELP_GEN ${HTML_HELP_COMPILER})
    else()
      message(STATUS "WARNING : Compressed help-files cannot be generated on this platform (hhc.exe is missing).")
    endif()
  else()
    find_package(Qt4)
    if(QT4_FOUND)
      message(STATUS "Qt4 package was found.")
      find_program(QT_QHELPGENERATOR_EXECUTABLE
                   NAMES qhelpgenerator-qt4 qhelpgenerator
                   PATHS ${QT_BINARY_DIR} NO_DEFAULT_PATH)
      if(QT_QHELPGENERATOR_EXECUTABLE)
        message(STATUS "QHelpGenerator executable was found: '${QT_QHELPGENERATOR_EXECUTABLE}'")
        set(DOXY_COMPRESSED_HELP_GEN ${QT_QHELPGENERATOR_EXECUTABLE})
      else()
        message(STATUS "WARNING : Compressed help-files cannot be generated on this platform (qhelpgenerator is missing).")
      endif()
    endif()
  endif()
  macro(add_doxygen_target TARGET_NAME)
    if(WIN32)
      if(HTML_HELP_COMPILER)
        set(DOXY_COMPRESSED_HELP_GEN ${HTML_HELP_COMPILER})
        set(DOXY_COMPRESSED_HELP_INPUT "${TARGET_NAME}/html/index.hhp")
      else()
        set(DOXY_COMPRESSED_HELP_GEN ${CMAKE_COMMAND} -E echo)
        set(DOXY_COMPRESSED_HELP_INPUT "Skipping compressed help-file generation for '${TARGET_NAME}'...")
      endif()
    else()
      if(QT_QHELPGENERATOR_EXECUTABLE)
        set(DOXY_COMPRESSED_HELP_GEN ${QT_QHELPGENERATOR_EXECUTABLE})
        set(DOXY_COMPRESSED_HELP_INPUT "${TARGET_NAME}/html/index.qhp" -o "${TARGET_NAME}_dox.qch")
      else()
        set(DOXY_COMPRESSED_HELP_GEN ${CMAKE_COMMAND} -E echo)
        set(DOXY_COMPRESSED_HELP_INPUT "Skipping compressed help-file generation for '${TARGET_NAME}'...")
      endif()
    endif()
    add_custom_target(${TARGET_NAME}_dox 
                      ${CMAKE_COMMAND} -E remove -f "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${CMAKE_COMMAND} -E copy "${TARGET_NAME}_Doxyfile" "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${PRINT_VERSION_COMMAND} >> "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${DOXYGEN_EXECUTABLE} "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${CMAKE_COMMAND} -E remove -f "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${DOXY_COMPRESSED_HELP_GEN} ${DOXY_COMPRESSED_HELP_INPUT}
              WORKING_DIRECTORY "${DOXROOT}" VERBATIM)
    add_dependencies(dox ${TARGET_NAME}_dox)
  endmacro()
else()
  message(STATUS "WARNING : Doxygen package is not available on this system!")
  macro(add_doxygen_target TARGET_NAME)
  endmacro()
endif()

# These are the main sub-directories
set(REAKDIR "/ReaK")

set(SUBDIRS ${REAKDIR})

# Add the sub-directories to cmake process.
foreach(CURRENT_SUBDIR ${SUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_SUBDIR}")
endforeach(CURRENT_SUBDIR)






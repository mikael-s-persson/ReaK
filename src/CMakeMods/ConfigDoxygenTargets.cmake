# Tries to find Doxygen and version control tools, and then create a 
# Doxygen configuration functions that allows the user to use an input
# Doxygen file and configure it to pull in the headers of a particular 
# cmake project. 
# 
# This module expects two cache variables to be set:
#  DOXROOT_PATH : the root path for the doxygen output.
#  DOXTOSRC_PATH : the relative path from the doxygen root output folder to the source (src) folder.
# 
#  N.B.: If either or both of these are not set prior to including this module, 
#  then it will be assumed that the CMAKE_SOURCE_DIR is the "src" folder, and that the 
#  doxygen output should be made to ${CMAKE_SOURCE_DIR}/../dox (i.e., alongside the "src" folder).
# 
# The custom functions provided by this file are:
# 
#  configure_doxyfile(TARGET_NAME TARGET_TITLE TARGET_BRIEF TARGET_INPUT_PATH)
#  
#    This function allows the user the create doxygen configuration file from a Doxyfile.in
#    that is present in the top-level source directory and a set of simple parameters for the 
#    documentation. The TARGET_NAME is the name to be given to the doxygen target that will be 
#    associated to this configured doxygen file. The TARGET_TITLE is the title that will appear
#    on the title page of the generated Doxygen documentation. The TARGET_BRIEF is the brief 
#    description of the overall library / documentation. The TARGET_INPUT_PATH parameter is a 
#    list of paths to be included into the documentation.
#  
#  add_doxygen_target(TARGET_NAME)
#  
#    This function registers a doxygen target named ${TARGET_NAME}_dox as well as registers 
#    it as a dependency for the global "dox" target that generates all the doxygen targets.
# 
# EXAMPLE:
#  
#  configure_doxyfile(ReaKcore 
#    "\"ReaK Library - Core libraries\"" 
#    "\"Core libraries for software architectural elements of the ReaK libraries.\"" 
#    "ReaK/core/base;ReaK/core/rtti;ReaK/core/serialization;ReaK/core/recorders;ReaK/core/sorting")
#  add_doxygen_target(ReaKcore)
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



if(NOT EXISTS "${DOXROOT_PATH}")
  set(DOXROOT_PATH "${CMAKE_SOURCE_DIR}/../dox")
endif()

if(NOT EXISTS "${DOXTOSRC_PATH}")
  set(DOXTOSRC_PATH "../src")
endif()

if(NOT EXISTS "${DOXROOT_PATH}")
  file(MAKE_DIRECTORY "${DOXROOT_PATH}")
endif()
message(STATUS "Configured to Doxygen output to: ${DOXROOT_PATH}")


# Git related macros and commands:
find_package(Git)
if(GIT_FOUND)
  message(STATUS "Git package was found: '${GIT_EXECUTABLE}'")
  
  execute_process(COMMAND ${GIT_EXECUTABLE} describe
                  RESULT_VARIABLE GIT_DESCRIBE_RESULT_VAL
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  ERROR_QUIET
                  OUTPUT_QUIET
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(${GIT_DESCRIBE_RESULT_VAL} EQUAL 0)
    message(STATUS "Detected the source directory is a Git repository.")
    set(CURRENTLY_IN_GIT_REPO "TRUE")
    set(CURRENTLY_IN_VERSION_CONTROL "TRUE")
    set(REPO_PRINT_VERSION_COMMAND ${GIT_EXECUTABLE} describe)
    set(REPO_PRINT_FILE_DATE_COMMAND_STR "${GIT_EXECUTABLE} log --pretty=format:'Author: %an <%ae>%nDate: %ad' --date=short -n 1")
  endif()
  unset(GIT_DESCRIBE_RESULT_VAL)
endif()  

if(NOT CURRENTLY_IN_VERSION_CONTROL)
  find_package(Subversion)
  if(Subversion_FOUND)
    execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} info ${CMAKE_SOURCE_DIR}
      RESULT_VARIABLE Subversion_INFO_RESULT_VAL
      ERROR_QUIET
      OUTPUT_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${Subversion_INFO_RESULT_VAL} EQUAL 0)
      message(STATUS "Detected the source directory is a subversion repository.")
      Subversion_WC_INFO(${CMAKE_SOURCE_DIR} project_svn_info)
      set(CURRENTLY_IN_SUBVERSION_REPO "TRUE")
      set(CURRENTLY_IN_VERSION_CONTROL "TRUE")
      set(REPO_PRINT_VERSION_COMMAND ${CMAKE_COMMAND} -E echo ${project_svn_info_WC_REVISION})
      set(REPO_PRINT_FILE_DATE_COMMAND_STR "${Subversion_SVN_EXECUTABLE} log -q -l 1")
      unset(project_svn_info_WC_URL)
      unset(project_svn_info_WC_ROOT)
      unset(project_svn_info_WC_REVISION)
      unset(project_svn_info_WC_LAST_CHANGED_AUTHOR)
      unset(project_svn_info_WC_LAST_CHANGED_DATE)
      unset(project_svn_info_WC_LAST_CHANGED_REV)
      unset(project_svn_info_WC_INFO)
    endif()
    unset(Subversion_INFO_RESULT_VAL)
  endif()
endif()

if(NOT CURRENTLY_IN_VERSION_CONTROL)
  find_package(CVS)
  if(CVS_FOUND)
    execute_process(COMMAND ${CVS_EXECUTABLE} status -n -d ${CMAKE_SOURCE_DIR}
      RESULT_VARIABLE CVS_STATUS_RESULT_VAL
      ERROR_QUIET
      OUTPUT_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${CVS_STATUS_RESULT_VAL} EQUAL 0)
      message(STATUS "Detected the source directory is a CVS repository.")
      set(CURRENTLY_IN_CVS_REPO "TRUE")
      set(CURRENTLY_IN_VERSION_CONTROL "TRUE")
      
      # TODO: Implement the CVS revision command.
      
      set(REPO_PRINT_VERSION_COMMAND ${CMAKE_COMMAND} -E echo unknown version)
      set(REPO_PRINT_FILE_DATE_COMMAND_STR "${CMAKE_COMMAND} -E echo")
    endif()
    unset(CVS_STATUS_RESULT_VAL)
  endif()
endif()

if(NOT CURRENTLY_IN_VERSION_CONTROL)
  set(REPO_PRINT_VERSION_COMMAND ${CMAKE_COMMAND} -E echo unknown version)
  set(REPO_PRINT_FILE_DATE_COMMAND_STR "${CMAKE_COMMAND} -E echo")
endif()


# Doxygen related macros and commands:
find_package(Doxygen)


macro(configure_doxyfile TARGET_NAME TARGET_TITLE TARGET_BRIEF TARGET_INPUT_PATHS)
  if(EXISTS "${CMAKE_SOURCE_DIR}/Doxyfile.in")
    set(DOXY_TARGET_OUTPUT_DIR "${TARGET_NAME}")
    set(DOXY_TARGET_ROOT_DIR "${DOXTOSRC_PATH}")
    set(DOXY_TARGET_NAME "${TARGET_NAME}")
    set(DOXY_TARGET_TITLE "${TARGET_TITLE}")
    set(DOXY_TARGET_BRIEF "${TARGET_BRIEF}")
    message(STATUS "Parsing input paths of Doxygen target ${TARGET_NAME} : '${TARGET_TITLE}' ...")
    foreach(INPUT_PATH ${TARGET_INPUT_PATHS})
      if(NOT EXISTS "${DOXROOT_PATH}/${DOXTOSRC_PATH}${INPUT_PATH}")
        message(STATUS "  Could NOT find path: ${DOXTOSRC_PATH}${INPUT_PATH}")
      else()
        set(DOXY_TARGET_INPUT_PATH "${DOXY_TARGET_INPUT_PATH} ${DOXTOSRC_PATH}${INPUT_PATH}")
        message(STATUS "  Found path: ${DOXTOSRC_PATH}${INPUT_PATH}")
      endif()
    endforeach(INPUT_PATH)
    set(DOXY_TARGET_FILE_VERSION_FILTER "\"${REPO_PRINT_FILE_DATE_COMMAND_STR} \"")
    if(NOT EXISTS "${DOXROOT_PATH}/${DOXY_TARGET_OUTPUT_DIR}")
      file(MAKE_DIRECTORY "${DOXROOT_PATH}/${DOXY_TARGET_OUTPUT_DIR}")
    endif()
    configure_file("${CMAKE_SOURCE_DIR}/Doxyfile.in" "${DOXROOT_PATH}/${TARGET_NAME}_Doxyfile" @ONLY)
    message(STATUS "Configured Doxyfile: '${CMAKE_SOURCE_DIR}/Doxyfile.in' --> '${DOXROOT_PATH}/${TARGET_NAME}_Doxyfile'.")
  else()
    message(STATUS "WARNING : The '${CMAKE_SOURCE_DIR}/Doxyfile.in' file does not exist!")
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
              COMMAND ${REPO_PRINT_VERSION_COMMAND} >> "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${DOXYGEN_EXECUTABLE} "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${CMAKE_COMMAND} -E remove -f "${TARGET_NAME}_Doxyfile.out"
              COMMAND ${DOXY_COMPRESSED_HELP_GEN} ${DOXY_COMPRESSED_HELP_INPUT}
              WORKING_DIRECTORY "${DOXROOT_PATH}" VERBATIM)
    add_dependencies(dox ${TARGET_NAME}_dox)
  endmacro()
else()
  message(STATUS "WARNING : Doxygen package is not available on this system!")
  macro(add_doxygen_target TARGET_NAME)
  endmacro()
endif()







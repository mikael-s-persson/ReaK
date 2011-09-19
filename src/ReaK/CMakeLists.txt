project(ReaK)

# The following definition of RK_VERBOSITY determines the depth-level of the debug option of the programs (usually between 0 (no debug output) and 10 (full debug output)).
add_definitions(-DRK_VERBOSITY=3)

# These are the main sub-directories
set(RKCOREDIR "${REAKDIR}/core")
set(RKCTRLDIR "${REAKDIR}/ctrl")
set(RKEXAMPLESDIR "${REAKDIR}/examples")
set(RKCOREOLDDIR "${REAKDIR}/core_old")
set(RKCTRLOLDDIR "${REAKDIR}/ctrl_old")
set(RKRENDERERDIR "${REAKDIR}/renderer")
set(RKAIDIR "${REAKDIR}/ai")

set(RKSUBDIRS ${RKCOREDIR} ${RKCTRLDIR} ${RKEXAMPLESDIR})
#set(RKSUBDIRS ${RKSUBDIRS} ${RKCOREOLDDIR} ${RKCTRLOLDDIR})

if (WIN32)
  set(BOOST_INCLUDE_DIRS "boost_1_46_0/")
  if (MSVC)
    set(BOOST_LIB_TOPDIR "C:/boost/boost_1_46_0/stage_msvc/lib/")
    link_directories("${BOOST_LIB_TOPDIR}")
    #    find_library(BOOST_THREAD_LIB NAMES boost_thread-vc100-mt boost_thread-vc80-mt PATHS BOOST_LIB_TOPDIR)
    #    find_library(BOOST_DATE_TIME_LIB NAMES boost_date_time-vc100-mt boost_date_time-vc80-mt boost_date_time-vc100 boost_date_time-vc80-mt PATHS BOOST_LIB_TOPDIR)
    set(BOOST_THREAD_LIB "")
    set(BOOST_DATE_TIME_LIB "")
    #    set(BOOST_THREAD_LIB "${BOOST_LIB_TOPDIR}boost_thread-vc100-mt.lib")
    #    set(BOOST_DATE_TIME_LIB "${BOOST_LIB_TOPDIR}boost_date-time-vc100-mt.lib")
  else()
    set(Boost_USE_STATIC_LIBS ON)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
    find_package( Boost 1.45.0 COMPONENTS thread date_time )
    if(Boost_FOUND)
      set(BOOST_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
      set(BOOST_THREAD_LIB ${Boost_LIBRARIES})
      set(BOOST_DATE_TIME_LIB ${Boost_LIBRARIES})
    else()
      set(BOOST_LIB_TOPDIR "boost_1_45_0/stage/lib/")
      set(BOOST_LIBRARYDIR "boost_1_45_0/stage/lib")
      set(BOOST_INCLUDEDIR "boost_1_45_0")
      find_library(BOOST_THREAD_LIB NAMES boost_thread-mt PATHS BOOST_LIBRARYDIR)
      set(BOOST_THREAD_LIB "${BOOST_LIB_TOPDIR}libboost_thread-mgw44-mt-1_45.a")
      find_library(BOOST_DATE_TIME_LIB NAMES boost_date_time-mt boost_date_time PATHS BOOST_LIBRARYDIR)
    endif()
  endif() 
else()
  set(BOOST_INCLUDE_DIRS "/usr/include/")
  set(BOOST_LIB_TOPDIR "/usr/lib/")
  find_library(BOOST_THREAD_LIB NAMES boost_thread PATHS BOOST_LIB_TOPDIR)
  find_library(BOOST_DATE_TIME_LIB NAMES boost_date_time-mt boost_date_time PATHS "/usr/lib/")
  #find_library(BOOST_SYSTEM_LIB NAMES boost_system-mt boost_system PATHS "/usr/lib/")
  #find_library(BOOST_SERIALIZATION_LIB NAMES boost_serialization-mt boost_serialization PATHS "/usr/lib/")
  #find_library(BOOST_REGEX_LIB NAMES boost_regex-mt boost_regex PATHS "/usr/lib/")
  #find_library(BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options-mt boost_program_options PATHS "/usr/lib/")
  #find_library(BOOST_FILESYSTEM_LIB NAMES boost_filesystem-mt boost_filesystem PATHS "/usr/lib/")
  #find_library(BOOST_IOSTREAMS_LIB NAMES boost_iostreams-mt boost_iostreams PATHS "/usr/lib/")
endif()




#if(NOT BOOST_THREAD_LIB)
  if(UNIX)
    set(BOOST_THREAD_LIB "/usr/lib/libboost_thread.so.1.42.0")
    #  else()
    #   message(FATAL_ERROR "Boost thread library is missing")
  endif()
  #endif()

  #if(NOT BOOST_DATE_TIME_LIB)
  #  message(FATAL_ERROR "Boost date time library is missing")
  #endif()

#if(NOT BOOST_SYSTEM_LIB)
#  message(FATAL_ERROR "Boost system library is missing")
#endif()

#if(NOT BOOST_SERIALIZATION_LIB)
#  message(FATAL_ERROR "Boost serialization library is missing")
#endif()

#if(NOT BOOST_REGEX_LIB)
#  message(FATAL_ERROR "Boost regex library is missing")
#endif()

#if(NOT BOOST_PROGRAM_OPTIONS_LIB)
#  message(FATAL_ERROR "Boost program options library is missing")
#endif()


#if(NOT BOOST_FILESYSTEM_LIB)
#  message(FATAL_ERROR "Boost filesystem library is missing")
#endif()

#if(NOT BOOST_IOSTREAMS_LIB)
#  message(FATAL_ERROR "Boost iostreams library is missing")
#endif()


# Add the sub-directories to cmake process.
foreach(CURRENT_RKSUBDIR ${RKSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKSUBDIR}")
endforeach(CURRENT_RKSUBDIR)



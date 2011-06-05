project(ReaK)

# The following definition of RK_VERBOSITY determines the depth-level of the debug option of the programs (usually between 0 (no debug output) and 10 (full debug output)).
add_definitions(-DRK_VERBOSITY=3)

# These are the main sub-directories
set(RKCOREDIR "${REAKDIR}/core")
set(RKCTRLDIR "${REAKDIR}/ctrl")
set(RKEXAMPLESDIR "${REAKDIR}/examples")

set(RKSUBDIRS ${RKCOREDIR} ${RKCTRLDIR} ${RKEXAMPLESDIR})

# Find all boost components needed.
set(BOOST_INCLUDE_DIRS "/usr/include/")

find_library(BOOST_THREAD_LIB NAMES boost_thread PATHS "/usr/lib/")
if(NOT BOOST_THREAD_LIB)
  message(FATAL_ERROR "Boost thread library is missing")
endif()
set(BOOST_THREAD_LIB "/usr/lib/libboost_thread.so.1.42.0")

find_library(BOOST_SYSTEM_LIB NAMES boost_system-mt boost_system PATHS "/usr/lib/")
if(NOT BOOST_SYSTEM_LIB)
  message(FATAL_ERROR "Boost system library is missing")
endif()

find_library(BOOST_SERIALIZATION_LIB NAMES boost_serialization-mt boost_serialization PATHS "/usr/lib/")
if(NOT BOOST_SERIALIZATION_LIB)
  message(FATAL_ERROR "Boost serialization library is missing")
endif()

#find_library(BOOST_REGEX_LIB NAMES boost_regex-mt boost_regex PATHS "/usr/lib/")
#if(NOT BOOST_REGEX_LIB)
#  message(FATAL_ERROR "Boost regex library is missing")
#endif()

#find_library(BOOST_PROGRAM_OPTIONS_LIB NAMES boost_program_options-mt boost_program_options PATHS "/usr/lib/")
#if(NOT BOOST_PROGRAM_OPTIONS_LIB)
#  message(FATAL_ERROR "Boost program options library is missing")
#endif()

find_library(BOOST_DATE_TIME_LIB NAMES boost_date_time-mt boost_date_time PATHS "/usr/lib/")
if(NOT BOOST_DATE_TIME_LIB)
  message(FATAL_ERROR "Boost date time library is missing")
endif()

find_library(BOOST_FILESYSTEM_LIB NAMES boost_filesystem-mt boost_filesystem PATHS "/usr/lib/")
if(NOT BOOST_FILESYSTEM_LIB)
  message(FATAL_ERROR "Boost filesystem library is missing")
endif()

#find_library(BOOST_IOSTREAMS_LIB NAMES boost_iostreams-mt boost_iostreams PATHS "/usr/lib/")
#if(NOT BOOST_IOSTREAMS_LIB)
#  message(FATAL_ERROR "Boost iostreams library is missing")
#endif()


# Add the sub-directories to cmake process.
foreach(CURRENT_RKSUBDIR ${RKSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKSUBDIR}")
endforeach(CURRENT_RKSUBDIR)



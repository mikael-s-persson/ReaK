project(ReaK)

# The following definition of RK_VERBOSITY determines the depth-level of the debug option of the programs (usually between 0 (no debug output) and 10 (full debug output)).
add_definitions(-DRK_VERBOSITY=3)

if( ENABLE_EXTERN_TEMPLATES )
  add_definitions(-DRK_ENABLE_EXTERN_TEMPLATES)
endif()

set(REAK_VERSION_NUMBER "0.26")

#set(EXTRA_SYSTEM_LIBS)

# These are the main sub-directories
set(RKCOREDIR "${REAKDIR}/core")
set(RKCTRLDIR "${REAKDIR}/ctrl")
set(RKEXAMPLESDIR "${REAKDIR}/examples")
set(RKGEOMETRYDIR "${REAKDIR}/geometry")
set(RKGUIQTDIR "${REAKDIR}/gui_qt")

set(RKSUBDIRS ${RKCOREDIR} ${RKCTRLDIR} ${RKEXAMPLESDIR} ${RKGEOMETRYDIR} ${RKGUIQTDIR})

configure_doxyfile(ReaK "\"ReaK Library - All libraries\"" "\"A complete library for robotics modeling, control and path-planning.\"" "${DOXTOSRCDIR}${RKCOREDIR} ${DOXTOSRCDIR}${RKCTRLDIR} ${DOXTOSRCDIR}${RKEXAMPLESDIR}")
add_doxygen_target(ReaK)

if( NOT WIN32 )
  find_package(OpenCV)

  if(OpenCV_FOUND)
    message(STATUS "OpenCV library components were found: '${OpenCV_LIBS}'")
    add_definitions( "-DREAK_HAS_OPENCV" )
  else()
    message(WARNING "OpenCV was not detected on this system, some test programs will not be built!")
  endif()

else()
  
  find_library( WINSOCK_LIB wsock32 PATHS "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib")
  if( NOT WINSOCK_LIB )
    message( FATAL_ERROR  "Could not located the winsock library (ws2_32)")
  else()
    set(EXTRA_SYSTEM_LIBS ${EXTRA_SYSTEM_LIBS} ${WINSOCK_LIB})
  endif()
  
endif()


set(Boost_ADDITIONAL_VERSIONS "1.45" "1.45.0" "1.46" "1.46.0" "1.46.1" "1.47" "1.47.0" "1.48" "1.48.0" "1.49" "1.49.0" "1.50" "1.50.0" "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0" "1.54" "1.54.0" "1.55" "1.55.0")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

# if a custom path for boost is provided, than use that (and suppress system paths).
if(EXISTS "${CUSTOM_BOOST_PATH}/include/boost")
  set(Boost_INCLUDE_DIR "${CUSTOM_BOOST_PATH}/include")
  set(Boost_LIBRARY_DIR "${CUSTOM_BOOST_PATH}/lib")
  set(Boost_NO_SYSTEM_PATHS TRUE)
endif()

if (NOT WIN32)

  # make sure that the *nix suffixes and prefixes are correct (some cmake installs of findBoost.cmake are wrong with this).
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES})
  if( Boost_USE_STATIC_LIBS )
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})

endif()

find_package(Boost 1.48 COMPONENTS thread chrono date_time system program_options unit_test_framework filesystem random REQUIRED)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
  add_definitions( "-DREAK_HAS_BOOST" )
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
  set(BOOST_THREAD_LIB ${Boost_thread_LIBRARY})
  set(BOOST_DATE_TIME_LIB ${Boost_date_time_LIBRARY})
  if(MSVC)
    #disable the libraries, since it uses automatic linking:
	set(Boost_LIBRARIES "")
  endif()
endif()

if( NOT WIN32 )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_PREFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES})
endif()


find_package(BGLExtra REQUIRED)
if(BGLEXTRA_FOUND)
  message(STATUS "BGL-Extra was detected on this system: '${BGLEXTRA_INCLUDE_DIR}'")
  include_directories(SYSTEM ${BGLEXTRA_INCLUDE_DIR})
endif()



find_package(Coin3D)

if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
  add_definitions( "-DREAK_HAS_COIN3D" )
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

find_package( Qt4 )

if(COIN3D_FOUND AND QT_FOUND)
  FIND_LIBRARY(SOQT4_LIBRARIES SoQt4)
endif()




include_directories(BEFORE SYSTEM ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCINCLUDEROOT}")


# Add the sub-directories to cmake process.
foreach(CURRENT_RKSUBDIR ${RKSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKSUBDIR}")
endforeach(CURRENT_RKSUBDIR)

INCLUDE(PyReaK.cmake)


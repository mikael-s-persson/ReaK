project(ReaK)

# The following definition of RK_VERBOSITY determines the depth-level of the debug option of the programs (usually between 0 (no debug output) and 10 (full debug output)).
add_definitions(-DRK_VERBOSITY=3)

if( ENABLE_EXTERN_TEMPLATES )
  add_definitions(-DRK_ENABLE_EXTERN_TEMPLATES)
endif()

set(REAK_VERSION_NUMBER "0.26")

# These are the main sub-directories
set(RKCOREDIR "${REAKDIR}/core")
set(RKCTRLDIR "${REAKDIR}/ctrl")
set(RKEXAMPLESDIR "${REAKDIR}/examples")
set(RKCOREOLDDIR "${REAKDIR}/core_old")
set(RKCTRLOLDDIR "${REAKDIR}/ctrl_old")
set(RKRENDERERDIR "${REAKDIR}/renderer")
set(RKAIDIR "${REAKDIR}/ai")
set(RKGEOMETRYDIR "${REAKDIR}/geometry")

set(RKSUBDIRS ${RKCOREDIR} ${RKCTRLDIR} ${RKEXAMPLESDIR} ${RKGEOMETRYDIR})
#set(RKSUBDIRS ${RKSUBDIRS} ${RKCOREOLDDIR} ${RKCTRLOLDDIR})

configure_doxyfile(ReaK "\"ReaK Library - All libraries\"" "\"A complete library for robotics modeling, control and path-planning.\"" "${DOXTOSRCDIR}${RKCOREDIR} ${DOXTOSRCDIR}${RKCTRLDIR} ${DOXTOSRCDIR}${RKEXAMPLESDIR}")
add_doxygen_target(ReaK)

if( NOT WIN32 )
  find_package(OpenCV REQUIRED core highgui imgproc)

  if(OpenCV_FOUND)
    message(STATUS "OpenCV library components were found: '${OpenCV_LIBS}'")
  else()
    message(WARNING "OpenCV was not detected on this system, some test programs will not be built!")
  endif()

  find_library(FreeImage_LIB NAMES freeimage)

  if(FreeImage_LIB)
    set(FreeImage_FOUND TRUE)
    message(STATUS "FreeImage library was found: '${FreeImage_LIB}'")
  else()
    set(FreeImage_FOUND FALSE)
    message(WARNING "FreeImage library was not detected on this system, some test programs will not be built!")
  endif()

endif()


set(Boost_ADDITIONAL_VERSIONS "1.45" "1.45.0" "1.46" "1.46.0" "1.46.1" "1.47" "1.47.0" "1.48" "1.48.0")
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

find_package(Boost 1.42 COMPONENTS thread date_time python REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  add_definitions( "-DHAS_BOOST" )
  set(BOOST_THREAD_LIB ${Boost_thread_LIBRARY})
  set(BOOST_DATE_TIME_LIB ${Boost_date_time_LIBRARY})
endif()

if( NOT WIN32 )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_PREFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES})
endif()




find_package(Coin3D)

if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE ${COIN3D_INCLUDE_DIRS})
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

find_package( Qt4 )

if(COIN3D_FOUND AND QT_FOUND)
  FIND_LIBRARY(SOQT4_LIBRARIES SoQt4)
endif()




include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")


# Add the sub-directories to cmake process.
foreach(CURRENT_RKSUBDIR ${RKSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKSUBDIR}")
endforeach(CURRENT_RKSUBDIR)

INCLUDE(PyReaK.cmake)


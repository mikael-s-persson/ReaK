# project(ReaK_planning)

find_package(Coin3D)
if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

find_package( Qt4 )
if(COIN3D_FOUND AND QT_FOUND)
  find_library(SOQT4_LIBRARIES SoQt4)
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../core")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../core/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.core module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../math")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../math/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.math module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../geometry")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../geometry/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.geometry module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../mbd")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../mbd/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.mbd module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../topologies")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../topologies/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.topologies module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../control")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../control/include")
else()
  message(FATAL_ERROR "The ReaK.planning module requires the presence of the ReaK.control module!")
endif()

set(REAK_CURRENT_MODULE planning)
set(REAK_CURRENT_MODULE_TITLE "ReaK Motion-Planning Library Components")
set(REAK_CURRENT_MODULE_DESCRIPTION "A suite of graph search algorithms and sampling-based motion planners, including much of its supporting functionalities such as graph connectors, nearest-neighbor query data-structures, and so on.")

ReaK_add_current_module()




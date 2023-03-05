# project(ReaK_topologies)

# find_package(Coin3D)
# if(COIN3D_FOUND)
#   message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
#   include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
# else()
#   message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
# endif()
# 
# find_package( Qt4 )
# if(COIN3D_FOUND AND QT_FOUND)
#   find_library(SOQT4_LIBRARIES SoQt4)
# endif()


if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../core")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../core/include")
else()
  message(FATAL_ERROR "The ReaK.geom module requires the presence of the ReaK.core module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../math")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../math/include")
else()
  message(FATAL_ERROR "The ReaK.geom module requires the presence of the ReaK.math module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../geometry")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../geometry/include")
else()
  message(FATAL_ERROR "The ReaK.geom module requires the presence of the ReaK.geometry module!")
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../mbd")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../mbd/include")
else()
  message(FATAL_ERROR "The ReaK.geom module requires the presence of the ReaK.mbd module!")
endif()

find_package(OpenCV)

if(OpenCV_FOUND)
  message(STATUS "OpenCV library components were found: '${OpenCV_LIBS}'")
  add_definitions( "-DREAK_HAS_OPENCV" )
else()
  message(WARNING "OpenCV was not detected on this system, some test programs will not be built!")
endif()


set(REAK_CURRENT_MODULE topologies)
set(REAK_CURRENT_MODULE_TITLE "ReaK Topologies Library Components")
set(REAK_CURRENT_MODULE_DESCRIPTION "A library of topological representations and related methods. Topologies represent algrebraic spaces such as joint-spaces or SO(3) / SE(3). Functionality provided on topologies include various interpolation methods.")

ReaK_add_current_module()




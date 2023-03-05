# project(ReaK_geom)

find_package(Coin3D)
if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

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


set(REAK_CURRENT_MODULE geom)
set(REAK_CURRENT_MODULE_TITLE "ReaK Geometry Library Components")
set(REAK_CURRENT_MODULE_DESCRIPTION "A tool-set for geometric representations, proximity queries and rendering for 2D and 3D objects.")

ReaK_add_current_module()




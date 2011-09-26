
set(PATH_PLANNING_SOURCES 
)

set(PATH_PLANNING_HEADERS 
  "${RKPATHPLANNINGDIR}/topological_search.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_search.hpp"
  "${RKPATHPLANNINGDIR}/temporal_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/reachability_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/reachability_sort.hpp"
  "${RKPATHPLANNINGDIR}/spatial_path_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_pt_to_pt.hpp"
  "${RKPATHPLANNINGDIR}/predicted_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/differentiable_space_concept.hpp"
)

#add_library(reak_path_planning STATIC ${PATH_PLANNING_SOURCES})
#setup_custom_target(reak_path_planning "${SRCROOT}${RKPATHPLANNINGDIR}")
setup_headers("${PATH_PLANNING_HEADERS}" "${RKPATHPLANNINGDIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}/test_vp_tree.cpp")
setup_custom_target(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")

#find_package( Qt4 REQUIRED )
#include( ${QT_USE_FILE} )
#set( PATHPLANNING_PPTESTER_UIS
#     "${SRCROOT}${RKPATHPLANNINGDIR}/pptester.ui" 
#     "${SRCROOT}${RKPATHPLANNINGDIR}/ppresultview.ui"
#     "${SRCROOT}${RKPATHPLANNINGDIR}/rrtproperties.ui"
#)
#QT4_WRAP_UI( PATHPLANNING_PPTESTER_UI_HDRS ${PATHPLANNING_PPTESTER_UIS} )
#set( PATHPLANNING_PPTESTER_MOC_HDRS
#     "${SRCROOT}${RKPATHPLANNINGDIR}/pptester_impl.h"
#)
#QT4_WRAP_CPP( PATHPLANNING_PPTESTER_MOC_SRCS ${PATHPLANNING_PPTESTER_MOC_HDRS} )
#include_directories( AFTER ${CMAKE_BINARY_DIR}${RKPATHPLANNINGDIR} )
#SET (CMAKE_MODULE_PATH "/usr/share/cmake/Modules")
#INCLUDE (UsePlayerC++)

#PLAYER_ADD_PLAYERCPP_CLIENT(pptester SOURCES "${SRCROOT}${RKPATHPLANNINGDIR}/pptester_impl.cpp" ${PATHPLANNING_PPTESTER_MOC_SRCS} ${PATHPLANNING_PPTESTER_UI_HDRS} LINKFLAGS ${replaceLib} ${QT_LIBRARIES} "-lcv" "-lcvaux" ${BOOST_THREAD_LIB} ${BOOST_DATE_TIME_LIB})
#setup_custom_target(pptester "${SRCROOT}${RKPATHPLANNINGDIR}")








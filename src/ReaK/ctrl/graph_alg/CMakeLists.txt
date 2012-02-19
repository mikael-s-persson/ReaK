
set(GRAPH_ALG_SOURCES 
)

set(GRAPH_ALG_HEADERS 
  "${RKGRAPHALGDIR}/adstar_search.hpp"
  "${RKGRAPHALGDIR}/rr_tree.hpp"
  "${RKGRAPHALGDIR}/fadprm.hpp"
  "${RKGRAPHALGDIR}/probabilistic_roadmap.hpp"
  
)

#add_library(reak_graph_alg STATIC ${GRAPH_ALG_SOURCES})
#setup_custom_target(reak_graph_alg "${SRCROOT}${RKGRAPHALGDIR}")
setup_headers("${GRAPH_ALG_HEADERS}" "${RKGRAPHALGDIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

if(NOT WIN32)

  if( FreeImage_FOUND )
    add_executable(test_adstar "${SRCROOT}${RKGRAPHALGDIR}/test_adstar.cpp")
    setup_custom_target(test_adstar "${SRCROOT}${RKGRAPHALGDIR}")
    #target_link_libraries(test_adstar ${Boost_LIBRARIES})
    target_link_libraries(test_adstar ${FreeImage_LIB})
    target_link_libraries(test_adstar reak_core)
  endif()

  if( OpenCV_FOUND )
    add_executable(test_rrt "${SRCROOT}${RKGRAPHALGDIR}/test_rrt.cpp")
    setup_custom_target(test_rrt "${SRCROOT}${RKGRAPHALGDIR}")
    target_link_libraries(test_rrt ${OpenCV_LIBS} reak_core)

    add_executable(test_prm "${SRCROOT}${RKGRAPHALGDIR}/test_prm.cpp")
    setup_custom_target(test_prm "${SRCROOT}${RKGRAPHALGDIR}")
    target_link_libraries(test_prm ${OpenCV_LIBS} reak_core)

    add_executable(test_fadprm "${SRCROOT}${RKGRAPHALGDIR}/test_fadprm.cpp")
    setup_custom_target(test_fadprm "${SRCROOT}${RKGRAPHALGDIR}")
    target_link_libraries(test_fadprm ${OpenCV_LIBS} reak_core)
  endif()
endif()

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")








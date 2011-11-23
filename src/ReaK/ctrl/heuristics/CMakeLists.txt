
set(HEURISTICS_SOURCES 
)

set(HEURISTICS_HEADERS 
  "${RKHEURISTICSDIR}/goal_reach_time_heuristic.hpp"
)

#add_library(reak_heuristics STATIC ${HEURISTICS_SOURCES})
#setup_custom_target(reak_heuristics "${SRCROOT}${RKHEURISTICSDIR}")
setup_headers("${HEURISTICS_HEADERS}" "${RKHEURISTICSDIR}")

# if(NOT WIN32)
#   add_executable(test_heuristics "${SRCROOT}${RKHEURISTICSDIR}/test_heuristics.cpp")
#   setup_custom_target(test_heuristics "${SRCROOT}${RKHEURISTICSDIR}")
#   #target_link_libraries(test_adstar ${BOOST_THREAD_LIB})
#   target_link_libraries(test_heuristics "-lcv" "-lcvaux")
#   target_link_libraries(test_heuristics reak_core)
# endif()









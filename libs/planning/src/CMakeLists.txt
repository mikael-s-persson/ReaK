
set(PLANNING_SOURCES
#   "something.cpp"
)

set(PLANNERS_SOURCES 
  "fadprm_manip_planners_dynamic.cpp"
  "fadprm_manip_planners_static.cpp"
  "prm_manip_planners_dynamic.cpp"
  "prm_manip_planners_static.cpp"
#   "rrt_manip_planners_dynamic.cpp"
  "rrt_manip_planners_static.cpp"
  "rrtstar_manip_planners_dynamic.cpp"
  "rrtstar_manip_planners_static.cpp"
  "sbastar_manip_planners_dynamic.cpp"
  "sbastar_manip_planners_static.cpp"
)

# foreach(PLANNERS_DIM_COUNTER 0 1 2 3 4 5 6 7 8 9 10)
foreach(PLANNERS_DIM_COUNTER 7)
  add_library(reakobj_planners_${PLANNERS_DIM_COUNTER} OBJECT ${PLANNERS_SOURCES})
  set_target_properties(reakobj_planners_${PLANNERS_DIM_COUNTER} PROPERTIES COMPILE_DEFINITIONS "RK_PLANNERS_DIM=${PLANNERS_DIM_COUNTER}")
  set(PLANNERS_OBJECT_FILES ${PLANNERS_OBJECT_FILES} $<TARGET_OBJECTS:reakobj_planners_${PLANNERS_DIM_COUNTER}>)
endforeach(PLANNERS_DIM_COUNTER)

add_library(reak_planning STATIC ${PLANNERS_OBJECT_FILES}) #${PLANNING_SOURCES})
ReaK_setup_static_library(reak_planning)
target_link_libraries(reak_planning reak_topologies reak_mbd)


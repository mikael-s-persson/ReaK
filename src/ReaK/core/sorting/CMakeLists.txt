
#set(SORTING_SOURCES "${SRCROOT}${RKROOTFINDERSDIR}/.cpp")

set(SORTING_HEADERS 
  "${RKSORTINGDIR}/bisection_method.hpp"
  "${RKSORTINGDIR}/broyden_method.hpp"
  "${RKSORTINGDIR}/newton_raphson_method.hpp"
  "${RKSORTINGDIR}/secant_method.hpp"
)

#add_library(reak_sorting STATIC ${SORTING_SOURCES})
#setup_custom_target(reak_sorting "${SRCROOT}${RKSORTINGDIR}")

setup_headers("${SORTING_HEADERS}" "${RKSORTINGDIR}")

add_executable(test_sorting_perf "${SRCROOT}${RKSORTINGDIR}/test_sorting_perf.cpp")
setup_custom_target(test_sorting_perf "${SRCROOT}${RKSORTINGDIR}")
target_link_libraries(test_sorting_perf reak_rtti ${EXTRA_SYSTEM_LIBS})








#set(SORTING_SOURCES "${SRCROOT}${RKROOTFINDERSDIR}/.cpp")

set(SORTING_HEADERS 
  "${RKSORTINGDIR}/bubble_sort.hpp"
  "${RKSORTINGDIR}/comb_sort.hpp"
  "${RKSORTINGDIR}/heap_sort.hpp"
  "${RKSORTINGDIR}/insertion_sort.hpp"
  "${RKSORTINGDIR}/intro_sort.hpp"
  "${RKSORTINGDIR}/merge_sort.hpp"
  "${RKSORTINGDIR}/quick_sort.hpp"
  "${RKSORTINGDIR}/selection_sort.hpp"
  "${RKSORTINGDIR}/shell_sort.hpp"
)

#add_library(reak_sorting STATIC ${SORTING_SOURCES})
#setup_custom_target(reak_sorting)

setup_headers("${SORTING_HEADERS}" "${RKSORTINGDIR}")


if(Boost_CHRONO_FOUND)

add_executable(test_sorting_perf "${SRCROOT}${RKSORTINGDIR}/test_sorting_perf.cpp")
setup_custom_target(test_sorting_perf)
target_link_libraries(test_sorting_perf reak_core)

endif()

add_executable(unit_test_sorting "${SRCROOT}${RKSORTINGDIR}/unit_test_sorting.cpp")
setup_custom_test_program(unit_test_sorting)
target_link_libraries(unit_test_sorting reak_core)







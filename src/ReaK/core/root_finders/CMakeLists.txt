
#set(ROOTFINDERS_SOURCES "${SRCROOT}${RKROOTFINDERSDIR}/.cpp")

set(ROOTFINDERS_HEADERS 
  "${RKROOTFINDERSDIR}/bisection_method.hpp"
  "${RKROOTFINDERSDIR}/broyden_method.hpp"
  "${RKROOTFINDERSDIR}/newton_raphson_method.hpp"
  "${RKROOTFINDERSDIR}/secant_method.hpp"
)

setup_headers("${ROOTFINDERS_HEADERS}" "${RKROOTFINDERSDIR}")

add_executable(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}/test_root_finders_perf.cpp")
setup_custom_target(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}")
target_link_libraries(test_root_finders_perf reak_core)

add_executable(unit_test_root_finders "${SRCROOT}${RKROOTFINDERSDIR}/unit_test_root_finders.cpp")
setup_custom_test_program(unit_test_root_finders "${SRCROOT}${RKROOTFINDERSDIR}")
target_link_libraries(unit_test_root_finders reak_core)







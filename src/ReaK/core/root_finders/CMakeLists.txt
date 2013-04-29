
#set(ROOTFINDERS_SOURCES "${SRCROOT}${RKROOTFINDERSDIR}/.cpp")

set(ROOTFINDERS_HEADERS 
  "${RKROOTFINDERSDIR}/bisection_method.hpp"
  "${RKROOTFINDERSDIR}/broyden_method.hpp"
  "${RKROOTFINDERSDIR}/newton_raphson_method.hpp"
  "${RKROOTFINDERSDIR}/secant_method.hpp"
)

#add_library(reak_lin_alg STATIC ${ROOTFINDERS_SOURCES})
#setup_custom_target(reak_lin_alg "${SRCROOT}${RKROOTFINDERSDIR}")

setup_headers("${ROOTFINDERS_HEADERS}" "${RKROOTFINDERSDIR}")

add_executable(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}/test_root_finders_perf.cpp")
setup_custom_target(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}")
target_link_libraries(test_root_finders_perf reak_rtti ${EXTRA_SYSTEM_LIBS})

add_executable(unit_test_root_finders "${SRCROOT}${RKROOTFINDERSDIR}/unit_test_root_finders.cpp")
setup_custom_test_program(unit_test_root_finders "${SRCROOT}${RKROOTFINDERSDIR}")
target_link_libraries(unit_test_root_finders reak_rtti)
target_link_libraries(unit_test_root_finders ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})







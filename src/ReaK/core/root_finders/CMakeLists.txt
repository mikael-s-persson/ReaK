
#set(ROOTFINDERS_SOURCES "${SRCROOT}${RKROOTFINDERSDIR}/.cpp")

set(ROOTFINDERS_HEADERS 
                 "${RKROOTFINDERSDIR}/bisection_method.hpp"
                 "${RKROOTFINDERSDIR}/secant_method.hpp"
                 "${RKROOTFINDERSDIR}/newton_raphson_method.hpp"
                 "${RKROOTFINDERSDIR}/broyden_method.hpp"
)

#add_library(reak_lin_alg STATIC ${ROOTFINDERS_SOURCES})
#setup_custom_target(reak_lin_alg "${SRCROOT}${RKROOTFINDERSDIR}")

setup_headers("${ROOTFINDERS_HEADERS}" "${RKROOTFINDERSDIR}")

add_executable(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}/test_root_finders_perf.cpp")
setup_custom_target(test_root_finders_perf "${SRCROOT}${RKROOTFINDERSDIR}")
target_link_libraries(test_root_finders_perf reak_rtti)








set(OPTIM_SOURCES "${SRCROOT}${RKOPTIMDIR}/line_search.cpp")

set(OPTIM_HEADERS "${RKOPTIMDIR}/line_search.hpp")
set(OPTIM_HEADERS "${RKOPTIMDIR}/function_types.hpp")

add_library(reak_optimization STATIC ${OPTIM_SOURCES})
setup_custom_target(reak_optimization "${SRCROOT}${RKOPTIMDIR}")

setup_headers("${OPTIM_HEADERS}" "${RKOPTIMDIR}")

add_executable(test_optim "${SRCROOT}${RKOPTIMDIR}/test.cpp")
setup_custom_target(test_optim "${SRCROOT}${RKOPTIMDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

target_link_libraries(test_optim reak_optimization reak_rtti)
target_link_libraries(test_optim ${BOOST_THREAD_LIB})





set(OPTIM_SOURCES 
)

set(OPTIM_HEADERS 
  "${RKOPTIMDIR}/line_search.hpp"
  "${RKOPTIMDIR}/function_types.hpp"
)

setup_headers("${OPTIM_HEADERS}" "${RKOPTIMDIR}")

add_executable(test_optim "${SRCROOT}${RKOPTIMDIR}/test.cpp")
setup_custom_target(test_optim "${SRCROOT}${RKOPTIMDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

target_link_libraries(test_optim reak_rtti)
target_link_libraries(test_optim ${BOOST_THREAD_LIB})




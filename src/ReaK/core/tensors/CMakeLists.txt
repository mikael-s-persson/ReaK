
# set(TENSORS_SOURCES "${SRCROOT}${RKTENSORSDIR}/.cpp")

set(TENSORS_HEADERS 
                 "${RKTENSORSDIR}/tensor_alg.hpp"
                 "${RKTENSORSDIR}/tensor_alg_general.hpp"
                 "${RKTENSORSDIR}/tensor_concepts.hpp"
                 "${RKTENSORSDIR}/tensor_traits.hpp"
)

# add_library(reak_tensors STATIC ${TENSORS_SOURCES})
# setup_custom_target(reak_tensors "${SRCROOT}${RKTENSORSDIR}")

setup_headers("${TENSORS_HEADERS}" "${RKTENSORSDIR}")

# add_executable(test_tensor_alg "${SRCROOT}${RKTENSORSDIR}/test_tensor_alg.cpp")
# setup_custom_target(test_tensor_alg "${SRCROOT}${RKTENSORSDIR}")
# target_link_libraries(test_tensor_alg reak_rtti)









# set(TENSORS_SOURCES "${SRCROOT}${RKTENSORSDIR}/.cpp")

set(TENSORS_HEADERS 
                 "${RKTENSORSDIR}/tensor_alg.hpp"
                 "${RKTENSORSDIR}/tensor_alg_general.hpp"
                 "${RKTENSORSDIR}/tensor_concepts.hpp"
                 "${RKTENSORSDIR}/tensor_traits.hpp"
)

setup_headers("${TENSORS_HEADERS}" "${RKTENSORSDIR}")








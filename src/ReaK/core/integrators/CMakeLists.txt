
set(INTEGRATORS_SOURCES )

set(INTEGRATORS_HEADERS 
  "${RKINTEGRATORSDIR}/fixed_step_integrators.hpp"
  "${RKINTEGRATORSDIR}/integration_exceptions.hpp"
  "${RKINTEGRATORSDIR}/integrator.hpp"
  "${RKINTEGRATORSDIR}/pred_corr_integrators.hpp"
  "${RKINTEGRATORSDIR}/variable_step_integrators.hpp"
)

#add_library(reak_integrators STATIC ${INTEGRATORS_SOURCES})
#setup_custom_target(reak_integrators)

setup_headers("${INTEGRATORS_HEADERS}" "${RKINTEGRATORSDIR}")

add_executable(unit_test_integrators "${SRCROOT}${RKINTEGRATORSDIR}/unit_test_integrators.cpp")
setup_custom_test_program(unit_test_integrators)
target_link_libraries(unit_test_integrators reak_core)
get_target_property(UNIT_TEST_INTEGRATORS_LOCATION unit_test_integrators LOCATION)
get_filename_component(UNIT_TEST_INTEGRATORS_PATH ${UNIT_TEST_INTEGRATORS_LOCATION} PATH)
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     make_directory "${UNIT_TEST_INTEGRATORS_PATH}/integ_records")
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     copy "${SRCROOT}${RKINTEGRATORSDIR}/integ_records/hires.pbuf" "${UNIT_TEST_INTEGRATORS_PATH}/integ_records/")
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     copy "${SRCROOT}${RKINTEGRATORSDIR}/integ_records/orego.pbuf" "${UNIT_TEST_INTEGRATORS_PATH}/integ_records/")
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     copy "${SRCROOT}${RKINTEGRATORSDIR}/integ_records/ringmod.pbuf" "${UNIT_TEST_INTEGRATORS_PATH}/integ_records/")
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     copy "${SRCROOT}${RKINTEGRATORSDIR}/integ_records/vanderpolmod.pbuf" "${UNIT_TEST_INTEGRATORS_PATH}/integ_records/")
add_custom_command(TARGET unit_test_integrators PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E 
                     copy "${SRCROOT}${RKINTEGRATORSDIR}/integ_records/vanderpol.pbuf" "${UNIT_TEST_INTEGRATORS_PATH}/integ_records/")





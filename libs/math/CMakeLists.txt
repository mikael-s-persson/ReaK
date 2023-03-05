# project(ReaK_math)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../core")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../core/include")
else()
  message(FATAL_ERROR "The ReaK.math module requires the presence of the ReaK.core module!")
endif()

set(REAK_CURRENT_MODULE math)
set(REAK_CURRENT_MODULE_TITLE "ReaK Math Library Components")
set(REAK_CURRENT_MODULE_DESCRIPTION "An advanced numerical analysis tool-set for linear algebra, matrix numerical methods, kinetostatics, and numerical integration.")

ReaK_add_current_module()




# project(ReaK_core)

find_package(BGLExtra REQUIRED)
if(BGLEXTRA_FOUND)
  message(STATUS "BGL-Extra was detected on this system: '${BGLEXTRA_INCLUDE_DIR}'")
  include_directories(SYSTEM ${BGLEXTRA_INCLUDE_DIR})
endif()

set(REAK_CURRENT_MODULE core)
set(REAK_CURRENT_MODULE_TITLE "ReaK Core Library Components")
set(REAK_CURRENT_MODULE_DESCRIPTION "Core software architectural elements of the ReaK library.")

ReaK_add_current_module()


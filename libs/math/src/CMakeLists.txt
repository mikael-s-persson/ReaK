
add_library(reak_math STATIC 
  "kinetostatics_externs.cpp"
  "mat_alg_externs.cpp"
  "mat_num_externs.cpp"
  "vect_alg_externs.cpp"
)
ReaK_setup_static_library(reak_math)
target_link_libraries(reak_math reak_core ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})


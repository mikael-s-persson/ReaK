

add_library(reak_py SHARED 
  "py_base.cpp"
  "py_kinetostatics.cpp"
  "py_kte_models.cpp"
  "py_mbd_kte.cpp"
  "py_vect_alg.cpp"
  "reak_py_bindings.cpp"
)
ReaK_setup_shared_library(reak_py)
target_link_libraries(reak_py reak_mbd reak_math ${PYTHON_LIBRARIES} ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})


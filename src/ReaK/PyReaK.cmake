

find_package(PythonLibs)

if(PYTHONLIBS_FOUND)

set(PYREAK_SOURCES 
  "${SRCROOT}${REAKDIR}/reak_py_bindings.cpp"
  "${SRCROOT}${RKCOREDIR}/base/py_base.cpp"
  "${SRCROOT}${RKCOREDIR}/lin_alg/py_vect_alg.cpp"
  "${SRCROOT}${RKCOREDIR}/kinetostatics/py_kinetostatics.cpp"
)


include_directories(BEFORE ${PYTHON_INCLUDE_DIRS})

add_library(reak_py SHARED ${PYREAK_SOURCES})
setup_custom_target(reak_py "${SRCROOT}${RKRECORDERSDIR}")
target_link_libraries(reak_py ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} reak_core)

endif()

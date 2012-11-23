

find_package(PythonLibs)

if(PYTHONLIBS_FOUND)

  set(PYREAK_SOURCES 
    "${SRCROOT}${REAKDIR}/reak_py_bindings.cpp"
    "${SRCROOT}${RKCOREDIR}/base/py_base.cpp"
    "${SRCROOT}${RKCOREDIR}/lin_alg/py_vect_alg.cpp"
    "${SRCROOT}${RKCOREDIR}/kinetostatics/py_kinetostatics.cpp"
    "${SRCROOT}${RKCTRLDIR}/mbd_kte/py_mbd_kte.cpp"
    "${SRCROOT}${RKCTRLDIR}/kte_models/py_kte_models.cpp"
  )
  
  set(original_Boost_FOUND ${Boost_FOUND})
  find_package(Boost 1.42 COMPONENTS python)
  if(Boost_FOUND)
    
    include_directories(BEFORE ${PYTHON_INCLUDE_DIRS})
    
    add_library(reak_py SHARED ${PYREAK_SOURCES})
    setup_custom_target(reak_py "${SRCROOT}${RKRECORDERSDIR}")
    target_link_libraries(reak_py ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} reak_mbd_kte reak_core)
  
  endif()
  set(Boost_FOUND ${original_Boost_FOUND})
  
endif()

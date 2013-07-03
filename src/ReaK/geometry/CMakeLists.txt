# These are the main sub-directories
set(RKSHAPESDIR "${RKGEOMETRYDIR}/shapes")
set(RKPROXIMITYDIR "${RKGEOMETRYDIR}/proximity")

set(RKGEOMETRYSUBDIRS ${RKSHAPESDIR} ${RKPROXIMITYDIR})

configure_doxyfile(ReaKgeom 
                   "\"ReaK Library - Geometry library\"" 
                   "\"A tool-set for geometric representations, proximity queries and rendering for 2D and 3D objects.\"" 
                   "${DOXTOSRCDIR}${RKSHAPESDIR} ${DOXTOSRCDIR}${RKPROXIMITYDIR}")
add_doxygen_target(ReaKgeom)


# Add the sub-directories to cmake process.
foreach(CURRENT_RKGEOMETRYSUBDIR ${RKGEOMETRYSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKGEOMETRYSUBDIR}")
endforeach(CURRENT_RKGEOMETRYSUBDIR)



add_library(reak_mbd STATIC $<TARGET_OBJECTS:reakobj_geom> $<TARGET_OBJECTS:reakobj_proximity> $<TARGET_OBJECTS:reakobj_mbd_kte> $<TARGET_OBJECTS:reakobj_kte_models>)
setup_custom_target(reak_mbd "${SRCROOT}${RKCTRLDIR}")
target_link_libraries(reak_mbd reak_core)



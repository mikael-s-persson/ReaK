# These are the main sub-directories
set(RKSHAPESDIR "${RKGEOMETRYDIR}/shapes")

set(RKGEOMETRYSUBDIRS ${RKSHAPESDIR})

configure_doxyfile(ReaKgeom 
                   "\"ReaK Library - Geometry library\"" 
                   "\"A tool-set for geometric representations, collision calculations and rendering for 2D and 3D objects.\"" 
                   "${DOXTOSRCDIR}${RKSHAPESDIR}")
add_doxygen_target(ReaKgeom)


# Add the sub-directories to cmake process.
foreach(CURRENT_RKGEOMETRYSUBDIR ${RKGEOMETRYSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKGEOMETRYSUBDIR}")
endforeach(CURRENT_RKGEOMETRYSUBDIR)





# These are the main sub-directories
set(RKGRAPHALGDIR "${RKCTRLDIR}/graph_alg")
set(RKPATHPLANNINGDIR "${RKCTRLDIR}/path_planning")
set(RKCTRLSYSDIR "${RKCTRLDIR}/ctrl_sys")
set(RKMBDKTEDIR "${RKCTRLDIR}/mbd_kte")
set(RKINTERPOLATIONDIR "${RKCTRLDIR}/interpolation")

set(RKCTRLSUBDIRS ${RKMBDKTEDIR} ${RKGRAPHALGDIR} ${RKPATHPLANNINGDIR} ${RKCTRLSYSDIR} ${RKINTERPOLATIONDIR})

# Add the sub-directories to cmake process.
foreach(CURRENT_RKCTRLSUBDIR ${RKCTRLSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKCTRLSUBDIR}")
endforeach(CURRENT_RKCTRLSUBDIR)





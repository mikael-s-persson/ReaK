# These are the main sub-directories
set(RKGRAPHALGDIR "${RKCTRLDIR}/graph_alg")
set(RKPATHPLANNINGDIR "${RKCTRLDIR}/path_planning")
set(RKCTRLSYSDIR "${RKCTRLDIR}/ctrl_sys")
set(RKMBDKTEDIR "${RKCTRLDIR}/mbd_kte")
set(RKKTEMODELSDIR "${RKCTRLDIR}/kte_models")
set(RKINTERPOLATIONDIR "${RKCTRLDIR}/interpolation")
set(RKSSSYSTEMSDIR "${RKCTRLDIR}/ss_systems")
set(RKSYSINTEGRATORSDIR "${RKCTRLDIR}/sys_integrators")
set(RKTOPOLOGIESDIR "${RKCTRLDIR}/topologies")
set(RKHEURISTICSDIR "${RKCTRLDIR}/heuristics")

set(RKCTRLSUBDIRS ${RKMBDKTEDIR} ${RKKTEMODELSDIR} ${RKGRAPHALGDIR} ${RKPATHPLANNINGDIR} ${RKCTRLSYSDIR} ${RKINTERPOLATIONDIR} ${RKTOPOLOGIESDIR} ${RKSSSYSTEMSDIR} ${RKHEURISTICSDIR} ${RKSYSINTEGRATORSDIR})

configure_doxyfile(ReaKctrl_sys 
                   "\"ReaK Library - Control and Estimation libraries\"" 
                   "\"Numerical library for control and estimation on dynamic systems.\"" 
                   "${DOXTOSRCDIR}${RKCTRLSYSDIR} ${DOXTOSRCDIR}${RKSSSYSTEMSDIR} ${DOXTOSRCDIR}${RKSYSINTEGRATORSDIR}")
add_doxygen_target(ReaKctrl_sys)

configure_doxyfile(ReaKpp 
                   "\"ReaK Library - Path-planning libraries\"" 
                   "\"A collection of path-planning methods, graph-search methods, interpolation and topology-based numerical analysis.\"" 
                   "${DOXTOSRCDIR}${RKGRAPHALGDIR} ${DOXTOSRCDIR}${RKPATHPLANNINGDIR} ${DOXTOSRCDIR}${RKINTERPOLATIONDIR} ${DOXTOSRCDIR}${RKTOPOLOGIESDIR} ${DOXTOSRCDIR}${RKHEURISTICSDIR}")
add_doxygen_target(ReaKpp)

configure_doxyfile(ReaKmbd_kte 
                   "\"ReaK Library - Multi-body dynamics library\"" 
                   "\"A complete tool-set for multi-body dynamics modeling based on Kinetostatic Transmission Elements (KTE).\"" 
                   "${DOXTOSRCDIR}${RKMBDKTEDIR} ${DOXTOSRCDIR}${RKKTEMODELSDIR}")
add_doxygen_target(ReaKmbd_kte)


# Add the sub-directories to cmake process.
foreach(CURRENT_RKCTRLSUBDIR ${RKCTRLSUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKCTRLSUBDIR}")
endforeach(CURRENT_RKCTRLSUBDIR)


add_library(reak_kte STATIC $<TARGET_OBJECTS:reakobj_mbd_kte> $<TARGET_OBJECTS:reakobj_kte_models>)
setup_custom_target(reak_kte "${SRCROOT}${RKCTRLDIR}")
target_link_libraries(reak_kte reak_core)






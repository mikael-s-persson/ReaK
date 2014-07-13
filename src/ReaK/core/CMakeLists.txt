# These are the main sub-directories
set(RKLINALGDIR "${RKCOREDIR}/lin_alg")
set(RKKINETOSTATICSDIR "${RKCOREDIR}/kinetostatics")
set(RKBASEDIR "${RKCOREDIR}/base")
set(RKOPTIMDIR "${RKCOREDIR}/optimization")
set(RKRTTIDIR "${RKCOREDIR}/rtti")
set(RKSERIALIZATIONDIR "${RKCOREDIR}/serialization")
set(RKRECORDERSDIR "${RKCOREDIR}/recorders")
set(RKINTEGRATORSDIR "${RKCOREDIR}/integrators")
set(RKROOTFINDERSDIR "${RKCOREDIR}/root_finders")
set(RKTENSORSDIR "${RKCOREDIR}/tensors")
set(RKSORTINGDIR "${RKCOREDIR}/sorting")

set(RKCORESUBDIRS 
  ${RKBASEDIR} 
  ${RKLINALGDIR} 
  ${RKKINETOSTATICSDIR} 
  ${RKRTTIDIR} 
  ${RKSERIALIZATIONDIR} 
  ${RKRECORDERSDIR}
  ${RKINTEGRATORSDIR} 
  ${RKOPTIMDIR}
  ${RKROOTFINDERSDIR}
  ${RKTENSORSDIR}
  ${RKSORTINGDIR}
)

configure_doxyfile(ReaKallcore 
                   "\"ReaK Library - All core libraries\"" 
                   "\"Core libraries for software architectural elements and advanced numerical analysis software.\"" 
                   "${RKBASEDIR};${RKLINALGDIR};${RKKINETOSTATICSDIR};${RKRTTIDIR};${RKSERIALIZATIONDIR};${RKRECORDERSDIR};${RKINTEGRATORSDIR};${RKROOTFINDERSDIR};${RKOPTIMDIR};${RKTENSORSDIR};${RKSORTINGDIR}")
add_doxygen_target(ReaKallcore)

configure_doxyfile(ReaKcore 
                   "\"ReaK Library - Core libraries\"" 
                   "\"Core software architectural elements of the ReaK library.\"" 
                   "${RKBASEDIR};${RKRTTIDIR};${RKSERIALIZATIONDIR};${RKRECORDERSDIR}")
add_doxygen_target(ReaKcore)

configure_doxyfile(ReaKmath 
                   "\"ReaK Library - Math libraries\"" 
                   "\"An advanced numerical analysis tool-set for linear algebra, matrix numerical methods, kinetostatics, and numerical integration.\"" 
                   "${RKBASEDIR};${RKLINALGDIR};${RKKINETOSTATICSDIR};${RKINTEGRATORSDIR};${RKROOTFINDERSDIR};${RKOPTIMDIR};${RKTENSORSDIR};${RKSORTINGDIR}")
add_doxygen_target(ReaKmath)

# Add the sub-directories to cmake process.
foreach(CURRENT_RKCORESUBDIR ${RKCORESUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKCORESUBDIR}")
endforeach(CURRENT_RKCORESUBDIR)

add_library(reak_core STATIC $<TARGET_OBJECTS:reakobj_rtti> $<TARGET_OBJECTS:reakobj_recorders> $<TARGET_OBJECTS:reakobj_serialization> $<TARGET_OBJECTS:reakobj_lin_alg> $<TARGET_OBJECTS:reakobj_kinetostatics>)
setup_custom_target(reak_core)
target_link_libraries(reak_core ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})


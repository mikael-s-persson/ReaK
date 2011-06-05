# These are the main sub-directories
set(RKMATHDIR "${RKCOREDIR}/math")
set(RKBASEDIR "${RKCOREDIR}/base")
set(RKOPTIMDIR "${RKCOREDIR}/optimization")
set(RKRTTIDIR "${RKCOREDIR}/rtti")
set(RKSERIALIZATIONDIR "${RKCOREDIR}/serialization")
set(RKRECORDERSDIR "${RKCOREDIR}/recorders")
set(RKINTEGRATORSDIR "${RKCOREDIR}/integrators")

set(RKCORESUBDIRS 
  ${RKBASEDIR} 
  ${RKMATHDIR} 
  ${RKRTTIDIR} 
  ${RKSERIALIZATIONDIR} 
  ${RKRECORDERSDIR} 
  ${RKCMDHDLDIR} 
  ${RKINTEGRATORSDIR} 
  ${RKOPTIMDIR}
)

# Add the sub-directories to cmake process.
foreach(CURRENT_RKCORESUBDIR ${RKCORESUBDIRS})
add_subdirectory("${SRCROOT}${CURRENT_RKCORESUBDIR}")
endforeach(CURRENT_RKCORESUBDIR)

get_target_property(RKRTTISOURCES reak_rtti SOURCES)
get_target_property(RKRECORDERSSOURCES reak_recorders SOURCES)
get_target_property(RKSERIALIZATIONSOURCES reak_serialization SOURCES)
add_library(reak_core STATIC ${RKRTTISOURCES} ${RKRECORDERSSOURCES} ${RKSERIALIZATIONSOURCES})
#target_link_libraries(reak_core reak_rtti reak_serialization reak_recorders)
setup_custom_target(reak_core "${SRCROOT}${RKCOREDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

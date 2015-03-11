

set(REAK_CORE_SOURCES 
  "archiver_factory.cpp"
  "ascii_recorder.cpp"
  "bin_archiver.cpp"
  "bin_recorder.cpp"
  "data_record.cpp"
  "data_record_options.cpp"
  "network_recorder.cpp"
  "objtree_archiver.cpp"
  "protobuf_archiver.cpp"
  "scheme_builder.cpp"
  "so_type.cpp"
  "so_type_repo.cpp"
  "vector_recorder.cpp"
  "xml_archiver.cpp"
)

if(_COMPILER_SUPPORTS_CXX11)
set(REAK_CORE_SOURCES ${REAK_CORE_SOURCES} 
  "disc_headers.cpp"
  "remote_function.cpp"
  "rpc_exceptions.cpp"
  "rpc_headers.cpp"
  "rpc_server.cpp"
  "rpc_version.cpp"
)
endif()

add_library(reak_core STATIC ${REAK_CORE_SOURCES})
ReaK_setup_static_library(reak_core)
target_link_libraries(reak_core ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})


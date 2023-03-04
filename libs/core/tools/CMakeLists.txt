

add_executable(convert_datastream "convert_datastream.cpp")
ReaK_setup_tool_program(convert_datastream)
target_link_libraries(convert_datastream reak_core)
target_link_libraries(convert_datastream absl::flags absl::flags_parse)


load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
  name = "boost_graph_extras",
  srcs = [],
  hdrs =
    glob([
      "include/boost/**",
    ]),
  strip_include_prefix = "include/",
  visibility = ["//visibility:public"],
  deps = [
    "@boost//:graph",
  ],
)

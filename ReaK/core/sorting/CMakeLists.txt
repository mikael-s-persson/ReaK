

reak_cc_library(
  NAME
    sorting
  HDRS
    "bubble_sort.h"
    "comb_sort.h"
    "heap_sort.h"
    "insertion_sort.h"
    "intro_sort.h"
    "merge_sort.h"
    "quick_sort.h"
    "selection_sort.h"
    "shell_sort.h"
  PUBLIC
)

reak_cc_test(
  NAME
    sorting_test
  SRCS
    "sorting_test.cc"
  DEPS
    ReaK::sorting
    GTest::gtest_main
)

reak_cc_binary(
  NAME
    sorting_perf
  SRCS
    "sorting_perf.cc"
  DEPS
    ReaK::sorting
)

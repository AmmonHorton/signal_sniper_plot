load("@rules_cc//cc:defs.bzl", "cc_library")

cc_library(
    name = "x11",
    hdrs = glob([
        "include/X11/**/*.h",
        "include/X11/extensions/*.h",
    ], allow_empty = True),
    includes = ["include"],
    linkopts = ["-lX11"],
    visibility = ["//visibility:public"],
)
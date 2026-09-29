"""Where the RPM and deb put the library, its pkg-config file and the headers."""

load("@rules_pkg//pkg:mappings.bzl", "pkg_attributes", "pkg_files", "pkg_mklink")

def install_files(name, libdir, library, lib, soversion, summary, version):
    """Library, dev symlink, pkg-config file and headers, with the library in `libdir`."""
    pkg_files(
        name = name + "_lib",
        srcs = [library],
        prefix = libdir,
        attributes = pkg_attributes(mode = "0755"),
    )
    pkg_mklink(
        name = name + "_devlink",
        link_name = libdir + "/" + lib,
        target = lib + "." + soversion,
    )
    native.genrule(
        name = name + "_pc_file",
        outs = [name + "/signal_sniper_plot.pc"],
        cmd = """cat > $@ <<'PC'
prefix=/usr
libdir={libdir}
includedir=$${{prefix}}/include/signal_sniper_plot

Name: signal_sniper_plot
Description: {summary}
Version: {version}
Libs: -L$${{libdir}} -lsignal_sniper_plot
Cflags: -I$${{includedir}}
PC""".format(libdir = libdir, summary = summary, version = version),
    )
    pkg_files(
        name = name + "_pc",
        srcs = [name + "_pc_file"],
        prefix = libdir + "/pkgconfig",
    )
    pkg_files(
        name = name + "_headers",
        srcs = ["//ssp:public_headers"],
        prefix = "/usr/include/signal_sniper_plot/ssp",
    )
    return [name + "_lib", name + "_devlink", name + "_pc", name + "_headers"]

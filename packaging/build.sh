#!/usr/bin/env bash
# Build release artifacts, each inside the oldest distribution it should install on, so it
# only needs that distribution's glibc (a package built on a newer host needs the newer glibc).
#
#   packaging/build.sh wheel   # manylinux_2_28 (AlmaLinux 8): pip, any Linux with glibc >= 2.28
#   packaging/build.sh rpm     # Rocky Linux 8: RHEL / CentOS Stream / Alma / Rocky 8 and newer
#   packaging/build.sh deb     # Ubuntu 20.04 and newer
#   packaging/build.sh all
#
# Needs docker. Artifacts land in dist/. Each image keeps its Bazel cache in a docker volume,
# so re-runs are quick. CI runs this same script.
set -euo pipefail

usage() { sed -n '2,11p' "$0"; exit 2; }
[[ $# -eq 1 ]] || usage
target=$1

if [[ $target == all ]]; then
    for t in wheel rpm deb; do "$0" "$t"; done
    exit 0
fi

case $target in
    wheel) image=quay.io/pypa/manylinux_2_28_x86_64 ;;
    rpm) image=rockylinux:8 ;;
    deb) image=ubuntu:20.04 ;;
    *) usage ;;
esac

if [[ "${SSP_IN_CONTAINER:-}" != 1 ]]; then
    root=$(cd "$(dirname "$0")/.." && pwd)
    exec docker run --rm \
        -e SSP_IN_CONTAINER=1 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
        -v "$root":/src -w /src \
        -v "ssp-cache-$target":/root/.cache \
        "$image" /src/packaging/build.sh "$target"
fi

# ── Inside the container ─────────────────────────────────────────────────────
case $target in
    wheel)
        dnf install -y -q libX11-devel
        ;;
    rpm)
        # Rocky 8's own GCC 8 is too old for the code; gcc-toolset builds binaries that still
        # run on the base system's libstdc++.
        dnf install -y -q gcc-toolset-12 libX11-devel rpm-build python3
        source /opt/rh/gcc-toolset-12/enable
        ;;
    deb)
        export DEBIAN_FRONTEND=noninteractive
        apt-get update -qq
        apt-get install -y -qq g++ libx11-dev curl ca-certificates python3 >/dev/null
        ;;
esac

BAZELISK=/root/.cache/bin/bazelisk
if [[ ! -x $BAZELISK ]]; then
    mkdir -p "$(dirname "$BAZELISK")"
    curl -fsSL -o "$BAZELISK" \
        https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
    chmod +x "$BAZELISK"
fi
bazel() { "$BAZELISK" "$@"; }

# --symlink_prefix=/ : leave the host's bazel-* convenience symlinks alone.
bazel build -c opt --symlink_prefix=/ "//packaging:$target"
out=$(bazel info -c opt execution_root 2>/dev/null)/$(bazel cquery -c opt --output=files "//packaging:$target" 2>/dev/null)

mkdir -p dist
if [[ $target == wheel ]]; then
    auditwheel repair --plat manylinux_2_28_x86_64 --wheel-dir dist "$out"
else
    install -m 0644 "$out" dist/
fi
chown -R "$HOST_UID:$HOST_GID" dist
ls -l dist

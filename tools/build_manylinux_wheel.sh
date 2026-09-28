#!/usr/bin/env bash
# Build the release wheel inside PyPA's manylinux_2_28 image, so it installs on any x86_64
# Linux with glibc >= 2.28 (RHEL / CentOS Stream / Alma / Rocky 8+, Ubuntu 20.04+, ...).
# Building on a newer host would tie the wheel to that host's glibc instead.
#
#   tools/build_manylinux_wheel.sh        # needs docker; writes dist/*.whl
#
# CI runs this same script. Bazel's cache lives in a docker volume, so re-runs are quick.
set -euo pipefail

IMAGE=quay.io/pypa/manylinux_2_28_x86_64
PLAT=manylinux_2_28_x86_64

if [[ "${SSP_IN_CONTAINER:-}" != 1 ]]; then
    root=$(cd "$(dirname "$0")/.." && pwd)
    exec docker run --rm \
        -e SSP_IN_CONTAINER=1 -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
        -v "$root":/src -w /src \
        -v ssp-manylinux-cache:/root/.cache \
        "$IMAGE" /src/tools/build_manylinux_wheel.sh
fi

# ── Inside the container ─────────────────────────────────────────────────────
dnf install -y -q libX11-devel
BAZELISK=/root/.cache/bin/bazelisk
if [[ ! -x $BAZELISK ]]; then
    mkdir -p "$(dirname "$BAZELISK")"
    curl -fsSL -o "$BAZELISK" \
        https://github.com/bazelbuild/bazelisk/releases/latest/download/bazelisk-linux-amd64
    chmod +x "$BAZELISK"
fi
bazel() { "$BAZELISK" "$@"; }

# --symlink_prefix=/ : don't replace the host's bazel-* convenience symlinks with container paths.
bazel build -c opt --symlink_prefix=/ //:signal_sniper_plot_wheel
wheel=$(bazel cquery -c opt --output=files //:signal_sniper_plot_wheel 2>/dev/null)
wheel=$(bazel info -c opt execution_root 2>/dev/null)/$wheel

mkdir -p dist
auditwheel repair --plat "$PLAT" --wheel-dir dist "$wheel"
chown -R "$HOST_UID:$HOST_GID" dist
ls -l dist

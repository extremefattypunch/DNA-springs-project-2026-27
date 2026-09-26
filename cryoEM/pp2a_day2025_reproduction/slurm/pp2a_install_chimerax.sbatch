#!/bin/bash -l
# Plan section 5A - ChimeraX 1.12 on a Rocky 8 cluster, without root.  (Recipe that worked on 2026-09-26.)
#
# UCSF ships ChimeraX 1.12 only for RHEL 9/10 and Ubuntu; Cannon nodes are Rocky 8 (glibc 2.28), and
# the "centos-9" build is made on CentOS Stream 9 (needs OpenSSL 3.3+, so a plain Rocky 9 image fails).
#  1. licence handshake + download of the el9 rpm (free for academic use), unpacked with rpm2cpio into
#     holylabs software/chimerax (1.8 GB);
#  2. a CentOS Stream 9 Singularity image (quay.io/centos/centos:stream9, 60 MB) supplies glibc/OpenSSL;
#  3. the X11/GL/OSMesa/Qt runtime libraries the minimal image lacks are resolved with `dnf repoquery
#     --location` inside that image, downloaded as rpms and unpacked into software/chimerax/sysroot
#     (~130 MB); the launcher adds sysroot/usr/lib64 to LD_LIBRARY_PATH;
#  4. smoke test: headless render of 9MIP to PNG via --nogui --offscreen (OSMesa).
# Every rpm is deleted after unpacking (hekstra-cluster rule 3).
#SBATCH -J pp2a-install-chimerax
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-01:30
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
CX=$SW/chimerax; SYS=$CX/sysroot; SIF=$SW/containers/centos-stream9.sif
CXBIN=$CX/usr/libexec/UCSF-ChimeraX/bin/ChimeraX
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP" "$SW/containers"; trap 'rm -rf "$TMP"' EXIT
export SINGULARITY_CACHEDIR=$TMP/sifcache SINGULARITY_TMPDIR=$TMP
VER=1.12; RPM="$VER/centos-9/ucsf-chimerax-$VER-1.el9.x86_64.rpm"

# 1. ChimeraX rpm
if [ ! -x "$CXBIN" ]; then
  cd "$TMP"
  curl -sSL -o accept.html --data-urlencode "file=$RPM" --data-urlencode "choice=Accept" \
       https://www.rbvi.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py
  next=$(grep -oE 'url=[^"]+' accept.html | head -1 | sed 's/^url=//; s/&amp;/\&/g')
  [ -n "$next" ] || { echo "no redirect URL in licence response:"; head -c 600 accept.html; exit 1; }
  curl -sSL -o chimerax.rpm "https://www.rbvi.ucsf.edu$next"
  file chimerax.rpm | grep -q RPM || { echo "download is not an rpm:"; head -c 400 chimerax.rpm; exit 1; }
  rm -rf "$CX"; mkdir -p "$CX"; cd "$CX"; rpm2cpio "$TMP/chimerax.rpm" | cpio -idm --quiet
  rm -f "$TMP/chimerax.rpm" "$TMP/accept.html"
fi

# 2. Stream 9 image
[ -f "$SIF" ] || singularity pull "$SIF" docker://quay.io/centos/centos:stream9

# 3. runtime libraries -> sysroot
if [ ! -f "$SYS/usr/lib64/libOSMesa.so.8" ]; then
  PKGS="mesa-libGL mesa-libGLU mesa-libglapi mesa-compat-libOSMesa mesa-libEGL mesa-libgbm mesa-dri-drivers
        libglvnd libglvnd-glx libglvnd-opengl libglvnd-egl libglvnd-gles llvm-libs
        libX11 libX11-xcb libxcb libXau libXdmcp libXext libXfixes libXxf86vm libXrender libXi libXrandr
        libXcursor libXcomposite libXdamage libXtst libXinerama libXft libXmu libXt libXpm libXv libXScrnSaver libXaw
        libSM libICE libuuid libdrm libxshmfence libxkbcommon libxkbcommon-x11 libxkbfile
        xcb-util xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util-wm xcb-util-cursor
        fontconfig freetype libpng harfbuzz graphite2 libbrotli expat bzip2-libs pcre2 libedit libffi libxml2
        libzstd ncurses-libs zlib xz-libs dbus-libs alsa-lib cups-libs nss nspr nss-util nss-softokn nss-softokn-freebl
        libwayland-client libwayland-server libwayland-egl libwayland-cursor pciutils-libs"
  singularity exec -B "$TMP" "$SIF" dnf -q --setopt=cachedir=$TMP/dnfcache --setopt=logdir=$TMP/dnflog \
      --setopt=persistdir=$TMP/dnfpersist repoquery --latest-limit 1 --location $PKGS \
      | grep -E "^http" | grep -vE "\.i686\.rpm|-devel-" | sort -u > "$TMP/urls.txt"
  wc -l "$TMP/urls.txt"
  mkdir -p "$TMP/rpms"; (cd "$TMP/rpms" && aria2c -q -j8 -i "$TMP/urls.txt" --allow-overwrite=true)
  mkdir -p "$SYS"; cd "$SYS"; for r in "$TMP"/rpms/*.rpm; do rpm2cpio "$r" | cpio -idm --quiet 2>/dev/null || true; done
  rm -rf "$TMP/rpms"
fi
du -sh "$CX" "$SYS" "$SIF"

# 4. smoke test
cat > "$TMP/test.cxc" <<EOF
open /n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025/em/models/9MIP.cif
color bychain
save $TMP/test.png width 600 height 450 supersample 2
exit
EOF
singularity exec -B /n,/scratch,/tmp --env LD_LIBRARY_PATH=$SYS/usr/lib64:$SYS/usr/lib64/dri --env FONTCONFIG_PATH=$SYS/etc/fonts \
    "$SIF" "$CXBIN" --nogui --offscreen --silent "$TMP/test.cxc" 2>&1 | tail -5
[ -s "$TMP/test.png" ] && echo "OFFSCREEN RENDER OK ($(stat -c %s "$TMP/test.png") bytes)" || { echo "render failed"; exit 1; }
find "$CX" -name "*.rpm" -delete 2>/dev/null || true
echo "holylabs:"; quota /n/holylabs/hekstra_lab | tail -1

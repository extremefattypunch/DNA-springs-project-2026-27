#!/bin/bash -l
# Plan section 5A - ChimeraX 1.12 on a Rocky 8 cluster.
#
# UCSF ships ChimeraX 1.12 only for RHEL 9/10 and Ubuntu (glibc >= 2.34); Cannon nodes are
# Rocky 8 (glibc 2.28).  So: download the el9 rpm (POST with the licence 'Accept' choice, as
# the download page's form does; free for academic use), unpack it with rpm2cpio into
# holylabs software/chimerax, and run it inside a Rocky 9 Singularity container that carries
# the runtime libraries.  Offscreen rendering (--nogui --offscreen) uses ChimeraX's bundled
# OSMesa, so no X server is needed.  ~470 MB rpm + ~700 MB container image on holylabs.
#SBATCH -J pp2a-install-chimerax
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 0-01:30
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
CX=$SW/chimerax; SIF=$SW/containers/rocky9-chimerax.sif
TMP=/scratch/${USER}_${SLURM_JOB_ID}; mkdir -p "$TMP" "$SW/containers"; trap 'rm -rf "$TMP"' EXIT
export SINGULARITY_CACHEDIR=$TMP/sifcache SINGULARITY_TMPDIR=$TMP APPTAINER_CACHEDIR=$TMP/sifcache APPTAINER_TMPDIR=$TMP
VER=1.12; RPM="$VER/centos-9/ucsf-chimerax-$VER-1.el9.x86_64.rpm"

if [ ! -x "$CX/opt/UCSF/ChimeraX/bin/ChimeraX" ]; then
  cd "$TMP"
  # Two-step handshake, exactly what the browser does: POST 'Accept' on the licence form returns a
  # page whose <meta refresh> carries a one-time ident token; GET that URL to receive the rpm.
  curl -sSL -o accept.html --data-urlencode "file=$RPM" --data-urlencode "choice=Accept" \
       https://www.rbvi.ucsf.edu/chimerax/cgi-bin/secure/chimerax-get.py
  next=$(grep -oE 'url=[^"]+' accept.html | head -1 | sed 's/^url=//; s/&amp;/\&/g')
  [ -n "$next" ] || { echo "no redirect URL in licence response:"; head -c 600 accept.html; exit 1; }
  curl -sSL -o chimerax.rpm "https://www.rbvi.ucsf.edu$next"
  file chimerax.rpm | grep -q RPM || { echo "download is not an rpm:"; head -c 400 chimerax.rpm; exit 1; }
  rm -f accept.html
  ls -la chimerax.rpm
  rm -rf "$CX"; mkdir -p "$CX"; cd "$CX"; rpm2cpio "$TMP/chimerax.rpm" | cpio -idm --quiet
  rm -f "$TMP/chimerax.rpm"
fi
ls "$CX/opt/UCSF/ChimeraX/bin" | head

# Rocky 9 runtime with the X/GL client libraries ChimeraX links against.
if [ ! -f "$SIF" ]; then
  cat > "$TMP/rocky9.def" <<'EOF'
Bootstrap: docker
From: rockylinux:9
%post
    dnf -y install epel-release || true
    dnf -y install libX11 libXext libXrender libXi libXrandr libXcursor libXcomposite libXdamage libXtst \
                   libxkbcommon libxkbcommon-x11 libglvnd-glx libglvnd-egl libglvnd-opengl mesa-libGL mesa-libEGL \
                   mesa-libgbm mesa-dri-drivers fontconfig freetype dbus-libs alsa-lib nss nspr cups-libs \
                   libxcb xcb-util xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util-wm \
                   xcb-util-cursor libSM libICE libXfixes libXxf86vm libXinerama pulseaudio-libs \
                   xdg-utils which findutils tar gzip && dnf clean all
EOF
  singularity build --fakeroot "$SIF" "$TMP/rocky9.def" 2>&1 | tail -15 || {
    echo "fakeroot build failed; falling back to the plain rockylinux:9 image"
    singularity pull "$SIF" docker://rockylinux:9
  }
fi
ls -la "$SIF"

# Smoke test: headless render of a deposited model to PNG.
cat > "$TMP/test.cxc" <<EOF
open /n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025/em/models/9MIP.cif
color bychain
save $TMP/test.png width 600 height 450 supersample 2
exit
EOF
set +e
singularity exec -B /n,/scratch "$SIF" "$CX/opt/UCSF/ChimeraX/bin/ChimeraX" --nogui --offscreen --silent "$TMP/test.cxc" 2>&1 | tail -20
rc=$?; set -e
if [ -s "$TMP/test.png" ]; then echo "OFFSCREEN RENDER OK: $(stat -c %s "$TMP/test.png") bytes"; cp "$TMP/test.png" /n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025/em/chimerax/smoke_test.png
else echo "offscreen render failed (rc=$rc); missing libraries:"; singularity exec -B /n "$SIF" ldd "$CX/opt/UCSF/ChimeraX/bin/ChimeraX" 2>/dev/null | grep "not found" | head -20 || true; fi
find "$CX" \( -name "*.rpm" -o -name "*.tar.gz" \) -delete 2>/dev/null || true
du -sh "$CX" "$SIF"; echo "holylabs:"; quota /n/holylabs/hekstra_lab | tail -1

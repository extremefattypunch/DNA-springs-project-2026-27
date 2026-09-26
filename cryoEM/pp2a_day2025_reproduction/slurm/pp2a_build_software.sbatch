#!/bin/bash -l
# Plan section 5A - build RELION 3.1 and GROMACS 2023.3 (CPU-only) into holylabs software/.
#
# Why these choices:
#  * test partition: CPU-only build, free (TRESBillingWeights all 0), 12 h cap.
#  * 8 cores / 16 GB: make -j8 on ~1M lines of C++ finishes RELION in ~20 min and
#    GROMACS in ~15 min; memory peak is the compiler (~2 GB/job).
#  * Build in node-local /scratch, install to holylabs, then delete the sources and
#    the FFTW tarballs both builds download (hekstra-cluster rule 3: no archives left).
#  * No CUDA: only relion_postprocess / relion_image_handler / relion_mask_create and
#    gmx rmsf / gmx trjconv are needed - all CPU code paths.
#SBATCH -J pp2a-build-sw
#SBATCH -p test
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 0-04:00
#SBATCH -o /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.out
#SBATCH -e /n/netscratch/hekstra_lab/Lab/ian_poon/logs/%x-%j.err
set -euo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.2.0-fasrc01 cmake/3.28.3-fasrc01

SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
TMP=/scratch/${USER}_${SLURM_JOB_ID}
mkdir -p "$TMP" "$SW"; trap 'rm -rf "$TMP"' EXIT
echo "holylabs before:"; quota /n/holylabs/hekstra_lab | tail -1

# ---------------- RELION 3.1 (paper's local-resolution tool) ----------------
if [ ! -x "$SW/relion-3.1/bin/relion_postprocess" ]; then
  cd "$TMP"; git clone -q --depth 1 -b ver3.1 https://github.com/3dem/relion.git relion31
  mkdir -p relion31/build && cd relion31/build
  set +e
  cmake -DCMAKE_INSTALL_PREFIX="$SW/relion-3.1" -DGUI=OFF -DCUDA=OFF -DALTCPU=OFF \
        -DFORCE_OWN_FFTW=ON -DCMAKE_BUILD_TYPE=Release .. && make -j8 install
  rc=$?; set -e
  if [ $rc -ne 0 ]; then
    echo "RELION 3.1 failed with gcc 12 - retrying with gcc 10.2"
    module swap gcc/12.2.0-fasrc01 gcc/10.2.0-fasrc01
    rm -rf "$TMP/relion31/build"/*; cd "$TMP/relion31/build"
    cmake -DCMAKE_INSTALL_PREFIX="$SW/relion-3.1" -DGUI=OFF -DCUDA=OFF -DALTCPU=OFF \
          -DFORCE_OWN_FFTW=ON -DCMAKE_BUILD_TYPE=Release .. && make -j8 install
    module swap gcc/10.2.0-fasrc01 gcc/12.2.0-fasrc01
  fi
  rm -rf "$TMP/relion31"
fi
"$SW/relion-3.1/bin/relion_postprocess" --version 2>&1 | head -2 || true

# ---------------- GROMACS 2023.3 (paper's RMSF tool) ----------------
if [ ! -x "$SW/gromacs-2023.3/bin/gmx" ]; then
  cd "$TMP"; curl -sSL -o gromacs-2023.3.tar.gz https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
  tar -xzf gromacs-2023.3.tar.gz && rm -f gromacs-2023.3.tar.gz
  mkdir -p gromacs-2023.3/build && cd gromacs-2023.3/build
  # AVX2_256 runs on every Intel/AMD node in test/shared/gpu_test (portable, still fast for analysis tools)
  cmake -DCMAKE_INSTALL_PREFIX="$SW/gromacs-2023.3" -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF \
        -DGMX_MPI=OFF -DGMX_SIMD=AVX2_256 -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_INSTALL_LEGACY_API=OFF \
        -DCMAKE_BUILD_TYPE=Release ..
  make -j8 install
  rm -rf "$TMP/gromacs-2023.3"
fi
source "$SW/gromacs-2023.3/bin/GMXRC"; gmx --version 2>&1 | grep -E "GROMACS version|SIMD" || true

# Hygiene: no tarballs anywhere under software/
find "$SW/relion-3.1" "$SW/gromacs-2023.3" \( -name "*.tar.gz" -o -name "*.tgz" -o -name "*.zip" \) -delete 2>/dev/null || true
du -sh "$SW/relion-3.1" "$SW/gromacs-2023.3"
echo "holylabs after:"; quota /n/holylabs/hekstra_lab | tail -1

# Environment for the PP2A (Day et al. 2025) reproduction.  `source` this file.
#
# Deliberately does NOT source "MD simulations/sfGFP/env/activate.sh" — that file
# exports MAMBA paths on /n/lab_storage, which is 100 % full.  Everything here is
# either in home (code), holylabs (software) or netscratch (data).

export PP2A_REPO=/n/home05/ianpoon/DNA-springs-project-2026-27/cryoEM/pp2a_day2025_reproduction
export PP2A_DATA=/n/netscratch/hekstra_lab/Lab/ian_poon/pp2a_day2025
export PP2A_SW=/n/holylabs/hekstra_lab/Everyone/ianpoon/software
export PP2A_LOGS=/n/netscratch/hekstra_lab/Lab/ian_poon/logs

# Python with OpenMM 8.6 / AmberTools (tleap) / MDAnalysis / mdtraj / mrcfile / PyMOL / ffmpeg.
export DNASPRING_ENV=/n/holylabs/hekstra_lab/Everyone/ianpoon/envs/dnaspring
export DNASPRING_PY=$DNASPRING_ENV/bin/python
export MPLBACKEND=Agg

# CryoSPARC's worker python + pip --target dir holding cryosparc-tools and gemmi.
export CS_WORKER_PY=/n/holylabs/hekstra_lab/Everyone/ianpoon/cryosparc/cryosparc_worker/.pixi/envs/worker/bin/python
export CS_PYLIBS=$PP2A_SW/cs_pylibs

# Phenix 2.2 (labmate's install on lab_storage; read-only use).  Coot is launched through Phenix.
PHENIX_ENV=/n/lab_storage/hekstra_lab/people/aldama/software/phenix2/phenix_env.sh
[ -r "$PHENIX_ENV" ] && source "$PHENIX_ENV"

# Tools built by slurm/pp2a_build_software.sbatch (each block is a no-op until installed).
if [ -d $PP2A_SW/relion-3.1/bin ]; then           # built with gcc 12 + OpenMPI: needs their runtime libs
  command -v module >/dev/null 2>&1 || source /etc/profile.d/lmod.sh
  module load gcc/12.2.0-fasrc01 openmpi/4.1.5-fasrc03 2>/dev/null
  export PATH=$PP2A_SW/relion-3.1/bin:$PATH
fi
[ -d $PP2A_SW/relion-4.0/bin ]      && export PATH=$PP2A_SW/relion-4.0/bin:$PATH
if [ -r $PP2A_SW/gromacs-2023.3/bin/GMXRC ]; then   # GMXRC is not nounset-clean; guard it
  _nounset=$(set -o | awk '$1=="nounset"{print $2}'); set +u
  source $PP2A_SW/gromacs-2023.3/bin/GMXRC
  [ "$_nounset" = on ] && set -u; unset _nounset
fi
# ChimeraX 1.12 (CentOS Stream 9 build) runs inside a Stream 9 container on these Rocky 8 nodes, with the
# X11/GL/OSMesa runtime libraries unpacked from Stream 9 rpms into chimerax/sysroot (no root needed).
# $CHIMERAX is the launcher: e.g.  $CHIMERAX --nogui --offscreen --silent render.cxc
if [ -x $PP2A_SW/chimerax/usr/libexec/UCSF-ChimeraX/bin/ChimeraX ] && [ -f $PP2A_SW/containers/centos-stream9.sif ]; then
  _CXSYS=$PP2A_SW/chimerax/sysroot
  export CHIMERAX="singularity exec -B /n,/scratch,/tmp --env LD_LIBRARY_PATH=$_CXSYS/usr/lib64:$_CXSYS/usr/lib64/dri --env FONTCONFIG_PATH=$_CXSYS/etc/fonts $PP2A_SW/containers/centos-stream9.sif $PP2A_SW/chimerax/usr/libexec/UCSF-ChimeraX/bin/ChimeraX"
fi
[ -d $PP2A_SW/emready ]             && export EMREADY_HOME=$PP2A_SW/emready
[ -d $PP2A_SW/namd-3.0.1 ]          && export PATH=$PP2A_SW/namd-3.0.1:$PATH
[ -d $PP2A_SW/vmd-1.9.4/bin ]       && export PATH=$PP2A_SW/vmd-1.9.4/bin:$PATH
export PATH=$DNASPRING_ENV/bin:$PATH   # ffmpeg, pymol, tleap, cpptraj

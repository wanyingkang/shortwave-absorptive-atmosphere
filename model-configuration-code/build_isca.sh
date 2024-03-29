#!/bin/bash -f
#PBS -A UHAR0008
#PBS -N fixqflux
#PBS -q regular
#PBS -l select=1:ncpus=16:mpiprocs=16:ompthreads=1
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -S /bin/csh -V
#PBS -M kangwanying1992@gmail.com
#PBS -m abe
#############
# environment:
# export TMPDIR=/glade/scratch/$USER/temp
# mkdir -p $TMPDIR
export GFDL_BASE=/glade/scratch/$USER/Isca1
source ${GFDL_BASE}/src/extra/env/${GFDL_ENV}

# passing info
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

export casename=${1:-"test"}
export ncores=${2:-"16"}
export todo=${3:-"compile"}
export mothercase=${4:-""}
echo casename=$casename, ncores=$ncores, todo=$todo, mothercase=$mothercase

# running model
python build_isca.py $casename $ncores $todo $mothercase


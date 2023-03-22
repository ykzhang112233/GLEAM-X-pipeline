#! /bin/bash

usage()
{
echo "drift_postmosaic.sh [-p project] [-d dep] [-t] [-m mosaicdir] [-l lowres_freq] [-r highres_freq] [-c comb_freq] -n mosaic_name
  -p project  : project, (must be specified, no default)
  -d dep      : job number for dependency (afterok)
  -t          : test. Don't submit job, just make the batch file
                and then return the submission command  
  -m mosaicdir: Directory name where mosaics stored (default=project/mosaic) 
  -l lowres   : Frequency range of the lowres image. (default=170-200MHz)
  -r highres  : Frequency range of the highres image (default=200-231MHz)
  -c comb     : Frequency range of combined image (default=170-231MHz)
  -n mosaicnm : Name of the mosaic, typically the drift name e.g. XG_D-27_20201015 (no default, must be specified)" 1>&2;
exit 1;
}

pipeuser=$(whoami)

# initial variables 
dep=
queue="-p ${GXSTANDARDQ}"
tst=
mosaicdir=
lowres_freq=
highres_freq=
comb_freq=

# parse args and set options
while getopts ':td:p:b:m:l:r:c:n:' OPTION
do
    case "$OPTION" in
    d)
        dep=${OPTARG} ;;
    p)
        project=${OPTARG} ;;
    m) 
        mosaicdir=${OPTARG} ;;
    t)
        tst=1 ;;
    l)
        lowres_freq=${OPTARG} ;;
    r)
        highres_freq=${OPTARG} ;;
    c)
        comb_freq=${OPTARG} ;;
    n)
        mosaicnm=${OPTARG} ;;
    ? | : | h)
            usage ;;
  esac
done

if [[ -z ${project} ]] || [[ -z ${mosaicnm} ]]
then
    usage
fi

if [[ ! -z ${dep} ]]
then
    depend="--dependency=afterok:${dep}"
fi

if [[ ! -z ${GXACCOUNT} ]]
then
    account="--account=${GXACCOUNT}"
fi

base="${GXSCRATCH}/${project}"
listbase=$(basename "${mosaicnm}")
listbase=${listbase%%.*}
script="${GXSCRIPT}/postmosaic_${listbase}.sh"


cat "${GXBASE}/templates/postmosaic.tmpl" | sed -e "s:BASEDIR:${base}:g" \
                                                -e "s:PIPEUSER:${pipeuser}:g" \
                                                -e "s:MOSAICNM:${mosaicnm}:g" \
                                                -e "s:MOSAICDIR:${mosaicdir}:g" \
                                                -e "s:LOWRES_FREQ:${lowres_freq}:g" \
                                                -e "s:HIGHRES_FREQ:${highres_freq}:g" \
                                                -e "s:COMB_FREQ:${comb_freq}:g" \


output="${GXLOG}/postmosaic_${listbase}.o%A_%a"
error="${GXLOG}/postmosaic_${listbase}.e%A_%a"

chmod 755 "${script}"

# sbatch submissions need to start with a shebang
echo '#!/bin/bash' > "${script}.sbatch"
echo "srun --cpus-per-task=${GXNCPUS} --ntasks=1 --ntasks-per-node=1 singularity run ${GXCONTAINER} ${script}" >> "${script}.sbatch"

# Automatically runs a job array for each sub-band
sub="sbatch  --begin=now+5minutes --export=ALL  --time=10:00:00 --mem=${GXABSMEMORY}G -M ${GXCOMPUTER} --output=${output} --error=${error}"
sub="${sub} ${GXNCPULINE} ${account} ${GXTASKLINE} ${depend} ${queue} ${script}.sbatch"
if [[ ! -z ${tst} ]]
then
    echo "script is ${script}"
    echo "submit via:"
    echo "${sub}"
    exit 0
fi

# submit job
jobid=($(${sub}))
jobid=${jobid[3]}

# rename the err/output files as we now know the jobid
error=${error//%A/"${jobid}"}
output=${output//%A/"${jobid}"}

echo "${output}"
echo "${error}"

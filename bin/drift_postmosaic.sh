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
                                                -e "s:COMB_FREQ:${comb_freq}:g" > ${script}


output="${GXLOG}/postmosaic_${listbase}.o%A"
error="${GXLOG}/postmosaic_${listbase}.e%A"

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

# Ok this is a bold choice, but we're going to try submitting an array job AFTER the postmosaic job (only to run once the postmosaic has been given the all clear) and then this next part will submit an array of jobs that is the priorized fitting based on the output of the postmosaic. Could I do them as separate scripts? Probably but I like having all the "post" stuff together so lets give it a whirl 

# OK SINCE YOU KNOW HOW MANY THERE ARE GOING TO BE YOU CAN JUST HAVE ARRAY= AND HARD CODE IT BECAUSE THE BUSINESS ACTUALLY COMES IN THE TEMPLATE AND I ALREADY SET THAT UP
script="${GXSCRIPT}/priorized_fitting_${listbase}.sh"

cat "${GXBASE}/templates/priorized_fitting.tmpl" | sed -e "s:BASEDIR:${base}:g" \
                                                -e "s:PIPEUSER:${pipeuser}:g" \
                                                -e "s:MOSAICNM:${mosaicnm}:g" \
                                                -e "s:MOSAICDIR:${mosaicdir}:g" \
                                                -e "s:COMB_FREQ:${comb_freq}:g" > ${script}

output="${GXLOG}/priorized_fitting_${listbase}.o%A_%a"
error="${GXLOG}/priorized_fitting_${listbase}.e%A_%a"

chmod 755 "${script}"

# sbatch submissions need to start with a shebang
echo '#!/bin/bash' > "${script}.sbatch"
echo "srun --cpus-per-task=${GXNCPUS} --ntasks=1 --ntasks-per-node=1 singularity run ${GXCONTAINER} ${script}" >> "${script}.sbatch"

# Automatically runs a job array for each sub-band
sub="sbatch  --begin=now+5minutes --dependency=afterok:${jobid} --array=25  --export=ALL  --time=10:00:00 --mem=${GXABSMEMORY}G -M ${GXCOMPUTER} --output=${output} --error=${error}"
sub="${sub} ${GXNCPULINE} ${account} ${GXTASKLINE} ${queue} ${script}.sbatch"
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



script="${GXSCRIPT}/join_cats_${listbase}.sh"

cat "${GXBASE}/templates/join_cats.tmpl" | sed -e "s:BASEDIR:${base}:g" \
                                                -e "s:PIPEUSER:${pipeuser}:g" \
                                                -e "s:MOSAICNM:${mosaicnm}:g" \
                                                -e "s:MOSAICDIR:${mosaicdir}:g" \
                                                -e "s:COMB_FREQ:${comb_freq}:g" > ${script}

output="${GXLOG}/join_cats_${listbase}.o%A"
error="${GXLOG}/join_cats_${listbase}.e%A"

chmod 755 "${script}"

# sbatch submissions need to start with a shebang
echo '#!/bin/bash' > "${script}.sbatch"
echo "srun --cpus-per-task=${GXNCPUS} --ntasks=1 --ntasks-per-node=1 singularity run ${GXCONTAINER} ${script}" >> "${script}.sbatch"

# Automatically runs a job array for each sub-band
sub="sbatch  --begin=now+5minutes --dependency=afterok:${jobid} --export=ALL  --time=02:00:00 --mem=${GXABSMEMORY}G -M ${GXCOMPUTER} --output=${output} --error=${error}"
sub="${sub} ${GXNCPULINE} ${account} ${GXTASKLINE} ${queue} ${script}.sbatch"
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



# TODO: Double check, the nextflow stuff has a move to rename things but I think I've already got things named how i want? check. 

# python join_catalogues.py --epochs "${mosaicnm}_catalogues.csv" --refcat "${mosaicnm}_${comb_freq}_comp_rescaled.fits" --out "${mosaicnm}_joined_comp.vot" --all 

# Note the nextflow has a line that makes a new rescaled .csv file, I already make it in the prep python script so no need 
# python join_cataloges.py --epochs "${mosaicnm}_catalogues_rescaled.csv" --refcat "${mosaicnm}_${comb_freq}_comp_rescaled.fits" --out "${mosaicnm}_joined_rescaled_comp.vot" --all 



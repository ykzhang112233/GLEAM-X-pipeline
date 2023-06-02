#!/bin/sh
#SBATCH --partition=copy
#SBATCH -M=setonix 
#SBATCH --time=48:00:00
#SBATCH --job-name=download_nimbus_XG_D-55_20201020
#SBATCH --output=${GXLOG}/download_nimbus_XG_D-55_20201020.o%A
#SBATCH --error=${GXLOG}/download_nimbus_XG_D-55_20201020.e%A


VMUSER=ubuntu
VMADDR="146.118.68.233"
VMPATH="/mnt/gxarchive/Archived_Obsids"
SAVEPATH="/astro/mwasci/kross/gleamx/XG_D-55_4Night/20201020/"
FILENAME="/astro/mwasci/kross/gleamx/XG_D-55_4Night/20201020/XG_D-55_20201020_redownload.txt"

DRYRUN=''

echo "We will be downloading ionosphere data from"
echo "  - USER: ${VMUSER}"
echo "  - HOST: ${VMADDR}"
echo "  - PATH: ${VMPATH}"

echo """

Some of these commands will take a long time to start. The directories stored on the virtual machine are mounted across a network (a network link that spans Australia!). There might be periods where this process takes a long time to complete. Rest assured the rsync should be doing the right thing.  Just make sure you've set the filename correctly and it should work :) 

"""
set -x


## touch includefile.txt 
readarray -t arr < ${FILENAME}

for element in ${arr[@]}
    do
        rsync -avh $DRYRUN --progress "${VMUSER}@${VMADDR}:${VMPATH}/${element}/*-image-pb_warp_rescaled.fits" "${SAVEPATH}/${element}/"
        rsync -avh $DRYRUN --progress "${VMUSER}@${VMADDR}:${VMPATH}/${element}/*-image-pb_warp_rescaled_weight.fits" "${SAVEPATH}/${element}/"
        rsync -avh $DRYRUN --progress "${VMUSER}@${VMADDR}:${VMPATH}/${element}/${element}_deep-MFS-image-pb_warp_rms.fits" "${SAVEPATH}/${element}/${element}_deep-MFS-image-pb_warp_rms.fits"
        rsync -avh $DRYRUN --progress "${VMUSER}@${VMADDR}:${VMPATH}/${element}/${element}_deep-MFS-image-pb_warp_comp.fits" "${SAVEPATH}/${element}/${element}_deep-MFS-image-pb_warp_comp.fits"
        wget -O ${element}.metafits http://ws.mwatelescope.org/metadata/fits?obs_id=${element}
            mv ${element}.metafits "${SAVEPATH}/${element}/${element}.metafits"
done

lastobsid=$(tail -n 1 ${FILENAME})

if [[ -e "${SAVEPATH}/${lastobsid}/${lastobsid}.metafits" ]]
then
        echo "Finished downloading everything!"
        exit 0
else
        echo "Something is wrong? Dont have the metafits for last obsid..." 
        exit 1
fi
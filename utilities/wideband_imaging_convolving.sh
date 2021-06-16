#!/bin/bash

# A small utility script to brutally coadd two mosaic images together. 
# It is only meant to be used to provide the initial source finding with 
# a deeper image from which the priorised fitting is based on. 

if [[ $# -ne 3 ]]
then 
    echo "USAGE: $0 LOWRES HIGHRES OUTTITLE"
    echo "LOWRES - The low resolution mosaic, i.e. XG_170-200MHz_ddmod.fits. "
    echo "LOWRES - The high resolution mosaic that will be convolved to LOWRES, i.e. XG_200-230MHz_ddmod.fits. "
    echo "OUTTITLE - String used as the basis of the new file created, i.e. XG_170-230MHz "
    exit
fi

SCRIPTSUFFIX=widebandmir

find . -iname "*.${SCRIPTSUFFIX}" -type d -exec rm -rf {} +

lowres_im=$1
highres_im=$2
combined_im=$3

# Reading both images into a miriad format
fits in="${lowres_im}" out="${lowres_im/fits/${SCRIPTSUFFIX}}" op=xyin
fits in="${highres_im}" out="${highres_im/fits/${SCRIPTSUFFIX}}" op=xyin

prthd in="${lowres_im/fits/${SCRIPTSUFFIX}}" 

lowres_fwhm_a=$(prthd in="${lowres_im/fits/${SCRIPTSUFFIX}}" | grep Beam | tr -s ' ' | cut -d ' ' -f3)
lowres_fwhm_b=$(prthd in="${lowres_im/fits/${SCRIPTSUFFIX}}" | grep Beam | tr -s ' ' | cut -d ' ' -f5)
lowres_pos_ang=$(prthd in="${lowres_im/fits/${SCRIPTSUFFIX}}" | grep Position | tr -s ' ' | cut -d ' ' -f3)

echo "Extracted FWHM of low-resolution image: ${lowres_fwhm_a}x${lowres_fwhm_b} and ${lowres_pos_ang}"

# Regriding the lowres to match highres
regrid in="${lowres_im/fits/${SCRIPTSUFFIX}}" \
        out="${lowres_im/fits/regrid.${SCRIPTSUFFIX}}" \
        tin="${highres_im/fits/${SCRIPTSUFFIX}}"

# Convolving the high res to low res 
convol map="${highres_im/fits/${SCRIPTSUFFIX}}" \
      fwhm="${lowres_fwhm_a},${lowres_fwhm_b}" \
      pa="${lowres_pos_ang}" \
      options=final \
      out="${highres_im/fits/convol.${SCRIPTSUFFIX}}" 


final_lowmir_im="${lowres_im/fits/regrid.${SCRIPTSUFFIX}}"
final_highmir_im="${highres_im/fits/convol.${SCRIPTSUFFIX}}" 
outfile="${combined_im}.${SCRIPTSUFFIX}"

# Averaging the two imgaes 
maths exp="'(<${final_lowmir_im}>+<${final_highmir_im}>)/2'" out="${outfile}"

# Exporting the miriad to a regular image 
fits in="${outfile}" out="${outfile/${SCRIPTSUFFIX}/fits}" op=xyout 

find . -iname "*.${SCRIPTSUFFIX}" -type d -exec rm -rf {} +

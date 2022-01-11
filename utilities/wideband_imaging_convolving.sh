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

lowpsf="${lowres_im/ddmod/projpsf_psf}"
highpsf="${highres_im/ddmod/projpsf_psf}"

# extract beamsizes
out=$(extract_lowhigh_psf_beam.py $lowpsf $highpsf -p)
echo "${out}"
low_maj=$(echo "${out}" | grep 'Low' | cut -d ' ' -f4)
low_min=$(echo "${out}" | grep 'Low' | cut -d ' ' -f5)

high_maj=$(echo "${out}" | grep 'High' | cut -d ' ' -f4)
high_min=$(echo "${out}" | grep 'High' | cut -d ' ' -f5)

# Reading both images into a miriad format
fits in="${lowres_im}" out="${lowres_im/fits/${SCRIPTSUFFIX}}" op=xyin
fits in="${highres_im}" out="${highres_im/fits/${SCRIPTSUFFIX}}" op=xyin

# Put the extracted values into the files, even the low frequency ones. Not entirely sure how often
# they are used by miriad, so lets just be sure they are in.
puthd in="${highres_im/fits/${SCRIPTSUFFIX}}/bmaj" value="${high_maj},arcseconds"
puthd in="${highres_im/fits/${SCRIPTSUFFIX}}/bmin" value="${high_min},arcseconds"
puthd in="${lowres_im/fits/${SCRIPTSUFFIX}}/bmaj" value="${low_maj},arcseconds"
puthd in="${lowres_im/fits/${SCRIPTSUFFIX}}/bmin" value="${low_min},arcseconds"

# Get out the values, we just put in, needlessly. See above message about why it is done this way.
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

puthd in="${outfile}"/freq value=200315000

# Exporting the miriad to a regular image
fits in="${outfile}" out="${outfile/${SCRIPTSUFFIX}/fits}" op=xyout

find . -iname "*.${SCRIPTSUFFIX}" -type d -exec rm -rf {} +

if [[ -z $GXNCPUS ]]
then
    export GXNCPUS=8
    echo "Environment variable GXNCPUS not found. Setting to $GXNCPUS"
fi

BANE --cores ${GXNCPUS} \
--compress \
--noclobber \
"${outfile/${SCRIPTSUFFIX}/fits}"

aegean \
--seedclip=10 \
--maxsummits=5 \
--cores 1 \
--progress \
--autoload \
--table="${outfile/.${SCRIPTSUFFIX}/_projpsf.fits}" \
"${outfile/${SCRIPTSUFFIX}/fits}"

psf_select.py --input="${combined_im}_projpsf_comp.fits"
psf_create.py --input="${combined_im}_projpsf_comp_psfcat.fits"

aegean \
--seedclip=4 \
--maxsummits=5 \
--cores 1 \
--autoload \
--progress \
--psf="${outfile/.${SCRIPTSUFFIX}/_projpsf_psf.fits}" \
--table="${combined_im}.fits" \
"${outfile/${SCRIPTSUFFIX}/fits}"


mosaic_global_rescale.py \
"${combined_im}_comp.fits" \
"${outfile/${SCRIPTSUFFIX}/fits}" \
"${GXBASE}/models/GGSM_sparse_unresolved.fits" \
--plot \
--verbose \
--apply

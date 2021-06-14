#!/bin/csh

project_dir=/Users/katross/Documents/PhD/MWA/GLEAM-X/wideband_imaging
cd ${project_dir}
# Remove any previous .mir files from the last run (maybe make less hacky)
rm -r *.mir 

ls

lowres_im=GS_170-200MHz_ddmod_cutout
highres_im=GS_200-231MHz_ddmod_cutout
combined_im=GS_170-231

lowres_fwhm_a=75.01355999978
lowres_fwhm_b=58.7584799998

echo ${lowres_im}.mir
# Reading both images into a miriad format
fits in=${lowres_im}.fits out=${lowres_im}.mir op=xyin
fits in=${highres_im}.fits out=${highres_im}.mir op=xyin

# Regriding the lowres to match highres
regrid in=${lowres_im}.mir out=${lowres_im}_regrid.mir tin=${highres_im}.mir

# Convolving the high res to low res 
convol map=${highres_im}.mir fwhm=${lowres_fwhm_a},${lowres_fwhm_b} options=final out=${highres_im}_convol.mir 


final_lowmir_im=GS_170-200MHz_ddmod_cutout_regrid.mir
final_highmir_im=GS_200-231MHz_ddmod_cutout_convol.mir
outfile=${combined_im}.mir

# Averaging the two imgaes 
maths exp="(<${final_lowmir_im}>+<${final_highmir_im}>)/2" out=${combined_im}.mir

# Exporting the miriad to a regular image 
fits in=${combined_im}.mir out=${combined_im}.fits op=xyout 


# BANE AND AEGEAN?? 
BANE ${combined_im}.fits 
BANE ${lowres_im}.fits 
BANE ${highres_im}.fits 

aegean --autoload --table ${combined_im}.fits --cores 3 ${combined_im}.fits
aegean --autoload --table ${highres_im}.fits --cores 3 ${highres_im}.fits
aegean --autoload --table ${lowres_im}.fits --cores 3 ${lowres_im}.fits

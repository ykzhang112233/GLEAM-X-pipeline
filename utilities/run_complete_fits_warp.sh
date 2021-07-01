#!/bin/bash

# A small utility script to run fits_warp.py on each obsid, producing just the
# fits_warp plot and XM for all sources found with the internal cross-match

base=$(pwd)
POS_REF=$GXBASE/models/NVSS_SUMSS_psfcal.fits

for obsnum in 12*
do
    cd $base/$obsnum
    for subchan in 0000 0001 0002 0003 MFS
    do
        fits_warp.py \
        --incat ${obsnum}_deep-${subchan}-image-pb_comp.fits \
        --refcat $POS_REF \
        --xm ${obsnum}_${subchan}_complete_sources_xm.fits \
        --plot \
        --ra1 ra \
        --dec1 dec \
        --ra2 RAJ2000 \
        --dec2 DEJ2000 \
        --infits ${obsnum}_deep-${subchan}-image-pb.fits
    done
done
#!/usr/bin/env python

""" Script to assess the quality of obsid images and identify high quality images to include in mosaic step. 

TODO: generalise to run for night rather than per channel. Make sure to check percentage of missing per night to log issue rather than just a raw number

TODO: add plotting options 

"""


import sys
import os
from argparse import ArgumentParser
from astropy.io import fits
from matplotlib import pyplot as plt
import astropy.units as u
import numpy as np
import logging 

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to assess the quality of images for obsids and return list of obsids that pass quality assurance to be included in moasics. Note: currently only works on obsids given per channel, not per night. "
    )
    parser.add_argument(
        '--flag_high_rms',
        default=True,
        help='Will only select obsids that have RMS no larger than the median RMS at that freq channel + the STD of the RMS of the night'
    )
    parser.add_argument(
        '--flag_bad_io',
        defualt=True, 
        help='Will run the selection cuts for bad ionosphere on individial obsids based on quality checks from Brandon'
    )
    parser.add_argument(
        '-save_bad_obsids',
        type=str,
        help="If defined, will make a .txt file with the bad obsids"
    )
    parser.add_argument(
        'save_missing_obsids',
        type=str,
        help="If defined, will make a .txt file in directory with all obsids with no *MFS-image-pb_warp_rms.fits file"
    )
    parser.add_argument(
        'obsids',
        type=str,
        help="Path to the .txt file with the new line separated obsids for the given channel"
    )
    parser.add_argument(
        'output',
        type=str,
        default=".",
        help="directory of .txt file for output, assumes a .txt file in current directory named after input obsids text file"
    )
    parser.add_argument(
        '-b',
        '--base_path',
        type=str,
        default=".",
        help="Path to folder containing obsid folders"
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help='Enable extra logging'
    )
    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)


    obsids =  np.loadtxt(args.obsids)


def read_rms(obsids):
    rms = []
    ra = [] #needed for the plotting, adding currently but not necessary
    missing_obsids = [] #badobsids is for the missing ones for reference of what to check up on 
    
    for obs in obsids: 
        rmsfile = f"{args.base_path}/{obs:10.0f}/{obs:10.0f}_deep-MFS-image-pb_warp_rms.fits"
        if os.path.exists(rmsfile):
            # plotobs.append(obs)
            hdu = fits.open(rmsfile)
            rms.append(1.e3*hdu[0].data[int(hdu[0].data.shape[0]/2), int(hdu[0].data.shape[1]/2)])
            ra.append(hdu[0].header["CRVAL1"])
            hdu.close()
        else:
            rms.append(np.nan)
            # plotobs.append(np.nan)
            ra.append(np.nan)
            missing_obsids.append(obs)

    if len(missing_obsids)>5:
        logger.warning(f"Large number of missing obsids: {len(missing_obsids)}/{len(obsids)}")
    if args.save_missing_obsids is not None:
        np.savetxt(args.obsids.replace(".txt", "_missing_obsids.txt"), missing_obsids, fmt="%10.0f")


    cutoff = np.nanmedian(rms)+np.nanstd(rms)
    obslist = obsids[rms < cutoff]
    
    # Just shouting out that many are high RMS
    # TODO: implement from here a QA check that potentially cans the night if all channels have this or it's some ridiculously large number
    frac_flagged = (len(obslist)/len(obsids))*100
    if frac_flagged > 15:
        logger.warning(f"Large number of obsids flagged for high RMS: {frac_flagged}%")


    # np.savetxt(args.obsids.replace(".txt", "_selected_by_RMS.txt"), obslist, fmt="%10.0f")

    return 
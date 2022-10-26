#!/usr/bin/env python

""" Script to assess the quality of obsid images and identify high quality images to include in mosaic step. 

TODO: generalise to run for night rather than per channel. Make sure to check percentage of missing per night to log issue rather than just a raw number

TODO: add plotting options 

"""


from ast import parse
from mmap import mmap
import sys
import os
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from matplotlib import pyplot as plt
import astropy.units as u
import numpy as np
import logging 

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)


def cut_high_rmsobsids(
    obsids
):#, base_dir=args.base_path):

    rms = []
    ra = [] #needed for the plotting, adding currently but not necessary
    missing_obsids = [] #badobsids is for the missing ones for reference of what to check up on 
    
    for obs in obsids: 
        rmsfile = f"{base_dir}/{obs:10.0f}/{obs:10.0f}_deep-MFS-image-pb_warp_rms.fits"
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
    # TODO: implement from here a QA check that potentially cans the night if all channels have this or it's some ridiculously large
    frac_flagged = (len(obslist)/len(obsids))*100
    if frac_flagged > 15:
        logger.warning(f"Large number of obsids flagged for high RMS: {frac_flagged}%")

    return obslist 

def crossmatch_cats(
    input_cat,
    ref_cat,
    sep,
):
    ref_cat_skycoords = SkyCoord(ref_cat.RAJ2000,ref_cat.DEJ2000, frame='fk5',unit=u.deg)
    input_cat_skycoords = SkyCoord(input_cat.ra, input_cat.dec, frame='fk5', unit=u.deg)
    
    # Crossmatching the two catalogues with sep 1arcmin default 
    # i.e. idx1 is list of indexes for inputcat that have corresponding crossmatch in refcat (but I don't actually care about refcat, just that it's there)
    idx1,idx2,sep2d,dist3d=search_around_sky(input_cat_skycoords,ref_cat_skycoords,sep*u.arcmin)

    output_cat = input_cat[idx1]

    return output_cat



def cut_cat_bright(
    base_dir, 
    obsid,
    refcat,
    # cut_type=cut_type, 
    # cut_level_int_peak=cut_level_int_peak,
    # cut_level_int_rms=cut_level_int_rms,
    # cut_level_snr=cut_level_snr,

):

    try:
        temp = fits.open(f"{base_dir}/{obsid}/{obsid}_deep-MFS-image-pb_warp_rescaled_comp.fits",mmap=True)
        temp_cat = temp[1].data
    except:
        logger.debug(f"No catalogue for io checks: {obsid}")
        return 

    int_over_peak = temp_cat["int_flux"]/temp_cat["peak_flux"]
    err_intoverrms = temp_cat["err_int_flux"]/temp_cat["local_rms"]
    snr = temp_cat["int_flux"]/temp_cat["local_rms"]
    blur = np.log10(temp_cat['int_flux']/temp_cat['peak_flux'])
    std_int_over_peak = np.nanstd(temp_cat['int_flux']/temp_cat['peak_flux'])

    temp_cat["int_over_peak"] = int_over_peak
    temp_cat["err_int_over_rms"] = err_intoverrms
    temp_cat["snr"] = snr
    temp_cat['blur'] = blur
    temp_cat['std_int_over_peak'] = std_int_over_peak

    if cut_type == "suggested":
        mask = np.where((int_over_peak<=cut_level_int_peak)&(err_intoverrms<=cut_level_int_rms))
        cat = temp_cat[mask]

    elif cut_type == "snr":
        mask = np.where((int_over_peak<=cut_level_int_peak)&(snr>=cut_level_snr))
        cat = temp_cat[mask]

    elif cut_type == "both":
        mask = np.where((int_over_peak<=cut_level_int_peak)&(err_intoverrms<=cut_level_int_rms)&(snr>=cut_level_snr))
        cat = temp_cat[mask]
    else: 
        logger.warning(f"No cut defined?!?! Carrying on with NO any snr or intoverpeak cuts")
        cat = temp_cat


    cat_xm = crossmatch_cats(cat,refcat)

    return cat_xm


def calc_intoverpeak_stats(
    base_dir,
    obslist,
):

    mean_int_over_peak = []
    std_int_over_peak = []
    for obs in obslist: 
        obsid_cat = cut_cat_bright(base_dir,obs)

        # TODO: add plotting option here to plot the obsid intoverpeak per source or something

        mean_int_over_peak.append(np.nanmedian(obsid_cat['int_over_peak']))
        std_int_over_peak.append(np.nanstd(obsid_cat['int_over_peak']))


    return mean_int_over_peak, std_int_over_peak




if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to assess the quality of images for obsids and return list of obsids that pass quality assurance to be included in moasics. Note: currently only works on obsids given per channel, not per night. "
    )
    parser.add_argument(
        '-b',
        '--base_path',
        type=str,
        default=".",
        help="Path to folder containing obsid folders"
    )
    parser.add_argument(
        'obsids',
        type=str,
        help="Path to the .txt file with the new line separated obsids for the given channel"
    )
    parser.add_argument(
        "--catalogue",
        type=str,
        dest="cat",
        default="GGSM_sparse_unresolved.fits",
        help="Filename of catalogue to use for comparison, assumed GGSM_sparse_unresolved from GLEAM-X-pipeline",
    )    
    parser.add_argument(
        '--sep',
        type=float,
        default=1,
        help="Separation from reference catalogue to that of obsid. (defualt=1arcmin)"
    )



    # parser.add_argument(
    #     "--racol",
    #     type=str,
    #     default="RAJ2000",
    #     help="The name of the RA column (in decimal degrees) in your catalogue (default = RAJ2000)",
    # )

    # parser.add_argument(
    #     "--decol",
    #     type=str,
    #     default="DEJ2000",
    #     help="The name of your Dec column (in decimal degrees) in your catalogue (default = DEJ2000)",
    # )
    parser.add_argument(
        'output',
        type=str,
        default=".",
        help="directory of .txt file for output, assumes a .txt file in current directory named after input obsids text file"
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
        '--cut_type',
        help='Can pick either "snr", "suggested" or "both". Suggested is based on two factors, int over peak ratio and int_err vs rms ratio. SNR is based on the SNR as well as the int over peak comparison ratio.',
        default='both',
        type=str
    )
    parser.add_argument(
        '--cut_level_int_over_peak',
        help='The ratio of int flux to peak flux below which the sources will be retained. (Default=2).',
        default=2,
        type=float
    )
    parser.add_argument(
        '--cut_level_err_int_rms',
        help='The cut level for the err_int and local_rms, below which sources will be kept. Default=2.'
        ,default=2,
        type=float
    )
    parser.add_argument(
        '--cut_level_snr',
        help='The SNR ratio at which the sources will be cut. Default=5',
        default=5,
        type=float
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
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help='Enable extra logging'
    )


    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Defining all the global variables based on inputs 
    obsids =  np.loadtxt(args.obsids)
    base_dir = args.base_path
    cut_type=args.cut_type
    cut_level_int_rms=args.cut_level_err_int_rms
    cut_level_int_peak=args.cut_level_int_over_peak
    cut_level_snr=args.cut_level_SNR


    # Double check: I don't think you need this anymore 
    ggsm = args.cat 
    ggsm_cat = ggsm[1].data
    coords = SkyCoord(ggsm_cat[args.racol], ggsm_cat[args.decol], unit=(u.deg, u.deg))


    # Actually implementing the cuts using functions above 
    # First up: checking to extract only obsids with nice RMS 
    if args.flag_high_rms is True: 
        obslist = cut_high_rmsobsids(obsids)
    else:
        logger.debug("Not running high RMS flagger")
        obslist = obsids

    # Establishing the stats needed to assess quality of night here, not actually doing any of the cuts here just assessing
    if args.flag_bad_io is True: 

        mean_int_over_peak = []
        std_int_over_peak = []
        for obs in obslist: 
            obsid_cat = cut_cat_bright(base_dir,obs)

            # TODO: add plotting option here to plot the obsid intoverpeak per source or something
            mean_int_over_peak.append(np.nanmedian(obsid_cat['int_over_peak']))
            std_int_over_peak.append(np.nanstd(obsid_cat['int_over_peak']))

    else: 
        logger.debug(f"Not running the ionosphere analysis")

    # TODO: add plotting option to plot the int/peak+/- std as function of obsid over the night. 
    
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
import matplotlib.ticker as ticker
from matplotlib import rcParams
import matplotlib.pyplot as plt



logger = logging.getLogger(__name__)
logging.basicConfig(format="%(module)s:%(levelname)s:%(lineno)d %(message)s")
logger.setLevel(logging.INFO)


rcParams['font.family'] = 'serif'
s, dt, axlab = 8, 1.1, 12.
plt.rcParams["xtick.major.size"] = s
plt.rcParams["xtick.minor.size"] = s
plt.rcParams["ytick.major.size"] = s
plt.rcParams["ytick.minor.size"] = s
plt.rcParams["xtick.major.width"] = dt
plt.rcParams["xtick.minor.width"] = dt
plt.rcParams["ytick.major.width"] = dt
plt.rcParams["ytick.minor.width"] = dt
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.major.pad"] = 5.
plt.rcParams["figure.figsize"] = [10., 4.5]


def cut_high_rmsobsids(
    obsids
):#, base_dir=args.base_path):

    rms = []
    ra = [] #needed for the plotting, adding currently but not necessary
    missing_obsids = [] #badobsids is for the missing ones for reference of what to check up on 
    
    for obs in obsids: 
        rmsfile = f"{base_dir}/{obs:10.0f}/{obs:10.0f}_deep-MFS-image-pb_warp_rms.fits"
        if os.path.exists(rmsfile):
            # TODO: plotobs.append(obs)
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
    if args.save_missing_obsids is True:
        logger.debug(f"Saving textfile with list of missing obsids: {args.save_missing_obsids}")
        np.savetxt(args.obsids.replace(".txt", "_missing_obsids.txt"), missing_obsids, fmt="%10.0f")


    cutoff = np.nanmedian(rms)+np.nanstd(rms)
    obslist = obsids[rms < cutoff]
    
    # Just shouting out that many are high RMS
    # TODO: implement from here a QA check that potentially cans the night if all channels have this or it's some ridiculously large
    frac_flagged = int((1-(len(obslist)/len(obsids)))*100)
    if frac_flagged > 15:
        logger.warning(f"Large number of obsids flagged for high RMS: {frac_flagged}%")

    # TODO: add save bad obsids list
    # if args.save_bad_obsids is not None:
    #     np.savetxt(args.obsids.replace(".txt", "_missing_obsids.txt"), missing_obsids, fmt="%10.0f")


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
    # refcat,
):

    try:
        temp = fits.open(f"{base_dir}/{obs:10.0f}/{obs:10.0f}_deep-MFS-image-pb_warp_rescaled_comp.fits",mmap=True)
        temp_cat = temp[1].data
        temp.close()

    except:
        logger.debug(f"No catalogue for io checks: {obs:10.0f}/{obs:10.0f}")
        return 

    int_over_peak = temp_cat["int_flux"]/temp_cat["peak_flux"]
    err_intoverrms = temp_cat["err_int_flux"]/temp_cat["local_rms"]
    snr = temp_cat["int_flux"]/temp_cat["local_rms"]
    blur = np.log10(temp_cat['int_flux']/temp_cat['peak_flux'])
    std_int_over_peak = np.nanstd(temp_cat['int_flux']/temp_cat['peak_flux'])

    if args.cut_type == "suggested":
        mask = np.where((int_over_peak<=args.cut_level_int_peak)&(err_intoverrms<=args.cut_level_int_rms))
        cat = temp_cat[mask]

    elif args.cut_type == "snr":
        mask = np.where((int_over_peak<=args.cut_level_int_peak)&(snr>=args.cut_level_snr))
        cat = temp_cat[mask]

    elif args.cut_type == "both":
        mask = np.where((int_over_peak<=args.cut_level_int_peak)&(err_intoverrms<=args.cut_level_err_int_rms)&(snr>=args.cut_level_snr))
        cat = temp_cat[mask]
    else: 
        logger.warning(f"No cut defined?!?! Carrying on with NO any snr or intoverpeak cuts")
        cat = temp_cat
    
    # cat_xm = crossmatch_cats(cat,refcat)

    return cat


def calc_intoverpeak_stats(
    base_dir,
    obslist,
):

    mean_int_over_peak = []
    std_int_over_peak = []
    for obs in obslist: 
        obsid_cat = cut_cat_bright(base_dir,obs)

        int_over_peak = obsid_cat["int_flux"]/obsid_cat["peak_flux"]
        err_intoverrms = obsid_cat["err_int_flux"]/obsid_cat["local_rms"]
        snr = obsid_cat["int_flux"]/obsid_cat["local_rms"]
        blur = np.log10(obsid_cat['int_flux']/obsid_cat['peak_flux'])
        std_int_over_peak = np.nanstd(obsid_cat['int_flux']/obsid_cat['peak_flux'])
        # TODO: add plotting option here to plot the obsid intoverpeak per source or something

        mean_int_over_peak.append(np.nanmedian(int_over_peak))
        std_int_over_peak.append(np.nanstd(int_over_peak))

    logger.debug(f"mean_int_over_peak: {mean_int_over_peak[0]}")
    return mean_int_over_peak, std_int_over_peak


def plot_blur_pernight(obsids, blur, ext="png"):

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obsids)):
        ax.errorbar(obsids[i], np.mean(blur[i]),yerr=np.nanstd(blur[i]),fmt="o",color="C6")
    ax.set_ylabel(f"obsid")
    ax.set_xlabel(f"blur")
    fig.suptitle(f"{night}: blur")
    
    plt.savefig(args.obsids.replace(".txt", f"_intoverpeak.{ext}"), overwrite=True, bbox_inches='tight')

    return 


def plot_intoverpeak_pernight(obslist, int_over_peak, std_int_over_peak, ext="png"):

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    for i in range(len(std_int_over_peak)):
        ax.scatter(int_over_peak[i],std_int_over_peak[i],s=75,color="C6")
    ax.axhline(0.15, color="k", alpha=0.3, linestyle="--")
    ax.axvline(1.15, color="k", alpha=0.3, linestyle="--")
    ax.set_ylabel(f"std(int/peak)")
    ax.set_xlabel(f"mean(int/peak)")
    fig.suptitle(f"{night}: Int/Peak")

    plt.savefig(args.obsids.replace(".txt", f"_intoverpeak.{ext}"), bbox_inches='tight')

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    ax.errorbar(obslist, int_over_peak,yerr=std_int_over_peak, color="C6")
    ax.axhline(np.nanmean(int_over_peak), color="k", alpha=0.3, linestyle="--")
    # ax.axvline(1.15, color="k", alpha=0.3, linestyle="--")
    ax.set_ylabel(f"mean(int/peak)")
    ax.set_xlabel(f"obsid")
    fig.suptitle(f"{night}: Int/Peak")

    plt.savefig(args.obsids.replace(".txt", f"_intoverpeak_perobsid.{ext}"), bbox_inches='tight')

    return 

def plot_intoverpeak_perobs(obsid, int_over_peak, std_int_over_peak, ext="png"):

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    ax.scatter(int_over_peak,std_int_over_peak,fmt="o",color="C6")
    ax.axhline(1.15, color="k", alpha=0.3, linestyle="--")
    ax.axvline(0.15, color="k", alpha=0.3, linestyle="--")
    ax.set_ylabel(f"integrated flux/peak flux")
    ax.set_xlabel(f"std(integrated flux/peak flux)")
    
    fig.suptitle(f"{obsid}: Int/Peak")
    plt.savefig(f"{base_dir}/{obsid}/{obsid}_intoverpeak.{ext}", bbox_inches='tight')
    return 


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to assess the quality of images for obsids and return list of obsids that pass quality assurance to be included in moasics. Note: currently only works on obsids given per channel, not per night. "
    )
    parser.add_argument(
        'obsids',
        type=str,
        help="Path to the .txt file with the new line separated obsids for the given channel, note: currently uses the name of obsid to split by _ and use that as title for night for plotting"
    )
    parser.add_argument(
        '-b',
        '--base_path',
        type=str,
        default=".",
        help="Path to folder containing obsid folders"
    )
    parser.add_argument(
        "--catalogue",
        type=str,
        dest="ref_cat",
        default="GGSM_sparse_unresolved.fits",
        help="Filename of catalogue to use for comparison, assumed ./GGSM_sparse_unresolved from GLEAM-X-pipeline",
    )   


    parser.add_argument(
        '--sep',
        type=float,
        default=1,
        help="Separation from reference catalogue to that of obsid. (default=1arcmin)"
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
        '--flag_high_rms',
        default=True,
        help='Will only select obsids that have RMS no larger than the median RMS at that freq channel + the STD of the RMS of the night'
    )
    parser.add_argument(
        '--flag_bad_io',
        default=True, 
        help='Will run the selection cuts for bad ionosphere on individial obsids based on quality checks from Brandon'
    )



    parser.add_argument(
        '--cut_type',
        help='Can pick either "snr", "suggested" or "both". Suggested is based on two factors, int over peak ratio and int_err vs rms ratio. SNR is based on the SNR as well as the int over peak comparison ratio.',
        default='both',
        type=str
    )
    parser.add_argument(
        '--cut_level_int_peak',
        help='The ratio of int flux to peak flux below which the sources will be retained. (Default=2).',
        default=2,
        type=float
    )
    parser.add_argument(
        '--cut_level_err_int_rms',
        help='The cut level for the err_int and local_rms, below which sources will be kept. Default=2.',
        default=2,
        type=float
    )
    parser.add_argument(
        '--cut_level_snr',
        help='The SNR ratio at which the sources will be cut. Default=5',
        default=5,
        type=float
    )

    parser.add_argument(
        '--plot',
        '-p',
        default=None,
        type=str,
        help="If True, will plot all plots and save to base_dir "
    )
    parser.add_argument(
        '--save_bad_obsids',
        default=False,
        help="If defined, will make a .txt file with the bad obsids"
    )
    parser.add_argument(
        '--save_missing_obsids',
        default=False,
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
    night = args.obsids.split("/")[0].split("_")[2]
    
    make_plots = args.plot
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
        blur = []
        int_over_peak = []
        for obs in obslist: 
            obsid_cat = cut_cat_bright(base_dir,obs)
            
            obs_int_over_peak = obsid_cat["int_flux"]/obsid_cat["peak_flux"]
            err_intoverrms = obsid_cat["err_int_flux"]/obsid_cat["local_rms"]
            snr = obsid_cat["int_flux"]/obsid_cat["local_rms"]
            blur.append(np.log10(obsid_cat['int_flux']/obsid_cat['peak_flux']))
            std_intpeak = np.nanstd(obsid_cat['int_flux']/obsid_cat['peak_flux'])
            std_int_over_peak.append(std_intpeak)
            int_over_peak.append(obs_int_over_peak)
            # TODO: add plotting option here to plot the obsid intoverpeak per source or something

            mean_int_over_peak.append(np.nanmedian(obs_int_over_peak))

        if make_plots is not None: 
            logger.debug(f"Plotting int over peak for night: {night}")
            plot_intoverpeak_pernight(obslist, mean_int_over_peak,std_int_over_peak)
        


            

        logger.debug(f"mean_int_over_peak: {mean_int_over_peak[0]}")
        logger.debug(f"ran cut bright but not crossmatch")
    else: 
        logger.debug(f"Not running the ionosphere analysis")

    # TODO: add plotting option to plot the int/peak+/- std as function of obsid over the night. 
    
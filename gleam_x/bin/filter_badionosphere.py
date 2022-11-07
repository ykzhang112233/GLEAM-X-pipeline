#!/usr/bin/env python

""" Script to assess the quality of obsid images and identify high quality images to include in mosaic step. 
"""


from ast import parse
from gzip import READ
from http.client import NON_AUTHORITATIVE_INFORMATION
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
# import cmasher as cmr
import matplotlib
import numpy.ma as ma

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

def read_obsids(filename):

    try: 
        channel_obsids = np.loadtxt(f"{args.project}/{filename}")
        channel_obsids = ma.array(channel_obsids)

    except FileNotFoundError:
        logger.warning(f"Cannot find txt file with obsids: {args.project}/{filename}")
        channel_obsids = np.array(())
    return channel_obsids


def remove_missing(obsids_chan):

    for i in range(len(obsids_chan)):
        obsid = obsids_chan[i]
        if os.path.exists(f"{base_dir}/{obsid:10.0f}/{obsid:10.0f}_deep-MFS-image-pb_warp_rescaled_comp.fits") is False: 
            logger.debug(f"No catalogue for io checks: {obsid:10.0f}")
            obsids_chan[i] = ma.masked
    
    num_missing = len(obsids_chan) - len(obsids_chan.compressed())

    logger.debug(f"Number of missing obsids: {num_missing}")
    if num_missing>5:
        logger.warning(f"Large number of missing obsids: {num_missing}/{len(obsids_chan)}")

    return obsids_chan


def cut_high_rms(obsids):

    # TODO: currently hardcoding limit for bad std(rms) per channel!!! FIX!
    rms_std_cutoff = [40, 15, 10, 10, 10]
    for i in range(len(obsids)):
        num_postmissing = len(obsids[i].compressed())
        rms_chan = ma.array([np.nan]*len(obsids[i])) 
        obs = obsids[i]
        for j in range(len(obs)):
            if obs[j] is not ma.masked: 
                rmsfile = f"{args.project}/{obsids[i][j]:10.0f}/{obsids[i][j]:10.0f}_deep-MFS-image-pb_warp_rms.fits"
                if os.path.exists(rmsfile):
                    hdu = fits.open(rmsfile)
                    rms_chan[j] = 1.e3*hdu[0].data[int(hdu[0].data.shape[0]/2), int(hdu[0].data.shape[1]/2)] 
                    hdu.close()
                else:
                    logger.debug(f"Found missing obsid while checking RMS, make sure ran check for missing earlier")

        
        logger.debug(f"The std of rms for chan {chans[i]}: {np.nanstd(rms_chan.compressed())}")
        if np.nanstd(rms_chan.compressed()) > rms_std_cutoff[i]:
            logger.warning(f"The std of rms for chan {chans[i]} is HUGE! Flagging all obsids for this channel! You need to come back and inspect everything further: {np.nanstd(rms_chan.compressed())}")
            cutoff = rms_std_cutoff[i]
        else: 
            cutoff = np.nanmean(rms_chan.compressed())+np.nanstd(rms_chan.compressed())

        rms_masked = ma.masked_greater_equal(rms_chan,cutoff) 
        obsids[i][rms_masked.mask] = ma.masked
        frac_flagged = int((1-(len(rms_masked.compressed())/num_postmissing))*100)
        if frac_flagged > 15:
            logger.warning(f"Large number of obsids flagged for high RMS chan {chans[i]}: {frac_flagged}%")

    return obsids

def crossmatch_cats(
    input_cat,
    ref_cat,
    sep=0.5,
):
    ref_cat_file = f"{ref_cat}"
    if os.path.exists(ref_cat_file):
        hdu = fits.open(ref_cat_file)
        refcat = hdu[1].data
        hdu.close()
    else:
        logger.warning(f"Can't find reference NVSS/SUMSS catalogue for xm! ")
        return input_cat

    ref_cat_skycoords = SkyCoord(refcat.RAJ2000,refcat.DEJ2000, frame='fk5',unit=u.deg)
    input_cat_skycoords = SkyCoord(input_cat.ra, input_cat.dec, frame='fk5', unit=u.deg)
    try: 
        idx1,idx2,sep2d,dist3d=search_around_sky(input_cat_skycoords,ref_cat_skycoords,1*u.arcmin)
    except: 
        logger.debug(f"Catalogue has missing or NaN ra/dec. Cutting.")
        clean_mask = (~np.isnan(input_cat.ra) & ~np.isnan(input_cat.dec))
        input_cat_clean = input_cat[clean_mask]
        input_cat_skycoords = SkyCoord(input_cat_clean.ra, input_cat_clean.dec, frame='fk5', unit=u.deg)
        idx1,idx2,sep2d,dist3d=search_around_sky(input_cat_skycoords,ref_cat_skycoords,1*u.arcmin)
    
    num_in_xm = len(idx1)
    frac_flagged = int(((len(idx1)/len(input_cat)))*100)
    if num_in_xm < 500: 
        logger.warning(f"Not many sources in obsid after xm! {num_in_xm} ({frac_flagged}%)")
        raise Exception(f"Too few sources!") 
         

    else: 
        output_cat = input_cat[idx1]
        return output_cat


def check_io(obsids):
    int_over_peak = []
    std_intoverpeak = []
    shape = []
    std_shape = []
    for i in range(len(obsids)):
        obs = obsids[i]
        int_over_peak_chan = ma.array([np.nan]*len(obs))
        std_intoverpeak_chan = ma.array([np.nan]*len(obs)) 
        shape_chan = ma.array([np.nan]*len(obs))
        std_shape_chan = ma.array([np.nan]*len(obs))
        for j in range(len(obs)):
            if obs[j] is not ma.masked: 
                catfile = f"{args.project}/{obs[j]:10.0f}/{obs[j]:10.0f}_deep-MFS-image-pb_warp_rescaled_comp.fits"
                if os.path.exists(catfile):
                    hdu = fits.open(catfile)
                    temp_cat = hdu[1].data
                    hdu.close()
                else:
                    logger.debug(f"Found missing obsid while checking src quality, make sure ran check for missing earlier") 
                    obs[j] = ma.masked
                    int_over_peak_chan[j] = ma.masked
                    std_intoverpeak_chan[j] = ma.masked
                    shape_chan[j] = ma.masked
                    std_shape_chan[j] = ma.masked
                    continue

                try:
                    cat_xm = crossmatch_cats(temp_cat, args.refcat)
                except Exception: 
                    logger.warning(f"Flagging {obs[j]:10.0f} for too few srcs")
                    int_over_peak_chan[j] = ma.masked
                    std_intoverpeak_chan[j] = ma.masked
                    shape_chan[j] = ma.masked
                    std_shape_chan[j] = ma.masked
                    obs[j] = ma.masked
                    continue 
                except: 
                    logger.warning(f"Couldnt do xm!!! Flagging {obs[j]:10.0f}")
                    int_over_peak_chan[j] = ma.masked
                    std_intoverpeak_chan[j] = ma.masked
                    shape_chan[j] = ma.masked
                    std_shape_chan[j] = ma.masked
                    obs[j] = ma.masked
                    continue 

                int_over_peak_obs = ma.array(cat_xm["int_flux"]/cat_xm["peak_flux"])
                err_int_rms_obs = ma.array(cat_xm["err_int_flux"]/cat_xm["local_rms"])
                snr_obs = ma.array(cat_xm["int_flux"]/cat_xm["local_rms"])
                shape_obs = ma.array(cat_xm["a"]/cat_xm["b"])
                num_srcs_precut = len(int_over_peak_obs)

                snr_mask = ma.masked_less(snr_obs,5).mask
                intoverpeak_mask = ma.masked_greater(int_over_peak_obs,2).mask
                err_int_rms_mask = ma.masked_greater(err_int_rms_obs,2).mask

                int_over_peak_obs[snr_mask] = ma.masked
                int_over_peak_obs[intoverpeak_mask] = ma.masked
                int_over_peak_obs[err_int_rms_mask] = ma.masked
                shape_obs[snr_mask] = ma.masked
                shape_obs[intoverpeak_mask] = ma.masked
                shape_obs[err_int_rms_mask] = ma.masked

                int_over_peak_chan[j] = np.nanmean(int_over_peak_obs.compressed())
                std_intoverpeak_chan[j] = np.nanstd(int_over_peak_obs.compressed())
                shape_chan[j] = np.nanmean(shape_obs.compressed())
                std_shape_chan[j] = np.nanstd(shape_obs.compressed())
                num_srcs_postcut = len(int_over_peak_obs.compressed())



                if args.plot == "all":
                    plt_io_obsid(int_over_peak_obs.compressed(), shape_obs.compressed(), f"{obs[j]:10.0f}", color=colors[i+3])
            else: 
                int_over_peak_chan[j] = ma.masked
                std_intoverpeak_chan[j] = ma.masked
                shape_chan[j] = ma.masked
                std_shape_chan[j] = ma.masked

        int_over_peak.append(int_over_peak_chan)
        shape.append(shape_chan)
        std_intoverpeak.append(std_intoverpeak_chan)
        std_shape.append(std_shape_chan)


    return int_over_peak, std_intoverpeak, shape, std_shape

def plt_io_obsid(
    intoverpeak,
    shape,
    obsid,
    ext="png",
    color="C6"
):

    # Just plotting the shape compared to int/flux 
    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    ax.scatter(intoverpeak,np.log(shape),s=50, color=color)
    ax.axhline(0,color="k",alpha=0.3, linestyle="--")

    ax.set_xlabel("int_flux/peak_flux")
    ax.set_ylabel("log(shape) (a/b)")
    fig.suptitle(f"{obsid}: Int/peak vs shape")

    plt.savefig(f"{args.project}/{obsid}/{obsid}_intoverpeak_shape.{ext}", bbox_inches='tight')
    plt.close(fig)
    return 

def plt_io_pernight(
    obsids,
    obsids_nomask,
    intoverpeak,
    intoverpeak_nomask,
    std_intoverpeak,
    std_intoverpeak_nomask,
    shape,
    shape_nomask,
    std_shape,
    std_shape_nomask,
    drift,
    chans,
    ext="png",
    cmap = plt.get_cmap("gnuplot2"),
):
    c_array = np.linspace(0,1,len(obsids)+4)
    colors = cmap(c_array)
    # colors=cmr.take_cmap_colors(
    #     "cmr.flamingo", len(obslist), cmap_range=(0.4, 0.7), return_fmt="hex"
    # )

 
    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    for i in range(len(obsids)):
        obs_chan = obsids[i]
        intoverpeak_chan = intoverpeak[i]
        std_intoverpeak_chan = std_intoverpeak[i]
        obs_chan_all = obsids_nomask[i]
        intoverpeak_chan_all = intoverpeak_nomask[i]
        std_intoverpeak_chan_all = std_intoverpeak_nomask[i]

        ax.errorbar(obs_chan_all, intoverpeak_chan_all,yerr=(std_intoverpeak_chan_all/np.sqrt(len(obs_chan_all))), fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(obs_chan, intoverpeak_chan,yerr=(std_intoverpeak_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3],markeredgecolor="k")
        ax.axhline(np.nanmean(intoverpeak[i]), color=colors[i+3], alpha=0.3, linestyle="--")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"mean(int/peak)")
    ax.set_xlabel(f"obsid")
    ax.legend()
    fig.suptitle(f"{drift}: Int/Peak")
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_intoverpeak_perobs.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obsids)):
        intoverpeak_chan = intoverpeak[i]
        shape_chan = shape[i]
        ax.errorbar(intoverpeak_nomask[i],shape_nomask[i], fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(intoverpeak_chan, shape_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"mean(a/b)")
    ax.set_xlabel(f"mean(int/peak)")
    fig.suptitle(f"{drift}: Int/Peak vs Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_intoverpeak_shape.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obsids)):
        stdshape_chan = std_shape[i]
        shape_chan = shape[i]

        ax.errorbar(shape_nomask[i],std_shape_nomask[i], fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(shape_chan, stdshape_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"std(a/b)")
    ax.set_yscale('log')
    ax.set_xlabel(f"mean(a/b)")
    fig.suptitle(f"{drift}: Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_shape.{ext}", bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obsids)):
        intoverpeak_chan = intoverpeak[i]
        std_intoverpeak_chan = std_intoverpeak[i]

        ax.errorbar(std_intoverpeak_nomask[i],intoverpeak_nomask[i], fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(std_intoverpeak_chan, intoverpeak_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.axhline(1.1, color="k", alpha=0.3, ls="--")
    ax.axvline(0.175, color="k", alpha=0.3, ls="--")
    ax.set_ylabel(f"mean(int/peak)")
    ax.set_xlabel(f"std(int/peak)")
    fig.suptitle(f"{drift}: Int/Peak")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_intoverpeak.{ext}", bbox_inches='tight')
    plt.close(fig)




    return 



if __name__ == "__main__":
    parser = ArgumentParser(
        description="Script to assess the quality of images for obsids and return list of obsids that pass quality assurance to be included in moasics. Note: currently only works on obsids given per channel, not per night. "
    )
    parser.add_argument(
        '--project',
        type=str,
        default=".",
        help="Path to project directory containing obsid folders, also where the drift scan folder is containing text files $project/$drift/*.txt (default= ./)"
    )
    parser.add_argument(
        'obsids',
        type=str,
        help="The text file of obsids to be processed. Will work out if its .txt or cenchan_chan.txt, at most plz $project/drift/drift.txt, no extra directories "
    )
    parser.add_argument(
        '--refcat',
        type=str,
        default="/models/NVSS_SUMSS_psfcal.fits",
        help="reference catalogue to crossmatch and get only bright, unresolved and sparse sources. (default=./models/NVSS_SUMSS_psfcal.fits)"
    )

    parser.add_argument(
        '--flag_high_rms',
        default=True,
        help='Will only select obsids that have RMS no larger than the mean RMS at that freq channel + the STD of the RMS of the night'
    )
    parser.add_argument(
        '--flag_bad_io',
        default=True,
        help="Will run cuts on the quality of sources in each obsid then calculate int/peak etc. to assess io per obsid and over night options (default=True)"
    )
    parser.add_argument(
        '--flag_bad_shape',
        default=False,
        help="Will flag obsids where the shape has a huge scatter (suggesting variable io within image) (default=False)"
    )
    parser.add_argument(
        "--flag_bad_intoverpeak",
        default=False,
        help="Run harsh check of int/peak to flag ones with huge scatter of int/peak (default=False)"
    )
    parser.add_argument(
        "--plot",
        default="min",
        type=str,
        help="Level of plotting to do: all, min, none",
    )

    parser.add_argument(
        '--save_missing_obsids',
        default=None,
        help="If defined, will make a .txt file in directory with all obsids with no *MFS-image-pb_warp_rms.fits file (default=None)"
    )
    parser.add_argument(
        '--save_bad_obsids',
        default=None,
        help="Will make a .txt file with the bad obsids (default=None) "
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
    

    base_dir = args.project
    txtfile = args.obsids 
    refcat=  args.refcat
    logger.debug(f"{txtfile}")
    cmap = plt.get_cmap("gnuplot2")
    c_array = np.linspace(0,1,9)
    colors = cmap(c_array)

    ref_cat_file = f"{args.refcat}"
    if os.path.exists(ref_cat_file):
        do_xm = True
    else:
        logger.warning(f"Can't find reference GGSM catalogue for xm! ")
        do_xm = False

    

    # Reading in the list of obsids: will deal with the cenchan or all based on what the input txt file is called 
    if "cenchan" in txtfile:
        logger.debug(f"Only detected one cenchan, proceeding with just 1")
        obs_txtfile = [txtfile]
        split_string = txtfile.split("/")
        if len(split_string) == 2: 
            split_string = split_string[-1].split("_cenchan_")
            drift = split_string[0]
            chans = [split_string[1].split(".")[0]]
        elif len(split_string)==1:
            split_string = split_string[0].split("_cenchan_")
            chans = [split_string[1].split(".")[0]]
            drift = split_string[0]
        logger.debug(f"drift: {drift}")
        logger.debug(f"cenchan: {chans[0]}")
    else: 
        split_string = txtfile.split("/")
        logger.debug(f"Detected no cenchan, proceeding with allchans")
        if len(split_string) == 2:
            drift = split_string[-1].replace(".txt", "")
        elif len(split_string) == 1:
            drift = split_string[0].replace(".txt","")
        logger.debug(f"drift: {drift}")
        if drift in ["XG_D-27_20201022", "XG_D-27_20201015", "XG_D-27_20201008", "XG_D-27_20201001"]:
            chans = ["69", "93", "121", "145", "169"]
        else: 
            chans = ["069", "093", "121", "145", "169"]
        obs_txtfile = []
        for chan in chans:
            obs_txtfile.append(txtfile.replace(".txt",f"_cenchan_{chan}.txt"))


    # Looking for any missing obsids so they're removed before assessing 
    logger.debug(f"{obs_txtfile}")
    # Reading in the obsids from txt files
    obsids = []
    for i in range(len(chans)):
        cenchan_obsids = read_obsids(obs_txtfile[i])
        obsids_chan = remove_missing(cenchan_obsids)
        obsids.append(obsids_chan)
        logger.debug(f"Number of obsids per in channel {chans[i]}: {len(obsids[i].compressed())}")
        if args.save_missing_obsids is not None: 
            logger.debug(f"Saving missing obsids")
            chan_missing_obsids = obsids_chan[obsids_chan.mask]
            chan_missing_obsids.mask = ma.nomask
            np.savetxt(obs_txtfile[i].replace(".txt", "_missing_obsids.txt"), chan_missing_obsids.compressed(), fmt="%10.0f")


    # Cutting obsids with high RMS in MFS image 
    if args.flag_high_rms is True: 
        obsids = cut_high_rms(obsids)
    else: 
        logger.debug(f"Not running the high RMS cut")

    # Running cut of bad sources to assess io
    # note: first cuts bad srcs, the xm with GGSM to find nice brihgt etc ones before doing the actually assessment 
    drift_intoverpeak, drift_stdintoverpeak, drift_shape, drift_stdshape = check_io(obsids)
    if args.flag_bad_io is True: 
        for i in range(len(chans)):
            color = colors[i+3]
            logger.debug(f"Running the io cut on chan {chans[i]}")
            num_obsids_precut = len(obsids[i].compressed())

            logger.debug(f"STD for drift at chan {chans[i]}: {np.nanstd(drift_intoverpeak[i])}")

            drift_intoverpeak[i] = ma.masked_greater(drift_intoverpeak[i], 1.1) 
            obsids[i][drift_intoverpeak[i].mask] = ma.masked
            drift_intoverpeak[i][drift_intoverpeak[i].mask] = ma.masked
            drift_stdintoverpeak[i][drift_intoverpeak[i].mask] = ma.masked
            drift_shape[i][drift_intoverpeak[i].mask] = ma.masked
            drift_stdshape[i][drift_intoverpeak[i].mask] = ma.masked

            drift_stdintoverpeak[i] = ma.masked_greater(drift_stdintoverpeak[i],0.175)
            obsids[i][drift_stdintoverpeak[i].mask] = ma.masked
            drift_intoverpeak[i][drift_stdintoverpeak[i].mask] = ma.masked
            drift_stdintoverpeak[i][drift_stdintoverpeak[i].mask] = ma.masked
            drift_shape[i][drift_stdintoverpeak[i].mask] = ma.masked
            drift_stdshape[i][drift_stdintoverpeak[i].mask] = ma.masked


            num_obsids_postio = len(obsids[i].compressed())

            logger.debug(f"Obsids for chan {chans[i]} before any cuts: {len(obsids[i])}")
            logger.debug(f"Obsids for chan {chans[i]} after io cuts: {len(obsids[i].compressed())}")



            frac_flagged = int((1-(num_obsids_postio/num_obsids_precut))*100)
            logger.debug(f"Number of obsids with *good* io for {chans[i]}: {num_obsids_postio} ({frac_flagged}%)")
            if frac_flagged > 20:
                logger.warning(f"Large number of obsids flagged for bad io at chan {chans[i]}!: {frac_flagged}%")

    if args.flag_bad_shape is True: 
        for i in range(len(chans)):
            # Now cutting obsids where shape has huge scatter too:
            # TODO: Also doing same style cut as RMS with harsh level, change
            harsh_stdshape_cutoff = [0.1, 0.01, 0.01, 0.01, 0.01]
            harsh_shape_cutoff = [0.3, 0.3, 0.2, 0.2, 0.2]
            if np.nanstd(drift_stdshape[i].compressed()) > harsh_stdshape_cutoff[i]:
                logger.warning(f"HUGE shape scatter for chan {chans[i]}: {np.nanstd(drift_stdshape[i].compressed())}")
                logger.warning(f"mean of std(shape) for chan {chans[i]}: {np.nanmean(drift_stdshape[i].compressed())}")
                shape_cutoff = harsh_shape_cutoff[i] + (3*harsh_stdshape_cutoff[i])
            else: 
                logger.debug(f"std of std(shape) for chan {chans[i]}: {np.nanstd(drift_stdshape[i].compressed())}")
                shape_cutoff = np.nanmean(drift_stdshape[i].compressed())+(3*np.nanstd(drift_stdshape[i].compressed()))

            num_preshapecut = len(drift_stdshape[i].compressed())
            shape_cut = ma.masked_greater(drift_stdshape[i], shape_cutoff)
            shape_mask = shape_cut.mask
            num_postshapecut = len(shape_cut.compressed())


            logger.debug(f"Cutting {num_preshapecut - num_postshapecut} based on big std of shape!!")
            frac_flagged = int((1-(num_postshapecut/num_preshapecut))*100)
            if frac_flagged > 10:
                logger.warning(f"Large number of obsids flagged for bad shape in {chans[i]}!: {num_preshapecut - num_postshapecut} ({frac_flagged})%")         

            obsids[i][shape_mask] = ma.masked
            drift_intoverpeak[i][shape_mask] = ma.masked
            drift_stdintoverpeak[i][shape_mask] = ma.masked
            drift_shape[i][shape_mask] = ma.masked
            drift_stdshape[i][shape_mask] = ma.masked

    if args.flag_bad_intoverpeak is True: 
        for i in range(len(chans)):
            # Now cutting obsids where std(int/peak) has huge scatter too:
            intoverpeak_cutoff = np.nanmean(drift_stdintoverpeak[i].compressed())+(3*np.nanstd(drift_stdintoverpeak[i].compressed()))


            num_preintcut = len(drift_stdintoverpeak[i].compressed())
            intoverpeak_cut = ma.masked_greater(drift_stdintoverpeak[i], intoverpeak_cutoff)
            intoverpeak_mask = intoverpeak_cut.mask
            num_postintcut = len(intoverpeak_cut.compressed())


            logger.debug(f"Cutting {num_preintcut - num_postintcut} based on big std of shape!!")
            frac_flagged = int((1-(num_postintcut/num_preintcut))*100)
            if frac_flagged > 10:
                logger.warning(f"Large number of obsids flagged for bad int/peak in {chans[i]}!: {num_preshapecut - num_postshapecut} ({frac_flagged})%")         

            obsids[i][intoverpeak_mask] = ma.masked
            drift_intoverpeak[i][intoverpeak_mask] = ma.masked
            drift_stdintoverpeak[i][intoverpeak_mask] = ma.masked
            drift_shape[i][intoverpeak_mask] = ma.masked
            drift_stdshape[i][intoverpeak_mask] = ma.masked


    if args.plot in ["all", "min"]:   
        logger.debug(f"Plotting for drift")     
        obsids_nomask = []
        intoverpeak_nomask = []
        std_intoverpeak_nomask = []
        shape_nomask = []
        std_shape_nomask = []
        obsids_good = []
        intoverpeak_good = []
        std_intoverpeak_good = []
        shape_good = []
        std_shape_good = []
        for i in range(len(obsids)):
            
            logger.debug(f"Making the non masked arrays")
            obs_chan_nomask = obsids[i]
            obsids_good.append(obsids[i].compressed())
            obs_chan_nomask.mask = ma.nomask 
            obsids_nomask.append(obs_chan_nomask)
            
            intoverpeak_nomask_chan = drift_intoverpeak[i]
            intoverpeak_good.append(drift_intoverpeak[i].compressed())
            intoverpeak_nomask_chan.mask = ma.nomask 
            intoverpeak_nomask.append(intoverpeak_nomask_chan)

            std_intoverpeak_nomask_chan = drift_stdintoverpeak[i]
            std_intoverpeak_good.append(drift_stdintoverpeak[i].compressed())
            std_intoverpeak_nomask_chan.mask = ma.nomask 
            std_intoverpeak_nomask.append(std_intoverpeak_nomask_chan)

            shape_nomask_chan = drift_shape[i]
            shape_good.append(drift_shape[i].compressed())
            shape_nomask_chan.mask = ma.nomask 
            shape_nomask.append(shape_nomask_chan)

            std_shape_nomask_chan = drift_stdshape[i]
            std_shape_good.append(drift_stdshape[i].compressed())
            std_shape_nomask_chan.mask = ma.nomask 
            std_shape_nomask.append(std_shape_nomask_chan)

            logger.debug(f"Num obsids per chan NO CUT: {len(obsids_nomask[i].compressed())}")

            total_frac_flagged = int((1-(len(obsids_good[i])/len(obsids_nomask[i].compressed())))*100)
            logger.warning(f"Num obsids per chan aftercut/beforecut: {len(obsids_good[i])}/{len(obsids_nomask[i].compressed())} ({total_frac_flagged}%)")
        plt_io_pernight(obsids_good, obsids_nomask, intoverpeak_good, intoverpeak_nomask, std_intoverpeak_good, std_intoverpeak_nomask, shape_good, shape_nomask, std_shape_good, std_shape_nomask, drift, chans)
        
    if args.save_bad_obsids is not None: 
        logger.debug(f"Saving missing obsids")
        for i in range(len(obsids)):
            obsids_chan = obsids[i]
            chan_bad_obsids = obsids_chan[obsids_chan.mask]
            chan_bad_obsids.mask = ma.nomask
            np.savetxt(obs_txtfile[i].replace(".txt", "_bad_obsids.txt"), chan_bad_obsids.compressed(), fmt="%10.0f")
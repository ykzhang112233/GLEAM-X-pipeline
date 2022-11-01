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


def remove_missing(obsids):

    for i in range(len(obsids)):
        obsid = obsids[i]
        if os.path.exists(f"{base_dir}/{obsid:10.0f}/{obsid:10.0f}_deep-MFS-image-pb_warp_rescaled_comp.fits") is False: 
            logger.debug(f"No catalogue for io checks: {obsid:10.0f}")
            obsids[i] = ma.masked
    
    num_missing = len(obsids) - len(obsids.compressed())

    logger.debug(f"Number of missing obsids: {num_missing}")
    if num_missing>5:
        logger.warning(f"Large number of missing obsids: {num_missing}/{len(obsids.compressed())}")

    return obsids


def cut_high_rms(obsids):

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

        cutoff = np.nanmedian(rms_chan)+np.nanstd(rms_chan)

        rms_masked = ma.masked_greater_equal(rms_chan,cutoff) 

        obsids[i][rms_masked.mask]

        frac_flagged = int((1-(len(obsids[i].compressed())/num_postmissing))*100)
        if frac_flagged > 15:
            logger.warning(f"Large number of obsids flagged for high RMS: {frac_flagged}%")

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

    output_cat = input_cat[idx1]

    num_in_xm = len(idx1)
    frac_flagged = int((1-(len(idx1)/len(input_cat)))*100)
    if num_in_xm < 500: 
        logger.debug(f"Not many sources in obsid after xm! {num_in_xm} ({frac_flagged}%)")
    else: 
        return output_cat


    return output_cat 

def check_io(obsid, do_xm, color):

    # for i in range(len(obsids_list)):
    catfile = f"{args.project}/{obsid:10.0f}/{obsid:10.0f}_deep-MFS-image-pb_warp_rescaled_comp.fits"
    if os.path.exists(catfile):
        hdu = fits.open(catfile)
        temp_cat = hdu[1].data
        hdu.close()
    else:
        logger.debug(f"Found missing obsid while checking src quality, make sure ran check for missing earlier")
        return 
    
    
    if do_xm == True:
        # logger.debug("Doing XM!")
        cat_xm = crossmatch_cats(temp_cat, args.refcat)
        logger.debug(f"Found {len(cat_xm)} sources for {obsid:10.0f}")
    else:
        logger.debug("NOT RUNNING XM!!")
        cat_xm = temp_cat
        # return , [np.nanmean(int_over_peak[mask]), np.nanstd(int_over_peak[mask])], [np.nanmean(shape[mask]),np.nanstd(shape[mask])]
    
        
    int_over_peak = cat_xm["int_flux"]/cat_xm["peak_flux"]
    err_intoverrms = cat_xm["err_int_flux"]/cat_xm["local_rms"]
    snr = cat_xm["int_flux"]/cat_xm["local_rms"]   
    shape = cat_xm["a"]/cat_xm["b"]

    mask = np.where((int_over_peak<=2)&(err_intoverrms<=2)&(snr>=5))
    cat = cat_xm[mask]

    if args.plot == "all":
        plt_io_obsid(int_over_peak[mask], shape[mask], f"{obsid:10.0f}", color=color)

    return cat, [np.nanmedian(int_over_peak[mask]), np.nanstd(int_over_peak[mask])], [np.nanmedian(shape[mask]),np.nanstd(shape[mask])]



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
    obslist,
    good_inds,
    intoverpeak,
    std_intoverpeak,
    shape,
    std_shape,
    drift,
    chans,
    ext="png",
    cmap = plt.get_cmap("gnuplot2"),
):



    c_array = np.linspace(0,1,len(obslist)+4)
    colors = cmap(c_array)
    # colors=cmr.take_cmap_colors(
    #     "cmr.flamingo", len(obslist), cmap_range=(0.4, 0.7), return_fmt="hex"
    # )
    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obslist)):
        obs_chan = obslist[i]
        intoverpeak_chan = intoverpeak[i]
        std_intoverpeak_chan = std_intoverpeak[i]
        good_inds_chan = good_inds[i]
        ax.errorbar(obs_chan, intoverpeak_chan,yerr=(std_intoverpeak_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(obs_chan[good_inds_chan], intoverpeak_chan[good_inds_chan],yerr=(std_intoverpeak_chan[good_inds_chan]/np.sqrt(len(obs_chan[good_inds_chan]))), fmt="o", color=colors[i+3],markeredgecolor="k")
        ax.axhline(np.nanmedian(intoverpeak[i]), color=colors[i+3], alpha=0.3, linestyle="--")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Good IO")
    ax.set_ylabel(f"median(int/peak)")
    ax.set_xlabel(f"obsid")
    ax.legend()
    fig.suptitle(f"{drift}: Int/Peak")
    plt.savefig(f"{args.project}/{drift}/{drift}_intoverpeak_perobs.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obslist)):
        chan_intoverpeak = intoverpeak[i]
        chan_shape = shape[i]
        good_inds_chan = good_inds[i]
        ax.errorbar(chan_intoverpeak,chan_shape, fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(chan_intoverpeak[good_inds_chan], chan_shape[good_inds_chan], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Good IO")
    ax.set_ylabel(f"median(a/b)")
    ax.set_xlabel(f"median(int/peak)")
    fig.suptitle(f"{drift}: Int/Peak vs Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_intoverpeak_shape.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obslist)):
        chan_stdshape = std_shape[i]
        chan_shape = shape[i]
        good_inds_chan = good_inds[i]
        ax.errorbar(chan_shape,chan_stdshape, fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(chan_shape[good_inds_chan], chan_stdshape[good_inds_chan], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Good IO")
    ax.set_ylabel(f"std(a/b)")
    ax.set_yscale('log')
    ax.set_xlabel(f"median(a/b)")
    fig.suptitle(f"{drift}: Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_shape.{ext}", bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(obslist)):
        chan_intoverpeak = intoverpeak[i]
        chan_std_intoverpeak = std_intoverpeak[i]
        good_inds_chan = good_inds[i]
        ax.errorbar(chan_intoverpeak,chan_std_intoverpeak, fmt="o", color=colors[i+3], label=chans[i])
        ax.errorbar(chan_intoverpeak[good_inds_chan], chan_std_intoverpeak[good_inds_chan], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Good IO")
    ax.set_ylabel(f"std(int/peak)")
    ax.set_xlabel(f"median(int/peak)")
    fig.suptitle(f"{drift}: Int/Peak")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_intoverpeak.{ext}", bbox_inches='tight')
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
        help='Will only select obsids that have RMS no larger than the median RMS at that freq channel + the STD of the RMS of the night'
    )
    parser.add_argument(
        '--flag_bad_io',
        default="drift",
        type=str,
        help="Will run cuts on the quality of sources in each obsid then calculate int/peak etc. to assess io per obsid and over night options [drift, obsid] (default=drift)"
    )
    parser.add_argument(
        "--plot",
        default="all",
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

        if args.save_missing_obsids is not None: 
            logger.debug(f"Saving missing obsids")
            chan_missing_obsids = obsids_chan[obsids_chan.mask]
            chan_missing_obsids.mask = ma.nomask
            np.savetxt(obs_txtfile[i].replace(".txt", "_missing_obsids.txt"), chan_missing_obsids.compressed(), fmt="%10.0f")

    
    # Cutting obsids with high RMS in MFS image 
    if args.flag_high_rms is True: 
        for i in range(len(chans)):
            logger.debug(f"Number of obsids per in channel {chans[i]}: {len(obsids[i].compressed())}")
            obsids = cut_high_rms(obsids)
    else: 
        logger.debug(f"Not running the high RMS cut")

    # Running cut of bad sources to assess io
    # note: first cuts bad srcs, the xm with GGSM to find nice brihgt etc ones before doing the actually assessment 
    if args.flag_bad_io in ["drift", "obsid"]: 
        drift_meanintoverpeak = []
        drift_stdintoverpeak = []
        drift_meanshape = []
        drift_stdshape = []

        drift_cut_obsids = []
        drift_cut_meanintoverpeak = []
        drift_cut_stdintoverpeak = []
        drift_cut_meanshape = []
        drift_cut_stdshape = []
        for i in range(len(chans)):
            color = colors[i+3]
            logger.debug(f"Running the io check on chan {chans[i]}")
            obslist_iochecks = obsids_postrms[i]
            num_obsids_precut = len(obslist_iochecks)
            obs_meanintoverpeak = []
            obs_stdintoverpeak = []
            obs_meanshape = []
            obs_stdshape = []
            for o in range(len(obslist_iochecks)):
                cat_xm, obs_intoverpeak_temp, obs_shape_temp = check_io(obslist_iochecks[o], do_xm, color)
                obs_meanintoverpeak.append(obs_intoverpeak_temp[0])
                obs_stdintoverpeak.append(obs_intoverpeak_temp[1])
                obs_meanshape.append(obs_shape_temp[0])
                obs_stdshape.append(obs_shape_temp[1])

            obs_meanintoverpeak = np.array(obs_meanintoverpeak)
            obs_stdintoverpeak = np.array(obs_stdintoverpeak)
            obs_meanshape = np.array(obs_meanshape)
            obs_stdshape = np.array(obs_stdshape)

            drift_meanintoverpeak.append(obs_meanintoverpeak)
            drift_stdintoverpeak.append(obs_stdintoverpeak)
            drift_meanshape.append(obs_meanshape)
            drift_stdshape.append(obs_stdshape)
            logger.debug(f"STD for drift at chan {chans[i]}: {np.nanstd(obs_meanintoverpeak)}")


            if np.nanstd(obs_meanintoverpeak)>=0.015:
                logger.warning(f"Io seems bad for chan{chans[i]}!! Cutting!: std = {np.nanstd(obs_meanintoverpeak)}")
                if args.flag_bad_io == "drift":

                    drift_cutoff_max = np.nanmean(obs_meanintoverpeak)+np.nanstd(obs_meanintoverpeak)
                    drift_cutoff_min = np.nanmean(obs_meanintoverpeak)-np.nanstd(obs_meanintoverpeak)
                    mask = np.where((obs_meanintoverpeak<drift_cutoff_max)&(obs_meanintoverpeak>drift_cutoff_min))
                elif args.flag_bad_io == "obsid":
                    obsid_cutoff = 0.015
                    mask = obs_stdintoverpeak>obsid_cutoff


                obslist_iocut = obslist_iochecks[mask]
                cut_meanintoverpeak = obs_meanintoverpeak[mask]
                cut_stdintoverpeak = obs_stdintoverpeak[mask]
                cut_meanshape = obs_meanshape[mask]
                cut_stdshape = obs_stdshape[mask]
                
            else: 
                logger.debug(f"std of int_over_peak for drift not too high, not cutting any obsids for chan: {chans[i]}")
                obslist_iocut = obslist_iochecks
                cut_meanintoverpeak = obs_meanintoverpeak
                cut_stdintoverpeak = obs_stdintoverpeak
                cut_meanshape = obs_meanshape
                cut_stdshape = obs_stdshape
            


            mask2 = np.where(cut_stdshape < (np.nanmedian(cut_stdshape)+np.nanstd(cut_stdshape)))

            obslist_iocut = obslist_iocut[mask2]
            num_obsids_postcut = len(obslist_iocut)
            cut_meanintoverpeak = cut_meanintoverpeak[mask2].tolist()
            cut_stdintoverpeak = cut_stdintoverpeak[mask2].tolist()
            cut_meanshape = cut_meanshape[mask2].tolist()
            cut_stdshape = cut_stdshape[mask2].tolist()

            drift_cut_obsids.append(obslist_iocut)
            drift_cut_meanintoverpeak.append(cut_meanintoverpeak)
            drift_cut_stdintoverpeak.append(cut_stdintoverpeak)
            drift_cut_meanshape.append(cut_meanshape)
            drift_cut_stdshape.append(cut_stdshape)

            logger.warning(f"Obsids for chan {chans[i]} before io cuts: {len(obslist_iochecks)}")
            logger.warning(f"Obsids for chan {chans[i]} after io cuts: {len(obslist_iocut)}")



            frac_flagged = int((1-(num_obsids_postcut/num_obsids_precut))*100)
            logger.debug(f"Number of obsids with *good* io for {chans[i]}: {len(obslist_iocut)} ({frac_flagged}%)")
            if frac_flagged > 15:
                logger.warning(f"Large number of obsids flagged for bad io!: {frac_flagged}%")

        if args.plot in ["all", "min"]:
            drift_good_inds = []
            for i in range(len(chans)):
                inds_good_io = np.where(np.in1d(drift_cut_obsids[i],obsids_postrms[i]))
                drift_good_inds.append(np.squeeze(inds_good_io))
            
            plt_io_pernight(obsids_postrms, drift_good_inds, drift_meanintoverpeak, drift_stdintoverpeak, drift_meanshape,drift_stdshape, drift,chans)


            # plt_io_pernight(drift_cut_obsids, drift_cut_meanintoverpeak, drift_cut_stdintoverpeak, drift_cut_meanshape,drift_cut_stdshape,drift,chans,name=f"postcut{args.flag_bad_io}")
        
    if args.save_bad_obsids is not None: 
        logger.debug(f"Saving missing obsids")
        for i in range(len(obsids)):
            obsids_chan = obsids[i]
            chan_bad_obsids = obsids_chan[obsids_chan.mask]
            chan_bad_obsids.mask = ma.nomask
            np.savetxt(obs_txtfile[i].replace(".txt", "_bad_obsids.txt"), chan_bad_obsids.compressed(), fmt="%10.0f")
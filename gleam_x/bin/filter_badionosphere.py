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
from astropy.table import Table, Column

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


def remove_missing(obsids_chan, extra = ""):

    for i in range(len(obsids_chan)):
        obsid = obsids_chan[i] 
        if os.path.exists(f"{base_dir}/{obsid:10.0f}/{obsid:10.0f}_deep-MFS-image-pb_warp_rescaled_comp{extra}.fits") is False: 
            logger.warning(f"No catalogue for io checks: {obsid:10.0f}")
            obsids_chan[i] = ma.masked
    
    num_missing = ma.count_masked(obsids_chan)
    missing_mask = obsids_chan.mask

    logger.debug(f"Number of missing obsids: {num_missing}")
    if num_missing>5:
        logger.warning(f"Large number of missing obsids: {num_missing}/{len(obsids_chan)}")

    return obsids_chan, missing_mask 


def cut_high_rms(obsids, extra = ""):

    # TODO: currently hardcoding limit for bad std(rms) per channel!!! FIX!
    rms_std_cutoff = [40, 15, 10, 10, 10]
    rms_mask = []
    for i in range(len(obsids)):
        num_postmissing = ma.count(obsids[i])
        total = ma.size(obsids[i])
        if total == num_postmissing is True: 
            logger.warning(f"No missing for this channel.. is that right? ")
        rms_chan = ma.array([np.nan]*len(obsids[i]), mask=obsids[i].mask) 
        obs = obsids[i]
        for j in range(len(obs)):
            if obs[j] is not ma.masked: 
                rmsfile = f"{args.project}/{obsids[i][j]:10.0f}/{obsids[i][j]:10.0f}_deep-MFS-image-pb_warp_rms{extra}.fits"
                if os.path.exists(rmsfile):
                    hdu = fits.open(rmsfile)
                    rms_chan[j] = 1.e3*hdu[0].data[int(hdu[0].data.shape[0]/2), int(hdu[0].data.shape[1]/2)] 
                    hdu.close()
                else:
                    logger.warning(f"Found missing obsid while checking RMS, make sure ran check for missing earlier: {obsids[i][j]:10.0f}")

        
        logger.debug(f"The std of rms for chan {chans[i]}: {np.nanstd(rms_chan.compressed())}")
        if np.nanstd(rms_chan.compressed()) > rms_std_cutoff[i]:
            logger.warning(f"The std of rms for chan {chans[i]} is HUGE! Flagging obsids for this channel using a hard cut, come back and inspect: {np.nanstd(rms_chan.compressed())}")
            cutoff = rms_std_cutoff[i]
        else: 
            cutoff = np.nanmean(rms_chan.compressed())+np.nanstd(rms_chan.compressed())

        rms_masked = ma.masked_greater_equal(rms_chan,cutoff)
        # obs[rms_masked.mask] = ma.masked
        rms_mask.append(rms_masked)
        frac_flagged = int((1-(len(rms_masked.compressed())/num_postmissing))*100)
        if frac_flagged > 15:
            logger.warning(f"Large number of obsids flagged for high RMS chan {chans[i]}: {frac_flagged}%")

    return rms_mask

def crossmatch_cats(
    input_cat,
    ref_cat,
    sep=1,
):


    ref_cat_skycoords = SkyCoord(ref_cat.RAJ2000,ref_cat.DEJ2000, frame='fk5',unit=u.deg)
    input_cat_skycoords = SkyCoord(input_cat.ra, input_cat.dec, frame='fk5', unit=u.deg)
    try: 
        idx1,idx2,sep2d,dist3d=search_around_sky(input_cat_skycoords,ref_cat_skycoords,sep*u.arcmin)
    except: 
        logger.debug(f"Catalogue has missing or NaN ra/dec. Cutting.")
        clean_mask = (~np.isnan(input_cat.ra) & ~np.isnan(input_cat.dec))
        input_cat_clean = input_cat[clean_mask]
        input_cat_skycoords = SkyCoord(input_cat_clean.ra, input_cat_clean.dec, frame='fk5', unit=u.deg)
        idx1,idx2,sep2d,dist3d=search_around_sky(input_cat_skycoords,ref_cat_skycoords,sep*u.arcmin)
    
    
    iso_nvss_sumss = idx1[np.unique(idx2,return_index=True)[1]]
    idx1_iso = np.unique(iso_nvss_sumss)
    # Adding only unique ones to make sure isolated 
    num_in_xm = len(idx1_iso)
    frac_flagged = int(((len(idx1_iso)/len(input_cat)))*100)
    if num_in_xm < 200: 
        logger.warning(f"Not many sources in obsid after xm! {num_in_xm} ({frac_flagged}%)")
        raise Exception(f"Too few sources!") 
         

    else: 
        # output_cat = input_cat[idx1]
        return idx1_iso


def check_io(obsids, missing_mask, xm_cat, extra = ""):
    int_over_peak = []
    std_intoverpeak = []
    shape = []
    std_shape = []
    for i in range(len(obsids)):
        obs = ma.array(obsids[i].data, mask=missing_mask[i])
        int_over_peak_chan = ma.array([np.nan]*len(obsids[i].data), mask=missing_mask[i])
        std_intoverpeak_chan = ma.array([np.nan]*len(obsids[i].data), mask=missing_mask[i])
        shape_chan = ma.array([np.nan]*len(obsids[i].data), mask=missing_mask[i])
        std_shape_chan = ma.array([np.nan]*len(obsids[i].data), mask=missing_mask[i])
        for j in range(len(obs)):
            if obs[j] is not ma.masked: 
                catfile = f"{args.project}/{obs[j]:10.0f}/{obs[j]:10.0f}_deep-MFS-image-pb_warp_rescaled_comp{extra}.fits"
                savefile = f"{args.project}/{obs[j]:10.0f}/{obs[j]:10.0f}_iocheck_comp{extra}.csv"
                if os.path.exists(catfile):
                    try: 
                        hdu=fits.open(catfile)
                        temp_cat = hdu[1].data
                    except:
                        hdu_temp = Table.read(catfile, format="votable")
                        hdu_temp.write(catfile, format="fits", overwrite=True)
                        hdu = fits.open(catfile)
                        temp_cat = hdu[1].data
                    save_cat = Table.read(catfile, format="fits")
                    uuid_save = save_cat["uuid"]
                    # hdu.close()
                    
                else:
                    logger.warning(f"Found missing obsid while checking src quality, make sure ran check for missing earlier") 
                    obs[j] = ma.masked
                    int_over_peak_chan[j] = ma.masked
                    std_intoverpeak_chan[j] = ma.masked
                    shape_chan[j] = ma.masked
                    std_shape_chan[j] = ma.masked
                    continue

                try:
                    idx1 = crossmatch_cats(temp_cat, xm_cat)
                    cat_xm = temp_cat[idx1]

                    
                    
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
                shape_obs = ma.array(cat_xm["a"]*cat_xm["b"])
                psf_obs = ma.array(cat_xm["psf_a"]*cat_xm["psf_b"])
                uuid_xm = ma.array(cat_xm["uuid"])
                num_srcs_precut = len(int_over_peak_obs)
                background_obs = ma.array(cat_xm["background"])
                rms_obs = ma.array(cat_xm["local_rms"])


                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="xm_srcs", data=srcs)
                    save_cat.add_column(io_srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("xm_srcs", srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)


                snr_mask = ma.masked_less(snr_obs,50).mask
                intoverpeak_mask = ma.masked_greater(int_over_peak_obs,2).mask
                err_int_rms_mask = ma.masked_greater(err_int_rms_obs,2).mask
                background_cut = ma.masked_greater(background_obs,np.nanmean(background_obs.data)+np.nanstd(background_obs.data)).mask
                rms_cut = ma.masked_greater(rms_obs,np.nanmean(rms_obs.data)+np.nanstd(rms_obs.data)).mask

                int_over_peak_obs[snr_mask] = ma.masked
                shape_obs[snr_mask] = ma.masked
                uuid_xm[snr_mask] = ma.masked

                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="io_snr_cut", data=srcs)
                    save_cat.add_column(io_srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("io_snr_cut", srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)

                int_over_peak_obs[intoverpeak_mask] = ma.masked
                uuid_xm[intoverpeak_mask] = ma.masked
                shape_obs[intoverpeak_mask] = ma.masked

                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="io_int_cut", data=srcs)
                    save_cat.add_column(io_srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("io_int_cut", srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)


                shape_obs[err_int_rms_mask] = ma.masked
                int_over_peak_obs[err_int_rms_mask] = ma.masked
                uuid_xm[err_int_rms_mask] = ma.masked

                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                num_srcs_postcut = len(uuid_xm.compressed())
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="io_errint_cut", data=srcs)
                    save_cat.add_column(io_srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("io_errint_cut", srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)    


                shape_obs[background_cut] = ma.masked
                int_over_peak_obs[background_cut] = ma.masked
                uuid_xm[background_cut] = ma.masked

                shape_obs[rms_cut] = ma.masked
                int_over_peak_obs[rms_cut] = ma.masked
                uuid_xm[rms_cut] = ma.masked

                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                num_srcs_postcut = len(uuid_xm.compressed())
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="io_rms_cut", data=srcs)
                    save_cat.add_column(io_srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("io_rms_cut", srcs)
                    # save_cat.write(catfile,format="fits",overwrite=True)    


                int_over_peak_chan[j] = np.nanmean(int_over_peak_obs.data)
                std_intoverpeak_chan[j] = np.nanstd(int_over_peak_obs.data)
                shape_chan[j] = np.nanmean(shape_obs)
                std_shape_chan[j] = np.nanstd(shape_obs.data)

                if args.plot == "all":
                    plt_io_obsid(int_over_peak_obs.compressed(), shape_obs.compressed(), f"{obs[j]:10.0f}", color=colors[i+3])
                
                frac_srcs_flagged = int((num_srcs_postcut/num_srcs_precut)*100)
                if num_srcs_postcut<200:
                    logger.warning(f"Only {num_srcs_postcut} srcs in field: flagging {obs[j]:10.0f}")
                    int_over_peak_chan[j] = ma.masked
                    std_intoverpeak_chan[j] = ma.masked
                    shape_chan[j] = ma.masked
                    std_shape_chan[j] = ma.masked


                srcs = np.zeros(len(save_cat))
                # srcs[idx1] = 1
                good_uuids = uuid_xm.compressed()
                for i in range(len(good_uuids)):
                    if good_uuids[i] in uuid_save:
                        index = np.where(uuid_save == good_uuids[i])
                        srcs[index] = 1                
                try:
                    io_srcs = Column(name="io_srcs", data=srcs)
                    save_cat.add_column(io_srcs)
                    save_cat.write(savefile,format="csv",overwrite=True)
                except: 
                    logger.debug(f"Already have cat, overwriting")
                    save_cat.replace_column("io_srcs", srcs)
                    save_cat.write(savefile,format="csv",overwrite=True)


                    

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
    intoverpeak,
    std_intoverpeak,
    shape,
    std_shape,
    mask,
    drift,
    chans,
    ext="png",
    cmap = plt.get_cmap("gnuplot2"),
    comparison = False,
):
    c_array = np.linspace(0,1,len(chans)+4)
    colors = cmap(c_array)
    # colors=cmr.take_cmap_colors(
    #     "cmr.flamingo", len(obslist), cmap_range=(0.4, 0.7), return_fmt="hex"
    # )

 
    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    for i in range(len(chans)):
        # if comparison is False:
        #     obs_chan = obsids[0]
        #     chan_mask = mask[0]
        # else:
        obs_chan = obsids[i].data
        chan_mask = mask[i]
        intoverpeak_chan = intoverpeak[i].data
        std_intoverpeak_chan = std_intoverpeak[i].data
        logger.debug(f"plotting for {chans[i]}")
            

        ax.errorbar(obs_chan, intoverpeak_chan,yerr=(std_intoverpeak_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3], label=chans[i])
        # if chan_mask.all() == False:
        #     ax.errorbar(obs_chan, intoverpeak_chan,yerr=(std_intoverpeak_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3],markeredgecolor="k")
        # else: 
        ax.errorbar(obs_chan[~chan_mask], intoverpeak_chan[~chan_mask],yerr=(std_intoverpeak_chan[~chan_mask]/np.sqrt(len(obs_chan[~chan_mask]))), fmt="o", color=colors[i+3],markeredgecolor="k")
        ax.axhline(np.nanmean(intoverpeak_chan), color=colors[i+3], alpha=0.3, linestyle="--")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"mean(int/peak)")
    ax.axhline(args.intoverpeak_cut, color="k", alpha=0.3, ls="--")
    ax.set_xlabel(f"obsid")
    ax.legend()
    fig.suptitle(f"{drift}: Int/Peak")
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_intoverpeak_perobs.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)

    for i in range(len(chans)):
        # if comparison is False:
        #     obs_chan = obsids[0].data
        #     chan_mask = mask[0]
        # else:
        obs_chan = obsids[i].data
        chan_mask = mask[i]
        shape_chan = shape[i].data
        std_shape_chan = std_shape[i].data

        ax.errorbar(obs_chan, shape_chan,yerr=(std_shape_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3], label=chans[i])
        # if chan_mask.all() == False:
        #     ax.errorbar(obs_chan, shape_chan,yerr=(std_shape_chan/np.sqrt(len(obs_chan))), fmt="o", color=colors[i+3],markeredgecolor="k")
        # else:
        ax.errorbar(obs_chan[~chan_mask], shape_chan[~chan_mask],yerr=(std_shape_chan[~chan_mask]/np.sqrt(len(obs_chan[~chan_mask]))), fmt="o", color=colors[i+3],markeredgecolor="k")
        ax.axhline(np.nanmean(shape_chan), color=colors[i+3], alpha=0.3, linestyle="--")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"mean(a*b/psf)")
    ax.set_xlabel(f"obsid")
    ax.legend()
    fig.suptitle(f"{drift}: Shape")
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_shape_perobs.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(chans)):
        intoverpeak_chan = intoverpeak[i].data
        shape_chan = std_shape[i].data
        # if args.comparison is False: 
        #     chan_mask = mask[0]
        # else:
        chan_mask = mask[i]
        ax.errorbar(intoverpeak_chan,shape_chan, fmt="o", color=colors[i+3], label=chans[i])
        # if chan_mask.all() == False:
        #     ax.errorbar(intoverpeak_chan, shape_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
        # else: 
        ax.errorbar(intoverpeak_chan[~chan_mask], shape_chan[~chan_mask], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"std(a*b/psf)")
    ax.set_xlabel(f"mean(int/peak)")
    ax.axvline(args.intoverpeak_cut, color="k", alpha=0.3, ls="--")
    fig.suptitle(f"{drift}: Int/Peak vs Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_intoverpeak_shape.{ext}", bbox_inches='tight')
    plt.close(fig)


    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(chans)):
        stdshape_chan = std_shape[i].data
        shape_chan = shape[i].data
        # if args.comparison is False: 
        #     chan_mask = mask[0]
        # else:
        chan_mask = mask[i]

        ax.errorbar(shape_chan,stdshape_chan, fmt="o", color=colors[i+3], label=chans[i])
        # if chan_mask.all() == False:
        #     ax.errorbar(shape_chan, stdshape_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
        # else: 
        ax.errorbar(shape_chan[~chan_mask], stdshape_chan[~chan_mask], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.set_ylabel(f"std(a*b/psf)")
    ax.set_yscale('log')
    ax.set_xlabel(f"mean(a*b/psf)")
    fig.suptitle(f"{drift}: Shape")
    ax.legend()
    plt.savefig(f"{args.project}/{drift}/{drift}_iocheck_shape.{ext}", bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure(dpi=plt.rcParams['figure.dpi']*4.0)
    ax = fig.add_subplot(1,1,1)
    for i in range(len(chans)):
        intoverpeak_chan = intoverpeak[i].data
        std_intoverpeak_chan = std_intoverpeak[i].data
        # if args.comparison is False: 
        #     chan_mask = mask[0]
        # else:
        chan_mask = mask[i]
        ax.errorbar(std_intoverpeak_chan,intoverpeak_chan, fmt="o", color=colors[i+3], label=chans[i])
        # if chan_mask.all() == False:
        #     ax.errorbar(std_intoverpeak_chan, intoverpeak_chan, fmt="o", color=colors[i+3],markeredgecolor="k")
        # else:
        ax.errorbar(std_intoverpeak_chan[~chan_mask], intoverpeak_chan[~chan_mask], fmt="o", color=colors[i+3],markeredgecolor="k")
    ax.errorbar(np.nan,np.nan, fmt="o", color='none',markeredgecolor="k", alpha=1, label="Selected")
    ax.axhline(args.intoverpeak_cut, color="k", alpha=0.3, ls="--")
    ax.axvline(args.std_cut, color="k", alpha=0.3, ls="--")
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
        '--intoverpeak_cut',
        type=float,
        default=1.1,
        help="Upper limit for mean(int/peak) for an obsid, all obsids with a value higher than this are cut. (default=1.1)"
    )
    parser.add_argument(
        '--std_cut',
        type=float,
        default=0.125,
        help="Upper limit for std(int/peak) for an obsid, all obsids with a value higher than this are cut. (default=0.125)"
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
        "--comparison",
        default=False,
        help="If not None, will look for sub, nosub, newcal images to compare"
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
        default=True,
        help="Will make a .txt file with the bad obsids (default=True) "
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
        if os.path.exists(ref_cat_file):
            hdu = fits.open(ref_cat_file)
            do_xm = hdu[1].data
            hdu.close()
    else:
        logger.warning(f"Can't find reference GGSM catalogue for xm! ")
        do_xm = False

    
    if args.comparison is False: 
        logger.warning(f"Running the comparison verison")
        obs_txtfile = [txtfile]
        extension = ["_sub", "_nosub", "_newcal", "_newmodel"]
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
            if drift in ["XG_D-27_20201022", "XG_D-27_20201008", "XG_D-27_20201001", "XG_D-40_20201014_Test"]:
                chans = ["69", "93", "121", "145", "169"]
            else: 
                chans = ["069", "093", "121", "145", "169"]
            obs_txtfile = []
            for chan in chans:
                obs_txtfile.append(txtfile.replace(".txt",f"_cenchan_{chan}.txt"))
            obs_txtfile.append(txtfile)
            if drift == "XG_D-40_20201014_Test":
                drift = "XG_D-40_20201014"
        
    # Looking for any missing obsids so they're removed before assessing 
    logger.debug(f"{obs_txtfile}")
    # Hacky just renaming this here so that everything saves with nice format but some drifts have no 0 for 69 and 93 
    chans = ["069", "093", "121", "145", "169"]
    # Reading in the obsids from txt files
    obsids = []
    missing_mask = []
    for i in range(len(chans)):
        cenchan_obsids = read_obsids(obs_txtfile[i])
        obsids_chan, missing_mask_chan = remove_missing(cenchan_obsids)
        obsids.append(obsids_chan)
        missing_mask.append(missing_mask_chan)
        logger.debug(f"Number of obsids per in channel {chans[i]}: {ma.count(obsids[i])}")
        if args.save_missing_obsids is not None: 
            logger.debug(f"Saving missing obsids")
            chan_missing_obsids = obsids_chan[obsids_chan.mask]
            chan_missing_obsids.mask = ma.nomask
            np.savetxt(obs_txtfile[i].replace(".txt", "_missing_obsids.txt"), chan_missing_obsids.compressed(), fmt="%10.0f")


    # Cutting obsids with high RMS in MFS image 
    if args.flag_high_rms is True: 
        rms_mask = cut_high_rms(obsids)
        for i in range(len(chans)):
            logger.warning(f"Number of obsids flagged for bad rms for {chans[i]}: {ma.count_masked(rms_mask[i])}")
    else: 
        logger.debug(f"Not running the high RMS cut")

    # Running cut of bad sources to assess io
    # note: first cuts bad srcs, the xm with GGSM to find nice brihgt etc ones before doing the actually assessment 

    if args.flag_bad_io is True: 
        drift_intoverpeak, drift_stdintoverpeak, drift_shape, drift_stdshape = check_io(obsids, missing_mask, do_xm)
        for i in range(len(chans)):
            color = colors[i+3]
            logger.debug(f"Running the io cut on chan {chans[i]}")
            num_obsids_precut = ma.count(obsids[i])
            logger.debug(f"STD for drift at chan {chans[i]}: {drift_intoverpeak[i].std()}")
            mask_intoverpeak = ma.masked_greater(drift_intoverpeak[i], args.intoverpeak_cut).mask
            mask_stdintoverpeak = ma.masked_greater(drift_stdintoverpeak[i],args.std_cut).mask


            obsids[i][mask_stdintoverpeak] = ma.masked
            obsids[i][mask_intoverpeak] = ma.masked

            num_obsids_postio = ma.count_masked(obsids[i])

            logger.debug(f"Obsids for chan {chans[i]} before any cuts: {ma.count(obsids[i])}")
            logger.debug(f"Number of flagged obsids for chan {chans[i]} after io cuts: {ma.count_masked(obsids[i])}")

            if args.flag_high_rms is True: 
                obsids[i][rms_mask[i].mask] = ma.masked

            frac_flagged = int(((num_obsids_postio/num_obsids_precut))*100)
            logger.debug(f"Number of obsids with flagged for bad io for {chans[i]}: {num_obsids_postio} ({frac_flagged}%)")
            if frac_flagged > 20:
                logger.warning(f"Large number of obsids flagged for bad io at chan {chans[i]}!: {frac_flagged}%")
    elif args.comparison is False: 
        logger.debug(f"Running iocheck but for comparison!")
        drift_intoverpeak = []
        drift_stdintoverpeak = []
        drift_shape = [] 
        drift_stdshape = []
        for i in range(len(extension)):
            drift_intoverpeak_ext, drift_stdintoverpeak_ext, drift_shape_ext, drift_stdshape_ext = check_io(obsids, missing_mask, do_xm, extra=extension[i])
            drift_intoverpeak.append(drift_intoverpeak_ext[0])
            drift_stdintoverpeak.append(drift_stdintoverpeak_ext[0])
            drift_shape.append(drift_shape_ext[0])
            drift_stdshape.append(drift_stdshape_ext[0])
    else: 
        drift_intoverpeak, drift_stdintoverpeak, drift_shape, drift_stdshape = check_io(obsids, missing_mask, do_xm)


    bad_io_mask = []
    all_obsids = ma.asanyarray(())
    for i in range(len(chans)):
        frac_flagged = int((ma.count_masked(obsids[i])/ma.size(obsids[i]))*100)
        logger.warning(f"OBSIDS for {chans[i]} after flagging!: {ma.count(obsids[i])}/{ma.size(obsids[i])} = {ma.count_masked(obsids[i])} ({frac_flagged}%) flagged")
        io_mask = obsids[i].mask
        bad_io_mask.append(io_mask)
        all_obsids = ma.append(all_obsids,obsids_chan)


    if args.save_bad_obsids is True: 
        logger.debug(f"Saving missing obsids")
        if len(chans) > 1:
            all_good_obsids = ma.asanyarray(())
            for i in range(len(chans)):
                obsids_chan = obsids[i]
                chan_bad_obsids = obsids_chan[obsids_chan.mask].data
                chan_good_obsids = obsids_chan[~obsids_chan.mask].data
                np.savetxt(obs_txtfile[i].replace(".txt", "_bad_obsids.txt"), chan_bad_obsids, fmt="%10.0f")
                np.savetxt(obs_txtfile[i].replace(".txt", "_bad_obsids.txt"), chan_bad_obsids, fmt="%10.0f")
                all_good_obsids = ma.append(all_good_obsids, chan_good_obsids)
            np.savetxt(obs_txtfile[-1].replace(".txt", "_bad_obsids.txt") ,all_obsids[all_obsids.mask].data, fmt="%10.0f") 
            np.savetxt(obs_txtfile[-1].replace(".txt", "_good_obsids.txt") ,all_good_obsids.data, fmt="%10.0f")   
        else: 
            np.savetxt(obs_txtfile[-1].replace(".txt", "_bad_obsids.txt") ,all_obsids[all_obsids.mask].data, fmt="%10.0f") 

    logger.warning(f"Total flagged for drift: {ma.count_masked(all_obsids)}")

    if args.plot in ["all", "min"]:   
        logger.debug(f"Plotting for drift")
        if args.comparison is False:  
            plt_io_pernight(obsids, drift_intoverpeak, drift_stdintoverpeak, drift_shape, drift_stdshape, bad_io_mask, drift, extension)
        else:
            plt_io_pernight(obsids, drift_intoverpeak, drift_stdintoverpeak, drift_shape, drift_stdshape, bad_io_mask, drift, chans)
        


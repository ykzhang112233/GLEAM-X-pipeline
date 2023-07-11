#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp2d, RectBivariateSpline, interpn, griddata
from argparse import ArgumentParser

def sigmoid(x, ref_x):
    return 1. / (1. + np.exp(-(x-ref_x)))

def create_sigweight(infits):

    img_fits = fits.open(infits)
    img_shape = img_fits[0].data.shape
    valid_mask = np.isfinite(img_fits[0].data)
    step = 50
    y, x = np.indices([s//step + 2 for s in img_fits[0].data.shape], dtype=np.int32)

    x_grid = np.clip(x*step, None, img_fits[0].data.shape[1]-1)
    y_grid = np.clip(y*step, None, img_fits[0].data.shape[0]-1)

    print('Image shape: ', img_shape)
    print('Max Indicies:', np.max(y_grid), np.max(x_grid))
    print('Min Indicies:', np.min(y_grid), np.min(x_grid))


    keep = valid_mask[y_grid.flatten(), x_grid.flatten()].reshape(y_grid.shape)
    # keep = np.ones_like(y_sparse)

    x_valid = x_grid[keep]
    y_valid = y_grid[keep]


    wcs = WCS(img_fits[0].header)

    sky_pos_grid = wcs.all_pix2world(x_grid, y_grid, 0)
    sky_pos = wcs.all_pix2world(x_valid, y_valid, 0)

    max_y = np.max(sky_pos[1])
    min_y = np.min(sky_pos[1])

    mask = sky_pos[1] > (max_y - 10)

    points = [(x_grid, y_grid, sky_pos_grid), (x_valid, y_valid, sky_pos)]

    for (xp, yp, sp) in points:

        print('Plotting sky_pos sigma')
        sig_col = np.max(
                (sigmoid(sp[1], max_y-4),
                1-sigmoid(sp[1], min_y+4)),
                axis=0
            ).astype(np.float16)
        print('Finished sigmoid')
    print('Finished')

    sp = sky_pos_grid
    xp = np.array(x_grid).astype(np.int32)
    yp = np.array(y_grid).astype(np.int32)

    sig_col = np.max(
            (sigmoid(sp[1], max_y-4),
            1-sigmoid(sp[1], min_y+4)),
            axis=0
        ).astype(np.float16)


    print(xp.shape, yp.shape, sig_col.shape)

    x_lin = np.arange(img_shape[1], dtype=np.int32)
    y_lin = np.arange(img_shape[0], dtype=np.int32)

    xx, yy = np.meshgrid(x_lin, y_lin)

    print(xx.shape, xx.dtype)

    stripes = 200
    strides = img_shape[1] // stripes

    sp = sky_pos_grid
    xp = np.array(x_grid).astype(np.int32)
    yp = np.array(y_grid).astype(np.int32)

    sig_col = np.max(
            (sigmoid(sp[1], max_y-4),
            1-sigmoid(sp[1], min_y+4)),
            axis=0
        ).astype(np.float16)

    results = []
    i=0
    while i <= (img_shape[1]-strides):
        x_lin = np.arange(img_shape[0], dtype=np.int32)
        maxy= min((i+strides),img_shape[1])
        y_lin = np.arange(
                i,maxy,
                dtype=np.int32)
        print(x_lin.max(), y_lin.max())
        yy, xx = np.meshgrid(y_lin, x_lin)


        print('interpolating')
        d = griddata(
                (xp.flatten(), yp.flatten()),
                sig_col.flatten(),
                (yy, xx),
        ).astype(np.float32)
        print(xx.shape, xx.dtype, d.dtype)

        results.append(d)
        i+=strides


    dd = np.hstack(results)

    dd[~valid_mask] = np.nan

    print(dd.shape, dd.dtype)


    ddd=1-dd

    return dd, ddd

def create_weightmap(sigmoidweight,rms):

    rms_fits = fits.open(rms)
    weightmap = sigmoidweight * (1/(rms_fits[0].data**2))

    return weightmap


if __name__ == "__main__":
    parser = ArgumentParser(
        description="create a weighting for a drift night based on a sigmoid weighting, will output the fits image for weighting that should be multiplied by the rms/weights "
    )
    parser.add_argument(
        "infits",
        help="The input image to use for pixels and sky wcs",
    )
    parser.add_argument(
        "rmsfits",
        help="RMS map to use the sigmoid weighting on and produce final weightmap thats input for swarp"
    )
    parser.add_argument(
        "outfits",
        help="Outfile to save the map to"
    )
    args = parser.parse_args()
    dd, ddd = create_sigweight(args.infits)

    # fits.writeto(args.outfits,ddd,overwrite=True)
    
    weightmap = create_weightmap(ddd, args.rmsfits)
    fits.writeto(args.outfits, weightmap, overwrite=True)

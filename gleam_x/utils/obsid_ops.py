#! /usr/bin/env python

import logging
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.stats.circstats import circmean
from astropy.coordinates import SkyCoord
from argparse import ArgumentParser
from typing import Iterable

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())


def calculate_mean(files: Iterable[str], print_pos: bool = False):

    positions = []
    for file in files:
        with fits.open(file, mode="readonly", memmap=True) as in_fits:
            ra, dec = None, None
            try:
                ra = in_fits[0].header["CRVAL1"]
                dec = in_fits[0].header["CRVAL2"]

                # print(ra, dec)
                positions.append((file, ra, dec))
            except:
                logger.warn(f"Position could not be read for {file}")

    ras = [i[1] for i in positions]
    decs = [i[2] for i in positions]

    sky_pos = SkyCoord(ras, decs, unit=(u.deg, u.deg))

    # The use of circmean is accurate only when delta(dec) is small. For the
    # GLEAM-X declination strips this is the case.
    mean_ra = circmean(sky_pos.ra)
    mean_dec = np.mean(sky_pos.dec)
    mean_pos = SkyCoord(mean_ra, mean_dec, unit=(u.deg, u.deg))

    if print_pos:
        print(f"{mean_pos.ra.deg} {mean_pos.dec.deg}")

    return mean_ra, mean_dec


if __name__ == "__main__":

    parser = ArgumentParser(
        description="Given a set of FITS files, extract their RA/Dec center and then calculate the mean position from this set"
    )

    parser.add_argument(
        "files",
        nargs="+",
        type=str,
        help="Path to the (many) FITS files to extract the position from ",
    )

    args = parser.parse_args()

    calculate_mean(args.files, print_pos=True)

#!/usr/bin/env python

__author__ = ("Kat Ross")
__date__ = "13/03/2023"

import csv
from argparse import ArgumentParser


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--mosaicnm",
        type=str,
        dest="mosaicnm",
        help="Name of the mosaic for the catalogue string"

    )
    args = parser.parse_args()

    mosaicnm = args.mosaicnm

    extensions = [
        "072-080",
        "072-103",
        "080-088",
        "088-095",
        "095-103",
        "103-111",
        "103-134",
        "111-118",
        "118-126",
        "126-134",
        "139-147",
        "139-170",
        "147-154",
        "154-162",
        "162-170",
        "170-177",
        "170-200",
        "177-185",
        "185-193",
        "193-200",
        "200-208",
        "200-231",
        "208-216",
        "216-223",
        "223-231"
    ]

    freq_suffixes = [
        "072_080",
        "072_103",
        "080_088",
        "088_095",
        "095_103",
        "103_111",
        "103_134",
        "111_118",
        "118_126",
        "126_134",
        "139_147",
        "139_170",
        "147_154",
        "154_162",
        "162_170",
        "170_177",
        "170_200",
        "177_185",
        "185_193",
        "193_200",
        "200_208",
        "200_231",
        "208_216",
        "216_223",
        "223_231"
    ]

    catalogues = []
    names = []
    images = []
    bkg = []
    rms = []
    psf = []
    suffixes = []
    prefixes = []
    rescaled_cats = []


    for i in range(len(extensions)):
        ext=extensions[i]
        catalogues.append(f"{mosaicnm}_{ext}MHz_ddmod_prior_comp.fits")
        rescaled_cats.append(f"{mosaicnm}_{ext}MHz_ddmod_prior_comp_rescaled.fits")
        # if i in [2, 7, 12, 17, 22]:
        #     suffixes.append(f"_W_{freq_suffixes[i]}MHz")
        # else: 
        suffixes.append(f"_{freq_suffixes[i]}MHz")
        prefixes.append("")
        names.append(f"{mosaicnm}_{ext}MHz")
        images.append(f"{mosaicnm}_{ext}MHz_ddmod.fits")
        bkg.append(f"{mosaicnm}_{ext}MHz_ddmod_bkg.fits")
        rms.append(f"{mosaicnm}_{ext}MHz_ddmod_rms.fits")
        psf.append(f"{mosaicnm}_{ext}MHz_projpsf_psf.fits")

    csv_catalogues_fields = ["#catalogue", "prefix", "suffix"]
    csv_catalogues = [catalogues, prefixes, suffixes]
    csv_cats = zip(*csv_catalogues)
    csv_catalogues_rescaled = [rescaled_cats, prefixes, suffixes]
    csv_cats_rescaled = zip(*csv_catalogues_rescaled)

    csv_image_fields = ["name", "image", "bkg", "rms", "psf"]
    csv_image = [names, images, bkg, rms, psf]    
    csv_ims = zip(*csv_image)

    with open(f"{mosaicnm}_catalogues.csv", "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csv_catalogues_fields)
        csvwriter.writerows(csv_cats)
    with open(f"{mosaicnm}_images.csv", "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csv_image_fields)
        csvwriter.writerows(csv_ims)
    with open(f"{mosaicnm}_catalogues_rescaled.csv", "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(csv_catalogues_fields)
        csvwriter.writerows(csv_cats_rescaled)
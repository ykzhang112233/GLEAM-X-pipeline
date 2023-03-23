#! /usr/bin/env python3

import astropy
from astropy.table import Table, Column
from astropy.io.votable import from_table, writeto
import numpy as np
import glob
import argparse
import sys

__author__ = "Paul Hancock"
__date__ = "2021-04-12"

def get_epoch_catalogues(epochs_file):
    """
    Read a file which contains a list of the catalogues to be read
    in the format of "filename,prefix,suffix"
    prefix/suffix may be empty/blank

    parameters
    ----------
    epochs_files : str
        A file which has a list of catalogues, one per line.

    returns
    -------
    files : list
        A list of filenames
    prefix : list
        A list of prefixes
    suffix : list
        A list of suffixes
    """
    lines = list(map(str.strip, open(epochs_file).readlines()))
    words = [l.split(',')[:3] for l in lines if not l.startswith('#')]
    files, prefix, suffix = zip(*words)
    return files, prefix, suffix


def join_catalogues(reference, epochs, prefix=None, suffix=None, all_cols=False):
    """
    Take a reference cataloge, strip all but the uuid/ra/dec columns and then
    join the flux/err_flux data from each of the epoch catalogues
    From each epoch we extract the peak_flux and err_peak_flux columns and rename
    them by appending _N where N is the epoch number

    parameters
    ----------
    reference : str
        Filename for the reference catalogue

    epochs : list
        A list of the individual epoch file names

    prefix : list
        Column name prefix for each epoch

    suffix : list
        Column name suffix for each epoch

    all_cols : bool
        If true then all columns are retained from all catalogues

    returns
    -------
    table : `astropy.table.Table`
        A joined table
    """
    # Read the reference catalogue and retain only the uuid and ra/dec columns
    # rename the ra/dec columns
    print("Using reference catalogue {0}".format(reference))
    if all_cols:
        ref = Table.read(reference)
        keep_cols = ref.colnames
    else:
        ref = Table.read(reference)['uuid', 'ra','dec']
        keep_cols = ['peak_flux','err_peak_flux','local_rms','background',
                     'int_flux', 'err_int_flux', 'a','b','pa','psf_a','psf_b','psf_pa',
                     'residual_mean','residual_std']
    dtypes = Table.read(reference).dtype
    ref.rename_column('ra', 'ref_ra')
    ref.rename_column('dec', 'ref_dec')
    ref.sort(keys='uuid')

    # create lists for when prefix/suffix are None
    if prefix is None:
        prefix = ['']*len(ref)
    if suffix is None:
        suffix = ['']*len(ref)

    # make the empty columns
    new_cols =[]
    #data = np.zeros(len(ref), dtype=np.float32)

    for i,(p,s) in enumerate(zip(prefix,suffix)):
        for a in keep_cols:
            colname = '{0}'+a+'{1}'
            new_cols.append(Column(length=len(ref), name=colname.format(p,s), dtype=dtypes[a]))
    print("ref table is {0} rows".format(len(ref)))
    ref.add_columns(new_cols)
    del new_cols

    read_cols = keep_cols[:]
    if 'uuid' not in read_cols:
        read_cols.append('uuid')
    # now put the data into the new big table
    for i,(f,p,s) in enumerate(zip(epochs,prefix,suffix)):
        print("Joining epoch {0} catalogue {1}".format(i,f))
        new_cols = Table.read(f)[read_cols]
        new_cols.sort(keys='uuid')
        # compute the order/presence
        ordering = np.argwhere(np.in1d(ref['uuid'], new_cols['uuid'], assume_unique=True))[:,0]
        for a in keep_cols :
            ref['{0}{1}{2}'.format(p,a,s)][ordering] = new_cols[a]

    return ref


def write_table(tab, filename):
    """
    Write a VOTable using a binary format.

    parameters
    ----------
    tab : `astropy.table.Table`
        The table to be written

    filename : str
        The output filename. Should end with .xml or .vot
    """
    vot = from_table(tab)
    vot.set_all_tables_format('binary')
    vot.to_xml(filename)
    print("Wrote {0}".format(filename))
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group1 = parser.add_argument_group("Combine catalogues into a flux table")
    group1.add_argument("--refcat", dest='ref', type=str, default=None,
                        help="The reference catalogue")
    group1.add_argument("--epochs", dest='epochs', type=str, default=None,
                        help="A file containing 'filename,prefix,suffix' for each catalogue.")
    group1.add_argument("--out", dest='outfile', type=str, default='flux_table.vot',
                        help="The output table name. Default = flux_table.vot")
    group1.add_argument("--all", dest='all', action='store_true', default=False,
                        help="Retain all columns from all catalogues. [Default false]")
    results = parser.parse_args()

    if None in (results.ref, results.epochs, results.outfile):
        parser.print_help()
        sys.exit()

    files, prefix, suffix = get_epoch_catalogues(results.epochs)
    table = join_catalogues(results.ref, files, prefix, suffix, all_cols=results.all)
    write_table(table, results.outfile)

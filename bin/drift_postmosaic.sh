#! /bin/bash

usage()
{
echo "drift_postmosaic.sh [-p project] [-d dep] [-q queue] [-a account] [-t] -o list_of_observations.txt
  -p project  : project, (must be specified, no default)
  -d dep      : job number for dependency (afterok)
  -t          : test. Don't submit job, just make the batch file
                and then return the submission command  
  -m mosaicdir: Directory name where mosaics stored (default=project/mosaic) 
  -o obslist  : the list of obsids to process" 1>&2;
exit 1;
}


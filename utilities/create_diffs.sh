#! /bin/bash -l

# A simple script to help work out the correct calibration files to apply
# Idea i to take the original calibration applied, and then decide on new
# calibrations with the segements method. 

# Actually, thinking about it this might not be the best way, as the original
# problem was in John's code. Will need to re-reconsider this again.

set -e

if [[ -z $GXBASE ]]
then
    echo "Looks like the GLEAM-X profile is not loaded."
    exit 1
fi

if [[ -z $1 ]]
then
    echo "No file passed in. Expected usuage"
    echo "create_diffs.sh OBSIDFILE"
    exit 3
fi

obs=$1

if [[ ! -f "$obs" ]]
then
   echo "It looks like $obs does not exist."
   exit 2
fi

echo "Checking binary files to apply calibration solutions"
echo "Original checks"
singularity exec $GXCONTAINER check_assign_solutions.py assign $obs | grep -v A  > original_calids.txt
wc -l original_calids.txt

# The incorrect calibration solutions ultimately slipped in after the XY/YX terms were set to zero. 
# The code to count the number of flagged channels/antenna counts the NANs in a solution file. Each
# solution is really a 2x2 Jones matrix of complex numbers. Code would count all NANs with the implicit
# assumption that a Jones is either entirily finite or completely masked. Setting the XY/YX to zero
# broke this. So, to figure out which obsids would need to be reprocessed, rather than re-calibrating
# with the now corrected antenna-referencing and XY/YX nuking script, we can use the previous one
echo "Four subvband checks against non-initial ref solutions (before XY zero terms introduces"
singularity exec $GXCONTAINER check_assign_solutions.py -s 4 assign $obs --suffix "_local_gleam_model_solutions_initial.bin" | grep -v A  > recent_calids.txt
wc -l recent_calids.txt

echo "Creating the diff"
diff original_calids.txt recent_calids.txt | tee diff_output.txt

echo "Extracting obsid-calid pairs to reprocess"
cat diff_output.txt | grep '>' | cut -d ' ' -f2- > new_obsid_calid.txt

echo "Items to fix"
wc -l new_obsid_calid.txt


data=new_cals
if [[ -d "$data" ]]
then
    echo "Directory $data exists. Exiting. "
    exit 2
fi

mkdir $data
cat new_obsid_calid.txt | cut -d ' ' -f2 > "${data}/Repair_${obs}"
mv new_obsid_calid.txt $data




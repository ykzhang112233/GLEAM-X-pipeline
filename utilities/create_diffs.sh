#! /bin/bash -l

# A simple script to help work out the correct calibration files to apply
# Idea i to take the original calibration applied, and then decide on new
# calibrations with the segements method. 

# Actually, thinking about it this might not be the best way, as the original
# problem was in John's code. Will need to re-reconsider this again.

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

echo "Checking binaru files to apply calibration solutions"
echo "Original checks"
singularity exec $GXCONTAINER check_assign_solutions.py assign $obs | grep -v A  > original_calids.txt
wc -l original_calids.txt

echo "Four subvband checks"
singularity exec $GXCONTAINER check_assign_solutions.py -s 4 assign $obs | grep -v A  > recent_calids.txt
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
cat new_obsid_calid.txt | cut -d ' ' -f3 > "${data}/Repair_${obs}"
mv new_obsid_calid.txt $data




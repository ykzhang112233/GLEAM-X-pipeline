# A deadset simple one-time use script meant to parse an apply_cal output
# log to extract the obsid and the corresponding calid applied. These
# are then inserted back into the GLEAM-X database, (hopefully) restoring
# it to a state before the great Nimbus database deletion. 
file=$1

echo $file
obsid=$(cat $file | grep 'apply_cal obsid' | cut -d ' ' -f3)
calid=$(cat $file | grep 'apply_cal calid' | cut -d ' ' -f3)

echo $obsid $calid

track_task.py obs_calibrator --obs_id $obsid --cal_id $calid

# A deadset simple, one-time use script to update values in the GLEAM-X
# database, restoring them to values that were present before the 
# nimbus instance was lost. This is meant to operate against a filename
# of the form: obsid_iono.csv
# Since these csv files were pulled of the Data-Central archive, we 
# can also update their observational status to archived as well. 
file=$1
obsid=$(echo $file | cut -d '_' -f1)

iono_update.py --ionocsv "$file"
track_task.py obs_status --obs_id="${obsid}" --status='archived'

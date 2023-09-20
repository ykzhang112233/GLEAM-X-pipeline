#!/bin/bash -l 
#SBATCH --partition=copy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --job-name=prep_nightcoadd_XG_D+03
#SBATCH --output=/astro/mwasci/kross/GLEAM-X-pipeline/log_garrawarla/acacia_send_XG_D+03_4Night.o%A
#SBATCH --error=/astro/mwasci/kross/GLEAM-X-pipeline/log_garrawarla/acacia_send_XG_D+03_4Night.e%A
#SBATCH --account=pawsey0272
#SBATCH --open-mode=append
#SBATCH --parsable
trap 'echo "Requeuing job."; scontrol requeue $SLURM_JOB_ID;' INT
module load rclone/1.59.2

cd /scratch/pawsey0272/kross/gleamx/XG_D+03_4Night/

srun rclone sync -P --transfers 4 --checkers 4 acacia:gleamxkross/XG_D+03_20201003.tar.gz/XG_D+03_20201003.tar.gz ./  &

srun rclone sync -P --transfers 4 --checkers 4 acacia:gleamxkross/XG_D+03_20201010.tar.gz/XG_D+03_20201010.tar.gz ./  &

srun rclone sync -P --transfers 4 --checkers 4 acacia:gleamxkross/XG_D+03_20201017.tar.gz/XG_D+03_20201017.tar.gz ./  &

srun rclone sync -P --transfers 4 --checkers 4 acacia:gleamxkross/XG_D+03_20201024.tar.gz/XG_D+03_20201024.tar.gz ./  &

# srun tar -xvzf XG_D+20_20201025.tar.gz ./mosaic/



wait $!

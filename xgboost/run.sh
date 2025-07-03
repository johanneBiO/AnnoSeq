#! /bin/bash
#SBATCH --job-name=xgboost
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=240:00:00
#SBATCH --nodelist=node06
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate xgboost

INPUT_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_processed/attention/attn_headQuan90_colQuan.pkl"
RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/results/xgboost_headQuan90_colQuan_01000_10iter/"

time python train_xgboost_models.py -f $INPUT_FILE -res $RES_DIR -n 10
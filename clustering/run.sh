#! /bin/bash
#SBATCH --job-name=getemb
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --nodelist=node06
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/clustering/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/clustering/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate emb

# File paths
ATT_DATA="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_processed/attention/attn_headQuan90_colQuan.pkl"
MOD_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/xgboost/results/xgboost_headQuan90_colQuan_01000_1094iter/models/"
EMB_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_output/esm_headQuan90_colQuan/embedding/embedding.json"
ACC_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/additional/accessions.txt"
RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_processed/region_embedding/"
ANNO_RDS="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/complete/annotations/processed/anno.rds"

python main.py -att $ATT_DATA -mod $MOD_DIR -emb $EMB_FILE -acc $ACC_FILE -res $RES_DIR #-a -a_rds $ANNO_RDS
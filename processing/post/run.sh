#! /bin/bash
#SBATCH --job-name=post
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --nodelist=node05
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/processing/post/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/processing/post/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate class

# File paths
ATTN_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_output_quan/attention/raw/"
ACC_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/additional/accessions.txt"
ANNO_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/complete/annotations/processed/anno_expanded.rds"
PART_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/additional/partition.csv"
RES_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_processed/attn_raw_q2.pkl"
RES_CSV_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/esm_processed/attn_raw_q2.csv"

python process_esm_attention.py -at $ATTN_FILE -ac $ACC_FILE -an $ANNO_FILE -p $PART_FILE -r $RES_FILE -csv -r_csv $RES_CSV_FILE -q 
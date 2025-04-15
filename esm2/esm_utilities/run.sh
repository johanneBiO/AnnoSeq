#! /bin/bash
#SBATCH --job-name=gpujob_esm2
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=shard:10
#SBATCH --mem=20G
#SBATCH --time=5:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate esm

INPUT_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/sequences/biological_seq/seq_00100.fasta"
RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/"
SCR_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/sequences/scrambled_seq/"

python main.py -seq $INPUT_FILE -res $RES_DIR -s -scr $SCR_DIR -n -bs 10
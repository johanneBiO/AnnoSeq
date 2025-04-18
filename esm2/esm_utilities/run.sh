#! /bin/bash
#SBATCH --job-name=gpujob_esm2
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=shard:23
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate esm

### Subset 100
#INPUT_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/sequences/biological_seq/seq_00100.fasta"
#RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/"
#SCR_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_00100/sequences/scrambled_seq/"

#python main.py -seq $INPUT_FILE -res $RES_DIR -s -scr $SCR_DIR -n

### Subset 1000
#INPUT_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/sequences/biological_seq/seq_01000.fasta"
#RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/"
#SCR_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/subset_01000/sequences/scrambled_seq/"

#python main.py -seq $INPUT_FILE -res $RES_DIR -s -scr $SCR_DIR -bs 10

### Complete
INPUT_FILE="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/complete/sequences/seq_complete_sp_cropped.fasta"
RES_DIR="/net/mimer/mnt/tank/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/data/complete/"

python main.py -seq $INPUT_FILE -res $RES_DIR -bs 10
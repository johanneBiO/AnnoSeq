#! /bin/bash
#SBATCH --job-name=gpujob_esm2
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=shard:10
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --nodelist=compute02,compute03,compute04,compute05
#SBATCH --output=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.out 
#SBATCH --error=/home/projects2/kvs_students/2025/jbo_unbiased_seq_annot/master_thesis/esm2/esm_utilities/slurm_logs/job-%j.err 

source /home/people/jobao/miniconda3/etc/profile.d/conda.sh
conda activate esm

INPUT_FILE="example_proteins/example_proteins.fasta"
RES_DIR="example_proteins/"

python esm_main.py -seq $INPUT_FILE -res $RES_DIR -n
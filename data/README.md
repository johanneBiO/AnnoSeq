# Data Folder

The data files are not included in this repository due to the size. 

## Source

This project utilizes the human reference proteome from UniProt (Proteome ID: UP000005640). The sequences were obtained from the [UniProt FTP server](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/) in FASTA format, and corresponding sequence annotation data was downloaded in XML. TrEMBL entries, which have not been reviewed, were excluded.

## Data Organization

XXX Needs refinement XXX

    .
    ├── _raw                                    # Directory to store the raw data 
    │   ├── UP000005640_9606.xml                # Sequence Annotations (features)
    │   └── UP000005640_9606.fasta              # Sequences
    │
    ├── complete                                # Complete dataset (processed) ~ 20.000 sequences
    │   ├── additional                   
    │   │   ├── seq_lengths_adj.csv             # Table for cropping sequences
    │   │   └── ...
    │   │
    │   ├── annotations        
    │   │   ├── features.rds
    │   │   ├── features_final.rds
    │   │   └── ...           
    │   │   
    │   ├── sequences
    │   │   ├── UP000005640_9606_sp.fasta
    │   │   └── UP000005640_9606_sp_cropped.fasta             
    │   │
    │   └── esm_output
    │       ├── attention
    │       │   └── raw
    │       │       ├── layer_01_attn_raw.json
    │       │       ├── layer_02_attn_raw.json
    │       │       ├── layer_03_attn_raw.json
    │       │       └── ...
    │       │
    │       └── embedding
    │           └── embedding.json
    │       
    ├── subset_00100 
    │                 
    └── subset_01000
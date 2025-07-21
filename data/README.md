# Data Folder

The data files are not included in this repository due to the size. 

## Source

This project utilizes the human reference proteome from UniProt (Proteome ID: UP000005640). The sequences were obtained from the [UniProt FTP server](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/) in FASTA format, and corresponding sequence annotation data was downloaded in XML. TrEMBL entries, which have not been reviewed, were excluded.

In addition, manually curated intrinsically disordered regions are downloaded from https://disprot.org/.

## Data Organization

The overall folder structure can be seen below. All three datasets follow the same organization as seen for the subset of 100 sequences. 

    .
    ├── _raw                                    # Directory to store the raw data 
    │   ├── UP000005640_9606.xml                # Sequence Annotations (features)
    │   └── UP000005640_9606.fasta              # Sequences
    │
    ├── complete                                # Complete dataset (20.000 sequences)
    │       
    ├── subset_01000                            # Subset of 1000 sequences
    │                 
    ├── subset_00100                            # Subset of 100 sequences
    │   ├── sequences                           # Selected sequences
    │   │
    │   ├── esm_output                          # Per-residue embeddings and raw attention signals after running ESM-2
    │   │   
    │   ├── esm_processed                       # Region embeddings and processed attention signals         
    │   │
    │   └── additional                          # Accessions and data partitioning
    │
    └── external                                # DisProt data 
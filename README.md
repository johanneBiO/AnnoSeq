# Universal Annotation of Biological Sequences using Large Language Models

Master Thesis by Johanne Badsberg Overgaard

## Project Overview

This repository contains the code and analyses for my master's thesis project focused on developing a new universal approach for protein sequence annotation. The goal was to explore how raw attention weights from the ESM-2 language model can be used to extract biologically relevant regions annotated by clustering the regions based on their ESM-2 embeddings. This approach aims to enhance our understanding of protein function through refinement of already annotated sequences and annotation of novel sequences. 

The project has been supervised by Associate Professor Kristoffer Vitting-Seerup, Group leader of The Isoform Analysis Group at the Department of Health Technology, Technical University of Denmark (DTU).

## Project Organization

.
├── /data/                    # Directory to store data (data not included in repo)
│   └── README.md             # Instructions on how to obtain the data and suggestion for data folder structure
│
├── /esm2/                    
│   ├── esm                   # ESM model code with modifications
│   │   ├── ...
│   │
│   └── esm_utilities         
│       ├── utils
│       ├── main.py
│       ├── run.sh
│       └── esm_evn.yml              
│
├── /processing/              # Scripts for ...
│   ├──              
│   ├── 
│   └── ...                  
│
├── /analysis/                # Scripts for ...
│   ├──              
│   ├── 
│   └── ...           
│
├── LICENSE                   # ???
└── README.md               

## Acknowledgement

This repository contains a modified version of the ESM-2 code, developed by the Meta Fundamental AI Research Protein Team (FAIR) and available at [https://github.com/facebookresearch/esm](https://github.com/facebookresearch/esm). The ESM-2 code is located in the `esm2/esm` folder and is licensed under the MIT License. See the `esm2/esm/LICENSE` file for details. The modifications are limited to the extraction and saving of raw and normalized attention scores, which are automatically averaged over attention heads by default. These changes cannot be disabled at this moment; the attention weights will always be saved in this manner. All other aspects of the original code remain unchanged.

The implementation of the model was inspired from the repositories [esm2_uilities from mnielLab](https://github.com/mnielLab/esm2_utilities) and [esm-utils from chihs-dtu](https://github.com/chihs-dtu/esm-utils/blob/main/utils/aggregate_results.py). Additionally, code shared by Isabella Østerlund from the Department of Health Technology, DTU contributed to optimization of the code.
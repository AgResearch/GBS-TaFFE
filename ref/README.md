# Directory for reference databases.
```ref/``` is currently kept in ```.gitignore``` due to the large files sizes. #TODOs include either: 
1. centralising these resources in a standard location and CI their versioning.
2. including rules to download and build these as needed in the smk workflow

## Current ```ref/``` directory structure 
as of: 20-10-2022
```
ref/
├── biobakery
├── GTDB
│   └── data.ace.uq.edu.au
│       ├── icons
│       └── public
│           └── gtdb
│               └── data
│                   └── releases
│                       ├── latest
│                       │   ├── auxillary_files
│                       │   ├── genomic_files_all
│                       │   └── genomic_files_reps
│                       ├── release202
│                       ├── release207
│                       ├── release80
│                       ├── release83
│                       ├── release86
│                       ├── release89
│                       └── release95
├── kraken
├── SILVA
│   └── chocophlan
└── zymoMCS
    └── ZymoBIOMICS.STD.refseq.v2
        ├── Genomes
        └── ssrRNAs
```
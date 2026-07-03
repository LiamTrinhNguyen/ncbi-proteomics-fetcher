## NCBI Proteomics Fetcher Pipeline
**An automated dry-lab utility script that interfaces with NCBI databases via Entrez APIs to find, download, and structure protein records for subsequent mass spectrometry analysis.**

## Input Modes
    1. **NCBI Query (Default):**
        Automatically queries the NCBI Protein database using Entrez E-utilities based on user-defined search terms (e.g., organism and protein name).
    2. **Interactive Refinement:**
        Supports custom search queries and result limits to ensure precise retrieval of specific proteoforms.

## Biological & Analytical Rationale
    This pipeline follows rigorous bioinformatics best practices for proteomic data preparation:
        - Dynamic Data Retrieval: Interfaces directly with NCBI Entrez APIs to maintain synchronization with the latest sequence annotations.
        - Statistical Precision: Employs Biopython to ensure FASTA record integrity and Polars for high-performance sequence structuring.
        - Enzymatic Simulation: Provides an in-silico digestion engine (Trypsin, Chymotrypsin, Pepsin, CNBr) to simulate wet-lab peptide generation.
        - Data Science Integration: Includes cleaning, log-transformation, and median normalization steps for experimental mass spectrometry intensity matrices.
        - Functional Mapping: Automatically links significantly altered proteins to functional interaction networks via STRING-DB.


## Prerequisites
    - Environment: Linux (recommended: GitHub Codespaces or Ubuntu)
    - Package Manager: Conda 
    - Core Dependencies: biopython, polars, tqdm
    - Storage: Optimized for minimal storage through isolated, timestamped run directories.


## Pipeline Steps Overview
| Step | Tool                  | Purpose                                      | Biological Reasoning                                      |
|------|-----------------------|----------------------------------------------|-----------------------------------------------------------|
| 1    | Workspace Setup       | Create isolated run directories              | Maintains reproducibility and prevents cross-run data mixing |
| 2    | Entrez.esearch        | Query NCBI Protein database                  | Retrieves high-confidence protein targets for analysis    |
| 3    | Entrez.efetch         | Download FASTA sequences                     | Streams full protein records directly from authoritative source |
| 4    | BioPython + Polars    | Parse FASTA into structured DataFrame        | Converts raw sequences into analysis-ready tabular format |
| 5    | Custom Cleavage Engine| In-silico enzymatic digestion                | Simulates trypsin (or other enzymes) to generate MS-ready peptides |
| 6    | Data Cleaning         | Filter contaminants & low-confidence hits    | Removes artifacts (keratin, decoys) and ensures high FDR confidence |
| 7    | Polars Analytics      | Normalization, Differential Expression & Pathway Mapping | Identifies biomarkers and maps biological networks |

## Getting Started
1. **Clone the repository**
```Bash
git clone https://github.com/yourusername/ncbi-proteomics-fetcher.git
cd ncbi-proteomics-fetcher
```

2. **Create and activate the environment**
```Bash
conda env create -f proteomics_analysis/environment.yaml
conda activate proteomics-fetcher
```
3. **Run the pipeline**
```Bash
python proteomics_analysis/fetch_proteome.py
```

## Outputs
    Each run creates a timestamped folder under result/run_YYYYMMDD_HHMMSS/:
        - proteomics_data.csv – Structured sequence metadata matrix.
        - digested_peptides.csv – In-silico generated peptide fragment library.
        - cleaned_proteomics.tsv – High-confidence filtered expression data.
        - differential_markers.tsv – Statistically isolated disease biomarker candidates.
        - pipeline_*.log – Full audit trail including statistical summaries of the run.

## Features
    1. Dynamic Cleavage Registry: Add or modify enzymatic cleavage rules (e.g., Trypsin vs Chymotrypsin) without altering core logic.
    2. Data Decomposition: Uses Polars to handle massive proteomic sequences with multi-threaded efficiency.
    3. Production-Ready: Includes logging, runtime error handling, and absolute path resolution to guarantee stability in cloud environments.
    4. Interactive Dashboard: Provides a single-step entry point followed by automated pipeline execution.

## Tools Used
    - Database Interface: Biopython (Entrez APIs)
    - Data Processing: Polars + tqdm
    - Analytical Workflow: Custom Enzymatic Registry Engine
    - Functional Mapping: STRING-DB API integration


## Author
Liam TrinhNguyen
Data Scientist, Wisconsin State Laboratory of Hygiene (WSLH)


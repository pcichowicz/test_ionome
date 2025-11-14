# < Ionome >

**Ionome** is an LC–MS data analysis pipeline designed for exploratory metabolomics, peak extraction, and automated chromatographic processing.  
The long-term goal of this project is to provide a modular, reproducible workflow for working with mzML files, extracting XICs, computing TIC/BPC, managing datasets, and generating a final report.

This repository currently focuses on:
- A clean, versioned Python codebase (starting with **Ver_1.0.0**)  
- Core readers/controllers for mzML processing  
- Chromatogram generation (TIC, BPC, XIC)
- Jupyter-based exploratory analysis

More components (peak deconvolution, automated annotation, QC workflows, dataset comparison utilities, and a full CLI) will be added as development continues.

---

## Installation

## Installation

**Status:** _Currently in development — not yet distributed as a package._

Current dependencies for `Ionome` are as follows:

pandas~=2.3.3  
pymzml~=2.5.11  
matplotlib~=3.10.7  
PyYAML~=6.0.3  
numpy~=2.3.4  
tqdm~=4.67.1  

## Repository Structure
ionome/  
│  
├── src/  
│   ├── correct_baseline.py              # MzMLReader class  
│   ├── extractors.py          # XICExtractor and helpers  
│   ├── controller.py          # SampleAnalysisController  
│   ├── extractors.py          # XICExtractor and helpers  
│   ├── controller.py          # SampleAnalysisController  
│   ├── extractors.py          # XICExtractor and helpers  
│   ├── controller.py          # SampleAnalysisController  
│   └── utils/                 # metadata, math helpers, io helpers  
│  
├── notebooks/  
│   └── 01_explore_single_file.ipynb  
│  
├── data/  
│   ├── raw/                   # raw mzML files (gitignored)  
│   └── processed/             # cleaned or exported tables  
│  
├── scripts/  
│   └── download_data.py       # pulls MassIVE files programmatically  
│  
├── results/                   # figures, exports (gitignored)  
├── config/                    # metadata dictionaries, settings  
├── docs/                      # expanded documentation  
├── requirements.txt  
├── README.md  
├── LICENSE  
└── .gitignore  


## Development Status
Development Status

Ver_1:

mzML loading
TIC/BPC computation
Basic XIC extraction
Preliminary metadata dictionaries
Jupyter notebook exploration

Next goals (Ver_2):
Peak detection and deconvolution
Compound annotation tables
Sample comparison (e.g., SL005 vs SL011)
Batch pipelines for entire datasets
Cleaner controller abstractions
Auto-QC metrics (noise, intensity drift, retention time shifts)

## LICENSE
License

MIT License — see LICENSE for details.

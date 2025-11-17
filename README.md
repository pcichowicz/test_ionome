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

## Installation
**Status:** _Currently in development — not yet distributed as a package._

Current dependencies for `Ionome` are as follows:

pandas~=2.3.3  
pymzml~=2.5.11  
matplotlib~=3.10.7  
PyYAML~=6.0.3  
numpy~=2.3.4  
tqdm~=4.67.1  
pyarrow
fastparquet

## Usage
Quick step-by-step for how the package works.    
### Repository Structure
```
ionome/  
│  
├── src/  
│   ├── correct_baseline.py     # BaselineCorrection class  
│   ├── helpers.py              # helper functions  
│   ├── ionome_core.py          # class Controller  
│   ├── main.py                 # XICExtractor and helpers  
│   ├── paths.py                # path module  
│   ├── peakDetection.py        # PeakDetection class  
│   ├── preprocess.py           # MzmkParser class  
│   ├── sampleMetaData.py       # SampleMetaData Controller
│   └── visualization.py        # Chromatogram class  
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
```

### Loading Chromatograms
In order for this to work, you should have two files located in the `config` directory:  
- config.yaml
- samples.yaml .

```yaml
# config.yaml

# General pipeline settings
rerun: false
save_parquet: true
verbose: true

# Paths
data_dir: "data/raw"
cached_dir: "data/processed"
results_dir: "results"
logs_dir: "logs"

# Metabolites for analysis
target_mz_list:
  p-coumerate: 163.04007
  urocanate: 137.03565
  caffeate: 179.03498
  cinnamate: 147.04460
  ferulate: 193.05009
  sinapate: 223.06065

...

# Peak detection parameters
peak_detection:
  window: 5
  precision: 9
  prominence: 0.01
  rel_height: 1
  buffer: 0

# Plotting options
plotting:
  tic: true
  bpc: true
  xic: true
  save_figures: true

```

```yaml
samples:
  - unique_id: SL011_EA_1
    id: SL011
    file: 031__20210415__SL011__EA_Blank__nosplit.mzML
    condition: Control
    description: Extraction Blank
    replicate: 1
    species: H. filiformis
  
  - unique_id: SL011_PooledQC_1
    id: SL011
    file: 033__20210415__SL011__PooledQC__nosplit.mzML
    condition: Control
    description: Pooled quality control
    replicate: 1
    species: H. filiformis
  
  - unique_id: SL011_HS_1
    id: SL011
    file: 045__20210415__SL011__HS1__nosplit.mzML
    condition: Treatment
    description: Sample
    replicate: 1
    species: H. filiformis

...
```

You can then run `Ionome` class object to initialize the workflow.
```python
from src.ionome_core import Ionome
example = Ionome(run_id="SL005", samples="samples_SL011.yaml")
```


## LICENSE
License

MIT License — see LICENSE for details.

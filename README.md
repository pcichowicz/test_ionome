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

### Running Ionome package
In order for this to work, you should have two files located in the `config` directory:  
- <details>
    <summary>see config example</summary>

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
    
    # tolerance for mz detection
    target_mz_params:
      tol: 0.01
    
    # Parsing and preprocessing of mzML files
    parser:
      ms_level: 1
    
    # Baseline correction parameters
    baseline:
      method: asls
      asls:
        lambda: 1e5
        p: 0.01
        niter: 10
        tol: 1e-6
      snip:
        window: 5
        precision: 9
        clip_negatives: true
    
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
  </details>

- <details>
    <summary>see sample list example</summary>
    
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
  
    ```
  </details>  

You can then run `Ionome` class object to initialize the workflow.

```python
from src.ionome_core import Ionome
example = Ionome(run_id="SL005", samples="samples_SL011.yaml")
```

    >[2025-11-17 06:11:16][Ionome.__init__]
    	Initializing run ID SL005 for analysis

Then you can load the spectra data from mzml files using `load.data()`:

```python
example.load_data()
```

    >[2025-11-17 06:11:44][Ionome.load_data]
    	> Loading spectra data for sample SL005_EA_1...
    	 Parsing 010__20201112__SL005__EA_Blank__nosplit.mzML ... 
    	 Caching spectra data
    	> Loading spectra data for sample SL005_EA_2...
    	 Parsing 011__20201112__SL005__EA_Blank__nosplit.mzML ... 
    	 Caching spectra data
             ...

```python
first_raw_df = example.data['SL005_EL_24hr']
first_raw_df
```

<div>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ms_level</th>
      <th>scan_id</th>
      <th>retention_time</th>
      <th>intensity</th>
      <th>mz</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>1</td>
      <td>4.293533</td>
      <td>24314.119141</td>
      <td>51.040001</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>1</td>
      <td>4.293533</td>
      <td>53341.289062</td>
      <td>52.060001</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>1</td>
      <td>4.293533</td>
      <td>6819.981934</td>
      <td>53.049999</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>1</td>
      <td>4.293533</td>
      <td>1238.951294</td>
      <td>53.980000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>1</td>
      <td>4.293533</td>
      <td>4658.718750</td>
      <td>55.020000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1792196</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990900</td>
      <td>1.795156</td>
      <td>596.169983</td>
    </tr>
    <tr>
      <th>1792197</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990900</td>
      <td>1.711735</td>
      <td>597.260010</td>
    </tr>
    <tr>
      <th>1792198</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990900</td>
      <td>3.206854</td>
      <td>597.969971</td>
    </tr>
    <tr>
      <th>1792199</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990900</td>
      <td>1.128073</td>
      <td>598.570007</td>
    </tr>
    <tr>
      <th>1792200</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990900</td>
      <td>3.180370</td>
      <td>599.280029</td>
    </tr>
  </tbody>
</table>
<p>1792201 rows × 5 columns</p>
</div>

If you have listed the metabolites (with mz target) you want to analyze in the config file, run `extract_xic()` to filter the mz

## LICENSE
License

MIT License — see LICENSE for details.

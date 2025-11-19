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

**Ver 1.0.0:**

mzML loading  
TIC/BPC computation  
Basic XIC extraction  
Quality control plots (Total Ion Chormatograms & Base Peak Chromatograms)  
Preliminary metadata dictionaries  
Jupyter notebook exploration  

**Next goals ( Ver 1.1.0 ) :**

Peak detection and deconvolution (currently working on)  
Compound annotation tables  
Sample comparison (e.g., SL005 vs SL011)  
Batch pipelines for entire datasets  
Cleaner controller abstractions  
Auto-QC metrics (noise, intensity drift, retention time shifts)  

### Dataset using to test the package:
> **Little AS, Younker IT, Schechter MS, Bernardino PN, Méheust R, Stemczynski J, Scorza K, Mullowney MW, Sharan D, Waligurski E, Smith R, Ramanswamy R, Leiter W, Moran D, McMillin M, Odenwald MA, Iavarone AT, Sidebottom AM, Sundararajan A, Pamer EG, Eren AM, Light SH.**  
*Dietary- and host-derived metabolites are used by diverse gut bacteria for anaerobic respiration.*  
Nat Microbiol. 2024 Jan;9(1):55-69. Epub 2024 Jan 4.  
https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?task=365b6b62ef34444c821aff7d6ff9e499#%7B%22table_sort_history%22%3A%22main.collection_asc%22%2C%22main.collection_input%22%3A%22peak%7C%7CEXACT%22%7D

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
example = Ionome(run_id="SL011", samples="samples_SL011.yaml")
```

    >[2025-11-17 06:11:16][Ionome.__init__]
    	Initializing run ID SL011 for analysis

Then you can load the spectra data from mzml files using `load.data()`:

```python
example.load_data()
```

    >[2025-11-19 17:49:44][Ionome.load_data]
    	> Loading spectra data for sample SL011_EA_1...
    	  Using cached file
    	> Loading spectra data for sample SL011_MB_1...
    	  Using cached file
    	> Loading spectra data for sample SL011_HS_1...
    	  Using cached file
    	...


You can also check each individual dataframes by `.data['sample_name']`

```python
ea_1_df = example.data['SL011_EA_1'] # extraction blank control
hs_1_df = example.data['SL011_HS_1'] # sample
ea_1_df
hs_1_df
```
<style>
/* Hide radio buttons */
.tabs input[type="radio"] {
  display: none;
}

/* Tab labels (buttons) */
.tabs label {
  padding: 8px 16px;
  border: 1px solid #ccc;
  margin-right: -1px;
  cursor: pointer;
}

/* Active tab label styling */
.tabs input[type="radio"]:checked + label {
  font-weight: bold;
}

/* Content boxes */
.tab-content {
  border: 1px;
  padding: 15px;
  display: none;
  border-radius: 0 4px 4px 4px;
}

/* Display the active tab content */
#tab1:checked ~ #content1 {
  display: block;
}
#tab2:checked ~ #content2 {
  display: block;
}
</style>

<div class="tabs">

  <!-- TAB SELECTORS -->
  <input type="radio" id="tab1" name="tabs" checked>
  <label for="tab1">Extraction Blank</label>

  <input type="radio" id="tab2" name="tabs">
  <label for="tab2">Sample</label>

  <!-- TAB CONTENT 1 -->
  <div id="content1" class="tab-content">

<pre>
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
      <td>4.29355</td>
      <td>430.933655</td>
      <td>50.119999</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>1</td>
      <td>4.29355</td>
      <td>836.703247</td>
      <td>50.910000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>1</td>
      <td>4.29355</td>
      <td>1448.419556</td>
      <td>52.040001</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>1</td>
      <td>4.29355</td>
      <td>293.509064</td>
      <td>52.959999</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>1</td>
      <td>4.29355</td>
      <td>341.261475</td>
      <td>54.090000</td>
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
      <th>1889741</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990917</td>
      <td>4.432630</td>
      <td>595.700012</td>
    </tr>
    <tr>
      <th>1889742</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990917</td>
      <td>3.410859</td>
      <td>596.929993</td>
    </tr>
    <tr>
      <th>1889743</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990917</td>
      <td>0.733822</td>
      <td>597.760010</td>
    </tr>
    <tr>
      <th>1889744</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990917</td>
      <td>1.921110</td>
      <td>598.359985</td>
    </tr>
    <tr>
      <th>1889745</th>
      <td>1</td>
      <td>2990</td>
      <td>22.990917</td>
      <td>0.810828</td>
      <td>599.460022</td>
    </tr>
  </tbody>
</table>
</pre>

  </div>

  <!-- TAB CONTENT 2 -->
  <div id="content2" class="tab-content">

<pre>
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
      <td>4.293383</td>
      <td>23233.423828</td>
      <td>50.040001</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>1</td>
      <td>4.293383</td>
      <td>35742.230469</td>
      <td>51.040001</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>1</td>
      <td>4.293383</td>
      <td>91149.109375</td>
      <td>52.040001</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>1</td>
      <td>4.293383</td>
      <td>11119.541992</td>
      <td>53.029999</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>1</td>
      <td>4.293383</td>
      <td>1989.252319</td>
      <td>54.029999</td>
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
      <th>1735420</th>
      <td>1</td>
      <td>2990</td>
      <td>22.99075</td>
      <td>2.737132</td>
      <td>595.630005</td>
    </tr>
    <tr>
      <th>1735421</th>
      <td>1</td>
      <td>2990</td>
      <td>22.99075</td>
      <td>1.623667</td>
      <td>596.559998</td>
    </tr>
    <tr>
      <th>1735422</th>
      <td>1</td>
      <td>2990</td>
      <td>22.99075</td>
      <td>3.173591</td>
      <td>597.719971</td>
    </tr>
    <tr>
      <th>1735423</th>
      <td>1</td>
      <td>2990</td>
      <td>22.99075</td>
      <td>2.930570</td>
      <td>598.669983</td>
    </tr>
    <tr>
      <th>1735424</th>
      <td>1</td>
      <td>2990</td>
      <td>22.99075</td>
      <td>1.037758</td>
      <td>599.260010</td>
    </tr>
  </tbody>
</table>
</pre>

  </div>

</div>

### Quality Control plots
To plot and visualize basic quality control plots, use `plot_chromatogram('tic')` to show Total Ion Chromatograms (TIC plots)
```python
example.plot_chromatogram("tic")
```
<img src="git_images/SL011_EA_1_tic_chrom.png" width="450">
<img src="git_images/SL011_MB_1_tic_chrom.png" width="450">
<img src="git_images/SL011_HS_1_tic_chrom.png" width="450">


## LICENSE

MIT License — see LICENSE for details.

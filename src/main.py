import pandas as pd
import pymzml
from matplotlib import pyplot as plt

import numpy as np
from src.ionome_core import Ionome

from scipy.sparse.linalg import spsolve
from scipy import sparse


def main():
    # ----- Init -----
    pd.set_option('display.max_columns', None)
    first = Ionome(run_id="SL005", samples="samples_SL011.yaml")

    #----------------------------------------------------------------

    # ----- Load data -----
    first.load_data()

    # ----------------------------------------------------------------

    # ----- Quality control -----
    first.extract_quality_control()

    # __ TIC __
    plt.figure(figsize=(10, 6))
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.quality_control["retention_time"], sampleData.quality_control["tic"], alpha=0.5, label=sampleData.unique_id)

    plt.title(f"Total Ion Chromatograms - Overlay")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Total Ion Intensity")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # __ Peaks per scan - richness/noise __
    plt.figure(figsize=(10, 3))
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.quality_control['retention_time'], sampleData.quality_control['peaks_per_scan'], label = sampleData.unique_id)

    plt.title("Peaks per scan")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Number of peaks")
    plt.grid(True, linestyle="-", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # __ BPC & bpc_mz overlay __
    fig, ax = plt.subplots(2,1, figsize=(10,6), sharex=True)
    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            ax[0].plot(sampleData.quality_control["retention_time"], sampleData.quality_control['bpc'], label = sampleData.unique_id)
            ax[0].set_ylabel("Base peak intensity")
            ax[1].scatter(sampleData.quality_control["retention_time"], sampleData.quality_control["bpc_mz"], c=sampleData.quality_control["bpc"], cmap='viridis', s=6)
            ax[1].set_ylabel("Base peak m/z")

            ax[1].set_xlabel("Retention time (min)")
    # fig.colorbar(im, ax=ax[1], label='BPC intensity')
    plt.legend()
    plt.tight_layout()
    plt.show()
    # ----------------------------------------------------------------

    # ----- XIC -----
    first.extract_ion_chromatograms()

    for sampleData in first.samples.values():
        if sampleData.condition == "Treatment":
            plt.plot(sampleData.xic['ferulate']["retention_time"], sampleData.xic['ferulate']['intensity'], alpha=0.5,
                     label=sampleData.unique_id)

    plt.title(f"XIC - Overlay")
    plt.xlabel("Retention time (min)")
    plt.ylabel("Intensity")
    plt.legend()
    plt.tight_layout()
    plt.show()
if __name__ == "__main__":
    main() 
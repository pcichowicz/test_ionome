import pandas as pd
import pymzml
from matplotlib import pyplot as plt

import numpy as np
from src.ionome_core import Ionome

from scipy.sparse.linalg import spsolve
from scipy import sparse


def main():
    pd.set_option('display.max_columns', None)
    first = Ionome(run_id="SL005", samples="samples_SL011.yaml")
    first.load_data()
    for sample, metadata in first.samples.items():
        print(f"\t> SampleData object - {sample}\n")
        print(f"\t> {metadata.unique_id}")

    # first.extract_xic()
    # hs_1 = first.samples["SL011_HS_1"].xic["ferulate"]

    # def als(y, lam=1e5, p=0.01, itermax=50):
    #     r"""
    #
    #     Inputs:
    #         y:
    #             input data (i.e. chromatogram of spectrum)
    #         lam:
    #             parameter that can be adjusted by user. The larger lambda is,
    #             the smoother the resulting background, z
    #         p:
    #             wheighting deviations. 0.5 = symmetric, <0.5: negative
    #             deviations are stronger suppressed
    #         itermax:
    #             number of iterations to perform
    #     Output:
    #         the fitted background vector
    #
    #     """
    #     L = len(y)
    #     #  D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    #     D = sparse.eye(L, format='csc')
    #     D = D[1:] - D[:-1]  # numpy.diff( ,2) does not work with sparse matrix. This is a workaround.
    #     D = D[1:] - D[:-1]
    #     D = D.T
    #     w = np.ones(L)
    #
    #     for i in range(itermax):
    #         W = sparse.diags(w, 0, shape=(L, L))
    #         Z = W + lam * D.dot(D.T)
    #         z = spsolve(Z, w * y)
    #         w = p * (y > z) + (1 - p) * (y < z)
    #     return z
    #
    # from scipy.signal import savgol_filter
    #
    # def smooth_savgol(y, window=15, poly=3):
    #     """
    #     Savitzky–Golay smoothing.
    #     window: must be odd
    #     poly: polynomial order
    #     """
    #     if len(y) < window:
    #         window = max(5, len(y) // 2 * 2 + 1)  # enforce minimum odd window
    #     return savgol_filter(y, window_length=window, polyorder=poly, mode='mirror')
    #
    # iten_y = hs_1["intensity"]
    #
    # from scipy.ndimage import gaussian_filter1d
    #
    #
    # baseline = als(iten_y)
    #
    # baseline = np.maximum(baseline, 0.0)
    # corrected = hs_1["intensity"] - baseline
    # corrected[corrected < 0] = 0
    # smoothed = gaussian_filter1d(corrected, sigma=2, mode='reflect')
    # # smoothed = smooth_savgol(corrected)
    #
    # print("baseline min, max:", baseline.min(), baseline.max())
    # print("corrected min, max:", corrected.min(), corrected.max())
    # # How many corrected positives where raw==0?
    # zero_raw_pos = np.sum((hs_1["intensity"] == 0) & (corrected > 0))
    # print("zero_raw_became_positive:", zero_raw_pos)
    #
    # hs_1["baseline"] = baseline
    # hs_1["corrected"] = corrected
    # hs_1["smoothed"] = smoothed
    #
    # print(hs_1)
    # plt.figure(figsize=(12, 8))
    # plt.plot(hs_1["retention_time"], hs_1["intensity"], color='black',label = "raw")
    # plt.plot(hs_1["retention_time"], hs_1["baseline"], color = 'red', label = "baseline")
    # plt.plot(hs_1["retention_time"], hs_1["corrected"], color = 'blue', label = "corrected")
    # plt.plot(hs_1["retention_time"], hs_1["smoothed"], color='green', label = "smoothed")
    # plt.title("Ferulate")
    # plt.xlabel("retention time")
    # plt.ylabel("intensity")
    # plt.legend()
    # plt.show()

    # fig, axes = plt.subplots(2, 2, figsize=(12,8),sharey=True, sharex=True)
    # fig.suptitle('Ferulate')
    # axes = axes.flatten()
    #
    # chroms = ["intensity", "baseline", "corrected", "smoothed"]
    # for ax, chrom in zip(axes, chroms):
    #     ax.plot(hs_1["retention_time"], hs_1[chrom], color="black")
    #     ax.set_title(chrom)
    #     ax.set_xlabel("retention time")
    #     ax.set_ylabel("intensity")
    #
    # plt.tight_layout()
    # plt.show()

    # fig, axes = plt.subplots(2, 3, figsize = (12,6), sharex= True)
    #
    # axes = axes.flatten()
    #
    # for ax, (metabolite, df) in zip(axes, hs_1.items()):
    #     ax.plot(df["retention_time"], df["intensity"])
    #     ax.set_title(metabolite)
    #     ax.set_xlabel("retention time")
    #     ax.set_ylabel("intensity")
    #
    # plt.tight_layout()
    # plt.show()
if __name__ == "__main__":
    main() 
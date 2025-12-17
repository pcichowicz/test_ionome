import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from src.sampleData import SampleData

from matplotlib import image as mpimg


class Chromatograms:
    def __init__(self, sampledata: SampleData):
        self.sampleData = sampledata

    def _plot_chromatogram(self,
                           x:str,
                           y:str | None = None,
                           *,
                           data_df = None,
                           series: list[dict] | None = None,
                           ax = None,
                           title="",
                           xlabel="",
                           ylabel="",
                           **plot_kwargs):

        if ax is None:
            raise ValueError(f"ax must be provided")

        # remove plot_kwargs from args
        xlim = plot_kwargs.pop("xlim", None)
        ylim = plot_kwargs.pop("ylim", None)

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle="-", alpha=0.3)

        x_vals = data_df[x].values

        if series is None:
            y_vals = data_df[y].values
            ax.plot(x_vals, y_vals, color="black")
        else:
            for s in series:
                ax.plot(data_df[s["y"]].values,
                        label=s.get("label", s["y"]),
                        color=s.get("color", None),
                        linestyle=s.get("linestyle", "-"),
                        alpha=s.get("alpha", 1  ),
                        )
                ax.legend()

        if xlim:
            ax.xlim(xlim)
        if ylim:
            ax.ylim(ylim)

    def plot_tic(self, save_path:str | None = None, **plot_kwargs):
        fig, ax = plt.subplots(figsize=(10, 6))
        self._plot_chromatogram("retention_time",
                                "tic",
                                data_df = self.sampleData.quality_control,
                                ax = ax,
                                title= f"{self.sampleData.unique_id} - Total Ion Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Total Ion Intensity",
                                save_path= save_path,
                                **plot_kwargs
                                )
        plt.tight_layout()

        plt.show()
        plt.close(fig)

    def plot_bpc(self, save_path:str | None = None, **plot_kwargs):
        fig, ax = plt.subplots(figsize=(10, 6))
        self._plot_chromatogram("retention_time",
                                "bpc",
                                ax = ax,
                                data_df = self.sampleData.quality_control,
                                title= f"{self.sampleData.unique_id} - Base Peak Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Ion Intensity",
                                save_path= save_path,
                                **plot_kwargs
                                )
        plt.tight_layout()

        plt.show()
        plt.close(fig)

    def plot_corrected(self,plot_types = ("tic","bpc"),save_path:str | None = None, **plot_kwargs):
        fig, ax = plt.subplots(figsize=(10, 6))

        if isinstance (plot_types, str):
            plot_types = [plot_types]

        for plot_type in plot_types:

            series = [
                {
                    "y": plot_type,  # raw TIC
                    "label": "Raw",
                    "color": "black",
                    "alpha": 0.6,
                },
                {
                    "y": f"{plot_type}_baseline",  # baseline
                    "label": "Baseline",
                    "color": "red",
                    "linestyle": "--",
                },
                {
                    "y": f"{plot_type}_corrected",  # corrected
                    "label": "Corrected",
                    "color": "blue",
                },
            ]
            col = f"{plot_type}_corrected"
            if col not in self.sampleData.quality_control.columns:
                print(f"\t\033[33m Column {col} not found in sample data. Run 'correct_baseline('{plot_type}'), skipping.\033[0m")
                continue

            self._plot_chromatogram("retention_time",
                                    series = series,
                                    data_df = self.sampleData.quality_control,
                                    ax = ax,
                                    title= f"{self.sampleData.unique_id}\n Total Ion Chromatogram (Corrected)",
                                    xlabel= "Retention Time (min)",
                                    ylabel= "Total Ion Intensity",
                                    save_path= save_path,
                                    **plot_kwargs
                                    )
        plt.tight_layout()

        plt.show()
        plt.close(fig)

    def plot_xic(self,metabolites = "all",save_path:str | None = None, **plot_kwargs):
        if self.sampleData.xic is None:
            print(f"\t No XIC data available, Run 'extract_ion_chromatograms(), skipping.")
            return

        if metabolites == "all":
            metabolites = list(self.sampleData.xic.keys())

        for metabolite in metabolites:
            if metabolite not in self.sampleData.xic:
                print(f"\t XIC data for {metabolite} not found, skipping.")

            fig, ax = plt.subplots(figsize=(10, 6))
            self._plot_chromatogram("retention_time",
                                    "intensity",
                                    data_df =self.sampleData.xic[metabolite],
                                    ax = ax,
                                    title= f"{self.sampleData.unique_id} - {metabolite} - Extracted Ion Chromatogram",
                                    xlabel= "Retention Time (min)",
                                    ylabel= "Ion Intensity",
                                    save_path= save_path,
                                    **plot_kwargs
                                    )
            plt.tight_layout()

            plt.show()
            plt.close(fig)

    def plot_tic_and_bpc(self,save_path: str | None = None,**plot_kwargs):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

        self._plot_chromatogram(
            "retention_time",
            "tic",
            ax = ax1,
            data_df = self.sampleData.quality_control,
            title="Total Ion Chromatogram",
            xlabel="",  # Let bottom plot set x-label
            ylabel="Total ion intensity",
            show = False,
            **plot_kwargs
        )

        self._plot_chromatogram(
            "retention_time",
            "bpc",
            data_df = self.sampleData.quality_control,
            ax=ax2,
            title="Base Peak Chromatogram",
            xlabel="Retention time (min)",
            ylabel="Base peak intensity",
            show = False,
            **plot_kwargs
        )

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)

        plt.show()
        plt.close(fig)

    def plot_deconvolution(self,
                           time_col="retention_time",
                           signal_col="corrected",
                           show_individual=True,
                           show_area=True,
                           fig_size=(12,6)
                           ):

        time = self.data_df[time_col]
        raw_signal = self.data_df[signal_col]
        reconstructed_signal = self.unmixed_array.sum(axis=1)

        fig, ax = plt.subplots(figsize=fig_size)

        #raw
        ax.plot(time, raw_signal,label = "raw (baseline corrected")

        # reconstructed
        ax.plot(time, reconstructed_signal,label = "reconstructed", color = "red", linestyle="--")

        #indiv
        if show_individual:
            for i in range(self.unmixed_array.shape[1]):
                peak_signal = self.unmixed_array[:, i]
                ax.plot(time, peak_signal,color=_colors[i % len(_colors)],alpha = 0.7,label = f"peak {i+1}")
                if show_area:
                    ax.fill_between(time,0, peak_signal,alpha = 0.2, color = _colors[i % len(_colors)])

        ax.set_title(f"{self.sample_id} - Peak Deconvolution")
        ax.set_xlabel("Retention time (min)")
        ax.set_ylabel("Intensity")
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.legend(bbox_to_anchor=(1.05, 1),loc='upper left',ncol = 2)

        plt.tight_layout()

        plt.show()
        plt.close(fig)

    def compute_tic_bpc(self):
        df = self.data_df

        # TIC = sum of all intensities per scan
        tic = df.groupby("scan_id")["intensity"].sum().reset_index(name="tic")

        # BPC = max intensity per scan
        bpc = df.groupby("scan_id")["intensity"].max().reset_index(name="bpc")

        # single df with RT, TIC, BPC
        rt = df[["scan_id", "retention_time"]].drop_duplicates()

        merged = rt.merge(tic, on="scan_id").merge(bpc, on="scan_id")

        return merged.sort_values("retention_time")

_colors = sns.color_palette("deep", 15)
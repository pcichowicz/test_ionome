import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path

from matplotlib import image as mpimg


class Chromatograms:
    def __init__(self,sample_id, data_df,data_corrected,df_xic, unmixed_chromatogram_array: np.ndarray or None):
        self.data_df = data_df
        self.data_corrected = data_corrected
        self.df_xic = df_xic
        self.unmixed_array = unmixed_chromatogram_array
        self.sample_id = sample_id

    def _plot_chromatogram(self,
                           x,
                           y,
                           data_df = None,
                           ax=None,
                           fig_size=(10, 5),
                           title="",
                           xlabel="",
                           ylabel="",
                           save_path=None,
                           show = True,
                           **plot_kwargs):

        if data_df is None:
            data_df = self.data_df

        # remove plot_kwargs from args
        xlim = plot_kwargs.pop("xlim", None)
        ylim = plot_kwargs.pop("ylim", None)

        own_figure = False
        if ax is None:
            fig, ax = plt.subplots(figsize=fig_size)
            own_figure = True

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, linestyle="-", alpha=0.3)

        x_vals = data_df[x].values
        y_vals = data_df[y].values


        ax.plot(x_vals, y_vals, color = "black")

        plt.tight_layout()
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if save_path:
            save_path = Path(save_path)
            # save_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=300)

        if own_figure:
            if show:
                plt.show()
            plt.close()

        return save_path


    def plot_tic(self, save_path, **plot_kwargs):

        df = self.compute_tic_bpc()

        self._plot_chromatogram("retention_time",
                                "tic",
                                data_df = df,
                                title= f"{self.sample_id} - Total Ion Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Total Ion Intensity",
                                save_path= save_path,
                                **plot_kwargs
                                )

    def plot_corrected(self,save_path, **plot_kwargs):
        df = self.data_corrected
        self._plot_chromatogram("retention_time",
                                "corrected",
                                data_df = df,
                                title= f"{self.sample_id} - Total Ion Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Total Ion Intensity (corrected)",
                                save_path= save_path,
                                **plot_kwargs
                                )

    def plot_xic(self,save_path, **plot_kwargs):
        df = self.df_xic
        for metabolite, data_df in df.items():
            new_save_path = Path(str(save_path).replace(".png", f"_{metabolite}.png"))
            self._plot_chromatogram("retention_time",
                                    "intensity",
                                    data_df = data_df,
                                    title= f"{self.sample_id} - {metabolite} - Extracted Ion Chromatogram",
                                    xlabel= "Retention Time (min)",
                                    ylabel= "Ion Intensity",
                                    save_path= new_save_path,
                                    **plot_kwargs
                                    )

    def plot_tic_and_bpc(self,save_path,**plot_kwargs):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        self._plot_chromatogram(
            "retention_time", "tic",
            ax=ax1,
            title="Total Ion Chromatogram",
            xlabel="",  # Let bottom plot set x-label
            ylabel="Total ion intensity",
            **plot_kwargs
        )

        self._plot_chromatogram(
            "retention_time", "bpc",
            ax=ax2,
            title="Base Peak Chromatogram",
            xlabel="Retention time (min)",
            ylabel="Base peak intensity",
            **plot_kwargs
        )

        plt.tight_layout()

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
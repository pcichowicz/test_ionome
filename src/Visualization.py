import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


class Chromatograms:
    def __init__(self,sample_id, data_df, unmixed_chromatogram_array: np.ndarray or None):
        self.data_df = data_df
        self.unmixed_array = unmixed_chromatogram_array
        self.sample_id = sample_id

    def _plot_chromatogram(self,
                           x,
                           y,
                           ax=None,
                           fig_size=(10, 5),
                           title="",
                           xlabel="",
                           ylabel="",
                           save_path=None,
                           **plot_kwargs):

        # fig,ax = plt.subplots(figsize=fig_size)

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

        x_vals = self.data_df[x].values
        y_vals = self.data_df[y].values


        ax.plot(x_vals, y_vals, color = "black")

        plt.tight_layout()
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if own_figure:

            plt.show()
            plt.close()


    def plot_tic(self, **plot_kwargs):

        self._plot_chromatogram("retention_time",
                                "tic",
                                title= f"{self.sample_id} - Total Ion Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Total Ion Intensity",
                                **plot_kwargs
                                )

    def plot_corrected(self, **plot_kwargs):

        self._plot_chromatogram("retention_time",
                                "corrected",
                                title= f"{self.sample_id} - Total Ion Chromatogram",
                                xlabel= "Retention Time (min)",
                                ylabel= "Total Ion Intensity (corrected)",
                                **plot_kwargs
                                )

    def plot_tic_and_bpc(self,**plot_kwargs):
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

_colors = sns.color_palette("deep", 15)
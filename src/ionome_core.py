# load Ionome classes
import time

import pandas as pd

from src.SampleMetaData import SampleMetaData
from src.preprocess import MzmlParser
from src.correct_baseline import BaselineCorrection
from src.PeakDetection import PeakDetection
from src.Visualization import Chromatograms
from src.paths import CACHED_DIR
from time import sleep
from src.helpers import log_method_entry

# libraries
from pathlib import Path
import yaml
import pickle

class Ionome:
    """Encapsulates the mzML file analysis"""

    def __init__(self,
                 run_id: str,
                 samples: str | Path,
                 config_file: str ="config.yaml",
                 ):

        # Load yaml configurations
        with open(f"../config/{config_file}", "r") as f:
            self.config = yaml.safe_load(f)

        self.run_id: str = run_id
        self.sample_list = SampleMetaData(samples).meta_by_unique
        self.target_mz_list = self.config.get("target_mz_list", [])
        self.rerun: bool = self.config.get("rerun", False)
        self._method = self.config["baseline"].get("method", "asls")

        self.data: dict[str,pd.DataFrame] = {} # long form dataframe
        self.df_xic: dict[str, dict[str, pd.DataFrame]] = {} # extracted dataframe
        self.chrom_data: dict[str, pd.DataFrame] = {}

        self.data_window_prop = None
        self.data_peak_prop = None
        self.data_unmixed_chromatogram = None

        log_method_entry()
        print(f"\tInitializing run ID {self.run_id} for analysis")

    def __repr__(self):
        return (f"<Ionome>\n(run_id = {self.run_id}"
                f"\nsample_list = {self.sample_list}"
                f"\ndata = {self.data}"
                f"\ndata_window_prop = {self.data_window_prop}"
                f"\ndata_unmixed_chromatogram = {self.data_unmixed_chromatogram}"
                f"\ntarget_mz_list = {self.target_mz_list}"
                f"\nrerun = {self.rerun}"
                f"\nmethod = {self._method})")

    def __str__(self):
        report = " Analysis Report \n"

        # for key in self.__dict__.keys():
        #     report += f"{key, getattr(self, key)}\n"
        return report

    def load_data(self, **kwargs):
        """Loads the mzML file, will parse mzML file if parquet file is not already cached,
        Will save cached parquet file upon first parse of mzML file"""

        parser_cfg = self.config.get("parser", {})
        log_method_entry()

        for sample in self.sample_list:
            print(f"\t> Loading spectra data for sample {sample}...")
            sleep(1)

            parser = MzmlParser(self.sample_list[sample], rerun = self.rerun, **parser_cfg)
            df_data = parser.parse_or_load_mzml(**kwargs)
            self.data[sample] = df_data

    def extract_xic(self):

        extract_xic_cfg = self.config.get("target_mz_params", {})
        tol = extract_xic_cfg["tol"]
        log_method_entry()

        self.df_xic = {}
        xic_df_cache = CACHED_DIR / f"{self.run_id}_mz_target.pkl"

        # load cached if exists
        if xic_df_cache.exists() and not self.rerun:
            print(f"\t> Loading XIC cached data for sample(s) {list(self.sample_list.keys())}")
            with open(xic_df_cache, "rb") as f:
                self.df_xic = pickle.load(f)
            return self

        #Extract XIC data if no cached
        print(f"\t>Extracting mz targets for analysis...")

        for sample in self.sample_list:
            self.df_xic[sample] = {}

            for name, mz_value in self.target_mz_list.items():

                xic_df = self.data[sample][
                    (self.data[sample]['mz'] >= mz_value - tol) &
                    (self.data[sample]['mz'] <= mz_value + tol)
                ].copy()

                xic_df_all = self.data[sample][['scan_id', 'retention_time']].drop_duplicates()
                xic_final = xic_df_all.merge(xic_df[['scan_id', 'intensity']],
                                             on='scan_id',
                                             how='left')
                xic_final['intensity'] = xic_final['intensity'].fillna(0)
                xic_final = xic_final.sort_values("retention_time").reset_index(drop=True)

                self.df_xic[sample][name] = xic_final

        with open(xic_df_cache, "wb") as f:
            pickle.dump(self.df_xic, f)

        print(f">XIC extracted and cached")
        return self

    def correct_baseline(self):
        baseline_cfg = self.config.get("baseline", {})
        log_method_entry()

        print(f"\t> Correcting baseline using '{baseline_cfg['method']}'")
        for sample, metabolite in self.df_xic.items():

            for metabolite_name in metabolite:
                y_arr = metabolite[metabolite_name]['intensity']

                correction_params = baseline_cfg.get(self._method)

                baseline, corrected= BaselineCorrection().asls(y_arr, **correction_params)

                metabolite[metabolite_name]['baseline'] = baseline
                metabolite[metabolite_name]['corrected'] = corrected

        # for sample in self.sample_list:
        #     sleep(0.2)
        #     print(f"\t> Correcting baseline for sample {sample}")
        #
        #     quality_control_df = self.data[sample].groupby("scan_id").agg(
        #         retention_time = ("retention_time", "first"),
        #         tic = ("intensity", "sum"),
        #         bpc = ("intensity", "max")
        #     )
        #
        #     # xic_df = self.data[sample]
        #
        #     if method == "asls" or "snip":
        #         self._method = method
        #         correction_parameters = baseline_cfg.get(self._method, {})
        #     else:
        #         print(f"\tMethod {method} not supported, defaulting to 'snip'.")
        #         self._method = baseline_cfg.get("method")
        #         correction_parameters = baseline_cfg.get(self._method, {})
        #
        #     if self._method == "snip":
        #         baseline, corrected = BaselineCorrection().snip(quality_control_df, **correction_parameters)
        #     elif self._method == "asls":
        #         baseline, corrected = BaselineCorrection().asls(quality_control_df, **correction_parameters)
        #     else:
        #         raise ValueError("Invalid method, only 'snip' or 'asls'")
        #
        #     quality_control_df["baseline"] = baseline
        #     quality_control_df["corrected"] = corrected
        #
        #     self.chrom_data[sample] = quality_control_df

    def peak_detection(self):

        self.data_window_prop = {}
        self.data_peak_prop = {}
        self.data_unmixed_chromatogram = {}

        log_method_entry()
        for sample in self.sample_list:

            cache_file = CACHED_DIR / f"{sample}_peak_detection.pkl"

            if cache_file.exists() and not self.rerun:
                print(f"\t> Cached peak detection results already exist for sample {sample}... skipping")

                with open(cache_file, "rb") as f:
                    cache_data = pickle.load(f)

                (self.data_window_prop[sample],
                 self.data_peak_prop[sample],
                 self.data_unmixed_chromatogram[sample]
                 ) = cache_data
                time.sleep(0.5)
                continue

            print(f"\t> Peak detection for sample {sample}...")
            sleep(1)

            peak_detection_cfg = self.config.get("peak_detection", {})
            # print(self.data_wins)
            peaks = PeakDetection(self.chrom_data[sample], **peak_detection_cfg)
            # print(self.data_wins)
            # print(peaks.new_df)
            window_prop, peak_prop, unmixed_chromatogram = peaks.detect_peaks()
            # print(peaks.new_df)

            self.data_window_prop[sample] = window_prop
            self.data_peak_prop[sample] = peak_prop
            self.data_unmixed_chromatogram[sample] = unmixed_chromatogram

            with open(cache_file, "wb") as file:
                pickle.dump(
                    (window_prop, peak_prop, unmixed_chromatogram),
                    file,
                    protocol=pickle.HIGHEST_PROTOCOL
                )

        return self

    def peak_detection_xic(self):
        results = []
        for sample, mets in self.df_xic.items():
            for metabolite, df in mets.items():

                peak = PeakDetection(df)
                peaks = peak.detect_xic_peaks(self.df_xic[sample][metabolite])
                if peaks is None:
                    results.append({
                        "sample": sample,
                        "metabolite": metabolite,
                        "peak_found": False
                    })
                else:
                    results.append({
                        "sample": sample,
                        "metabolite": metabolite,
                        "peak_found": True,
                        **peaks
                    })

        results_df = pd.DataFrame(results)
        print(results_df)

    def plot_chromatogram(self, type_plot: str = "tic"):
        """
        Default is TIC plot.
        """

        chrom_plots: list = []

        for sample in self.sample_list:

            # plot_chromatogram = Chromatograms(sample,self.chrom_data[sample], self.data_unmixed_chromatogram[sample])
            plot_chromatogram = Chromatograms(sample,self.data, self.data_unmixed_chromatogram[sample])
            plot_map = {
                "tic": plot_chromatogram.plot_tic,
                "corrected": plot_chromatogram.plot_corrected,
                "tic_bpc": plot_chromatogram.plot_tic_and_bpc,
                "decon": plot_chromatogram.plot_deconvolution,

            }
            chrom_plots.append(plot_map[type_plot]())

        return chrom_plots



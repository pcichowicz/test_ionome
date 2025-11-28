# load Ionome classes
import time

import pandas as pd

from src.sampleData import SampleData
from src.preprocess import MzmlParser
from src.correct_baseline import BaselineCorrection
from src.PeakDetection import PeakDetection
from src.Visualization import Chromatograms
from src.paths import CACHED_DIR, RESULTS_DIR
from time import sleep
from src.helpers import log_method_entry

# libraries
from pathlib import Path
import yaml
import pickle

class Ionome:
    """
    Encapsulates the mzML file analysis
    orchestrator class
    """

    def __init__(self,
                 run_id: str,
                 samples: str | Path,
                 config_file: str ="config.yaml",
                 ):
        log_method_entry()

        # Load yaml configurations
        with open(f"../config/{config_file}", "r") as f:
            self.config = yaml.safe_load(f)

        self.run_id: str = run_id

        # config yaml settings
        self.target_mz_list = self.config.get("target_mz_list", [])
        self.rerun: bool = self.config.get("rerun", False)
        self._method = self.config["baseline"].get("method", "asls")

        # sample yaml metadata
        self.sample_metadata = self._load_sample_yaml(samples)

        # SampleData objects for each sample
        # {unique id: SampleData}
        self.samples = {}
        for uid, meta in self.sample_metadata.items():
            sample = SampleData(
                unique_id = uid,
                batch_id = meta["id"],
                file = meta["file"],
                condition = meta["condition"],
                description=meta["description"],
                replicate = meta["replicate"],
                species = meta["species"]
                )
            self.samples[uid] = sample

        print(f"\tInitializing run ID {self.run_id} with {len(self.samples)} samples:")
        for s in self.samples.values():
            print(f"\t → {s.unique_id} ({s.description} | rep {s.replicate} | species {s.species})")

    def _load_sample_yaml(self,file):
        with open(f"../config/{file}", "r") as f:
            y = yaml.safe_load(f)

        meta_by_unique_id = {}
        for sample in y["samples"]:
            unique_id = sample["unique_id"]
            meta_by_unique_id[unique_id] = sample
        return meta_by_unique_id

    def load_data(self, **kwargs):
        """Loads the mzML file, will parse mzML file if parquet file is not already cached,
        Will save cached parquet file upon first parse of mzML file"""
        log_method_entry()

        parser_cfg = self.config.get("parser", {})

        for uid, sampleData in self.samples.items():
            print(f"\t> Loading spectra data for sample {uid} ...")
            # sleep(0.2)

            parser = MzmlParser(sampleData.file, rerun=self.rerun, **parser_cfg)
            sampleData.raw = parser.parse_or_load_mzml(**kwargs)

    def extract_quality_control(self):
        log_method_entry()

        for uid, sampleData in self.samples.items():
            print(f"\t> Extracting TIC BPC data for sample {uid} ...")
            sampleData.qc_df()

    def extract_ion_chromatograms(self, tolerance_type: str = "tol"):
        log_method_entry()

        extract_xic_cfg = self.config.get("target_mz_params", {})

        # low resolution -> 'tol'
        # high resolution -> 'ppm'
        if tolerance_type == "tol":
            tolerance = extract_xic_cfg["tol"]
        else:
            tolerance = extract_xic_cfg["ppm"]

        for uid, sampleData in self.samples.items():
            print(f"\t> Extracting XIC chromatogram data for sample {uid} ...")
            sampleData.xic_df(target_list=self.target_mz_list, tol=tolerance)

    def correct_baseline(self):
        baseline_cfg = self.config.get("baseline", {})
        log_method_entry()

        for sample, ddf in self.data.items():
            print(f"\t> Correcting baseline for sample {sample}...")
            sleep(1)

            skeleton_df = ddf[['scan_id', 'retention_time']].drop_duplicates()
            tic = ddf.groupby('scan_id')['intensity'].sum().reset_index(name='tic')
            y = tic['tic'].to_numpy()
            correction_params = baseline_cfg.get(self._method)
            baseline, corrected = BaselineCorrection().asls(y, **correction_params)
            clean_df = skeleton_df.merge(tic, on='scan_id')
            clean_df['baseline'] = baseline
            clean_df['corrected'] = corrected

            self.data_corrected[sample] = clean_df

        # print(f"\t> Correcting baseline using '{baseline_cfg['method']}'")
        # for sample, metabolite in self.df_xic.items():
        #     print(sample, metabolite)
        #     for metabolite_name in metabolite:
        #         print(metabolite_name, metabolite[metabolite_name])
                # tic = df.groupby("scan_id")["intensity"].sum().reset_index(name="tic")
                # y_arr = metabolite[metabolite_name]['intensity']

                # correction_params = baseline_cfg.get(self._method)
                #
                # baseline, corrected= BaselineCorrection().asls(y_arr, **correction_params)
                #
                # metabolite[metabolite_name]['baseline'] = baseline
                # metabolite[metabolite_name]['corrected'] = corrected

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

        output_path = {}
        chrom_plots: list = []

        for sample in self.sample_list:
            print(f"\t> Plotting chromatogram for sample {sample}...")

            chrom = Chromatograms(
                sample,
                data_df=self.data[sample],
                data_corrected=self.data_corrected[sample],
                df_xic = self.df_xic[sample],
                unmixed_chromatogram_array=self.data_unmixed_chromatogram,
            )

            plot_map = {
                "tic": chrom.plot_tic,
                "corrected": chrom.plot_corrected,
                "tic_bpc": chrom.plot_tic_and_bpc,
                "decon": chrom.plot_deconvolution,
                "xic": chrom.plot_xic,

            }

            out_dir = RESULTS_DIR / f"{sample}_{type_plot}_chrom.png"

            chrom_plots.append(plot_map[type_plot](out_dir))

        return chrom_plots

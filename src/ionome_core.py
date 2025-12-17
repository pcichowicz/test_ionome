# load Ionome classes
import time

import pandas as pd

from src.sampleData import SampleData
from src.preprocess import MzmlParser
from src.correct_baseline import BaselineCorrection
from src.PeakDetection import PeakDetection
from src.Visualization import Chromatograms
from src.paths import output_path
from time import sleep
from src.helpers import log_method_entry
from src.detectPeaks import DetectPeaks

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
                 samples: str | Path
                 ):
        log_method_entry()

        config_file = Path(__file__).parent.parent / "projects" / run_id / f"config_{run_id}.yaml"
        # Load yaml configurations
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        self.run_id: str = run_id
        self.project_path: Path = Path(__file__).parent.parent / "projects" / run_id
        self.raw_data = self.project_path / "raw_data"

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

        print(f"\t> Initializing run ID {self.run_id} with {len(self.samples)} samples:")
        for s in self.samples.values():
            print(f"\t  → {s.unique_id} ({s.description} | rep {s.replicate} | species {s.species})")

    def _load_sample_yaml(self,file):
        with open(file, "r") as f:
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

        print(f"\t> Loading spectra data:")

        for uid, sampleData in self.samples.items():
            mzml_path = self.raw_data / sampleData.file
            parser = MzmlParser(mzml_path,run_id=self.run_id, rerun=self.rerun, **parser_cfg)

            sampleData.raw = parser.parse_or_load_mzml(**kwargs)

    def extract_quality_control(self):
        log_method_entry()

        print(f"\t> Extracting quality control data (TIC,BPC):")
        for uid, sampleData in self.samples.items():

            sampleData.qc_df()

    def extract_ion_chromatograms(self, tolerance_type: str = "ppm"):
        log_method_entry()

        extract_xic_cfg = self.config.get("target_mz_params", {})

        # low resolution -> 'da'
        # high resolution -> 'ppm'
        if tolerance_type == "da":
            tolerance = extract_xic_cfg["da"]
        else:
            tolerance = extract_xic_cfg["ppm"]

        print(f"\t> Extracting XIC chromatogram data:")
        for uid, sampleData in self.samples.items():

            sampleData.xic_df(target_list=self.target_mz_list, tol=tolerance, tol_type=tolerance_type)

    def correct_baseline(self, chromatogram: str):
        baseline_cfg = self.config.get("baseline", {})
        correction_params = baseline_cfg.get(self._method)
        log_method_entry()

        bc = BaselineCorrection()
        print(f"\t> Correcting {chromatogram} chromatogram baseline using '{self._method}' method:")
        for uid, sampleData in self.samples.items():

            chroms = sampleData.get_chromatograms(chromatogram)
            for name, df in chroms.items():
                # print(f"\t> Correcting {name} chromatogram")

                if chromatogram == "xic":
                    col_name = "intensity"
                else:
                    col_name = chromatogram
                baseline, corrected = bc.asls(df[col_name], **correction_params)

                df = df.copy()
                df["baseline"] = baseline
                df["corrected"] = corrected

                sampleData.add_chromatograms(chromatogram, name, df)
            print(f"\t \033[32m ✓ \033[0m{sampleData.unique_id}")

    def peak_d(self):
        log_method_entry()
        peak_detect_cfg = self.config.get("peak_detection", {})

        print(f"\t> Detecting peaks in XIC Chromatograms:")
        for uid, sampleData in self.samples.items():
            if uid == "SL2031_EL-cat_MS1_neg_1":
                for metabolite, xic_df in sampleData.xic.items():
                    peak_detector = DetectPeaks(metabolite, xic_df, **peak_detect_cfg)
                    if metabolite == "catechin":
                        # print(f"\t> Peak detection for sample {sampleData.unique_id} | metabolite '{metabolite}'")
                        peak_detector.detect_peaks()
                        window_df_props, peak_properties, unmixed_chromatogram = peak_detector.detect_peaks()
                        # print(window_df_props, peak_properties, unmixed_chromatogram, sep="\n\n")
                        sampleData.window_df_properties[metabolite] = window_df_props
                        sampleData.peaks_properties[metabolite] = peak_properties
                        sampleData.unmixed_chromatograms[metabolite] = unmixed_chromatogram

    def plot_chromatogram(self, type_plot: str = "tic"):
        """
        Default is ....?.
        """
        log_method_entry()
        print(f"Plotting Chromatograms")
        plot_cfg = self.config.get("plotting", {})
        plotting_params = plot_cfg.get("plotting_params", {})

        for uid, sampleData in self.samples.items():
            if uid != "SL2031_EL-cat_MS1_neg_1":
                continue

            chrom = Chromatograms(sampleData)

            plot_map = {
                "tic": chrom.plot_tic,
                "bpc": chrom.plot_bpc,
                "tic_and_bpc": chrom.plot_tic_and_bpc,
                "corrected": chrom.plot_corrected,
                "decon": chrom.plot_deconvolution,
                "xic": chrom.plot_xic,

            }

            for plot, enabled in plot_cfg.items():
                if not enabled:
                    continue

                plot_func = plot_map.get(plot)
                if plot_func is None:
                    print(f"\t\033[32m Unknown plot \033[0m{plot}, skipping \033[0m")
                    continue

                print(f"\t Plotting {plot} for {sampleData.unique_id}")

                params = plotting_params.get(plot, {})
                plot_func(**params)

            # for plot in plot_cfg.keys():
            #     if plot_cfg[plot] is True: # next is to also add in check if data is available/previous methods ran
            #         print(f"\t> Plotting {plot} for {sampleData.unique_id}")
            #         plot_map[plot]()
















        # plot_cfg["type"] = type_plot

        # for uid, sampleData in self.samples.items():
        #     if uid == "SL2031_EL-cat_MS1_neg_1":
        #
        #         chromatogram = Chromatograms(sampleData)
        #
        #         chromatogram.plot_tic()
        #         chromatogram.plot_bpc()
        #         chromatogram.plot_corrected(plot_type=type_plot)
        #         chromatogram.plot_tic_and_bpc()
        #
        #         for metabolite in chromatogram.sampleData.xic.keys():
        #             print(metabolite)
        #             chromatogram.plot_xic(metabolite)


        # output_path = {}
        # chrom_plots: list = []
        #
        # for sample in self.sample_list:
        #     print(f"\t> Plotting chromatogram for sample {sample}...")
        #
        #     chrom = Chromatograms(
        #         sample,
        #         data_df=self.data[sample],
        #         data_corrected=self.data_corrected[sample],
        #         df_xic = self.df_xic[sample],
        #         unmixed_chromatogram_array=self.data_unmixed_chromatogram,
        #     )
        #
        #     plot_map = {
        #         "tic": chrom.plot_tic,
        #         "corrected": chrom.plot_corrected,
        #         "tic_bpc": chrom.plot_tic_and_bpc,
        #         "decon": chrom.plot_deconvolution,
        #         "xic": chrom.plot_xic,
        #
        #     }
        #
        #     out_dir = output_path(self.run_id,"results_dir") / f"{sample}_{type_plot}_chrom.png"
        #
        #     chrom_plots.append(plot_map[type_plot](out_dir))



# class RunPaths:
#     def __init__(self, run_id: str, base: Path = Path("projects")
#         self.config = config
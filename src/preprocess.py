"""
Class method for parsing mzml file format of LCMS data.
Returns dataframe
"""
import pandas as pd
import pymzml
from src.paths import output_path
from pathlib import Path

class MzmlParser:
    def __init__(self,
                 mzml_file: str | Path,
                 rerun,
                 run_id,
                 **parser_cfg):
        self.mzml_file = mzml_file
        self._rerun: bool = rerun
        self._settings = parser_cfg
        self.run_id = run_id

    def parse_mzml_file(self, **kwargs):

        rows = []
        reader = pymzml.run.Reader(self.mzml_file)
        for spec in reader:
            if spec.ms_level != self._settings.get("ms_level"):
                continue

            rt = spec.scan_time_in_minutes()
            scan_id = spec.ID

            for intensity, mz in zip(spec.i, spec.mz):
                rows.append({
                    "ms_level": spec.ms_level,
                    "scan_id": scan_id,
                    "retention_time": rt,
                    "intensity": intensity,
                    "mz": mz,
                })
        scans_df = pd.DataFrame(rows)

        return scans_df

    def parse_or_load_mzml(self,**kwargs):
        """
        Parse mzML file into a master DataFrame or load cached Parquet version.

        Returns
        -------
        pd.DataFrame
            Parsed chromatogram and peak data
        """
        parquet_path = output_path(self.run_id, "cached_dir") / self.mzml_file.name.replace(".mzML",".parquet")

        # parquet_path = self.mzml_file.with_suffix(".parquet")
        # If cached file exists and rerun=False → load from parquet

        if parquet_path.exists() and not self._rerun:
            print(f"\t \033[32m ✓ \033[0m{self.mzml_file.name}")
            return pd.read_parquet(parquet_path)


        # Otherwise, parse mzML → DataFrame
        # print(f"\t Parsing {self.mzml_file.name} ... ")
        master_df = self.parse_mzml_file(**kwargs)
        print(f"\t ✓ {self.mzml_file.name}")

        # Save to Parquet, reuse if 'rerun' is False
        master_df.to_parquet(parquet_path, index=False)

        return master_df
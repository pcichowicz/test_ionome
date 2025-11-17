"""
Class method for parsing mzml file format of LCMS data.
Returns dataframe
"""
import pandas as pd
import pymzml
from src.paths import DATA_DIR, CACHED_DIR

class MzmlParser:
    def __init__(self, mzml_file: dict, rerun, **parser_cfg):
        self.mzml_file = mzml_file["file"]
        self._rerun: bool = rerun
        self._settings = parser_cfg

    def parse_mzml_file(self, **kwargs):
        mzml_file_path = DATA_DIR / self.mzml_file
        rows = []
        reader = pymzml.run.Reader(mzml_file_path)
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
        parquet_path = CACHED_DIR / self.mzml_file.replace(".mzML", ".parquet")

        # If cached file exists and rerun=False → load from parquet
        if parquet_path.exists() and not self._rerun:
            print(f"\t Using cached file")
            return pd.read_parquet(parquet_path)

        # Otherwise, parse mzML → DataFrame
        print(f"\t Parsing {self.mzml_file} ... ")
        master_df = self.parse_mzml_file(**kwargs)

        # Save to Parquet, reuse if 'rerun' is False
        print(f"\t Caching spectra data")
        master_df.to_parquet(parquet_path, index=False)

        return master_df
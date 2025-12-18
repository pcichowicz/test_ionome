#sampleData.py
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Any, Optional

# from src.scratch import sampleData


def format_data_summary(data: pd.DataFrame | np.ndarray | Dict[str, Any] | None) -> Optional[str]:
    """
    Formats DataFrame/ndarray data into clean summary lines (top 5 values).
    Returns a list of lines (strings) to be joined/indented later.
    Returns None if input data is None or empty.
    """
    if data is None:
        return None

    summary_lines = []
    indent_size = f"\t\t"
    if isinstance(data, pd.DataFrame):
        if data.empty:
            return None # Indicate empty

        df_head = data.head(5)
        for col_name in data.columns:
            values = df_head[col_name].tolist()
            values_str = str(tuple(values)).replace("'", "").replace('"', '')
            # Return just the line item format, indentation handled by the caller
            summary_lines.append(f"{col_name:<20}: {values_str}")
        return "\n".join(summary_lines)

    elif isinstance(data, np.ndarray):
        if data.size == 0:
            return None # Indicate empty
        values = data.flatten()[:5].tolist()
        values_str = str(tuple(values)).replace("'", "").replace('"', '')
        # Return just the line item format, indentation handled by the caller
        summary_lines.append(f"Values (first 5)     : {values_str}")
        return "\n".join(summary_lines)

    elif isinstance(data, dict):
        # If a dictionary reaches this function unexpectedly, we return None and let the dict handler manage it
        return None

    else:
        return f"Unsupported type: {type(data).__name__}"

# -----------------------------------------------------------
# Helper 2: Handles nested Dictionaries recursively
# -----------------------------------------------------------
def format_dict_structure(data_dict: Dict[str, Any], indent_level: int = 1) -> str:
    """Recursively formats a dictionary structure."""
    if data_dict is None or not isinstance(data_dict, dict):
        return f"\t\t Input is not a valid dictionary or is None"

    if not data_dict:
        return f"\t\t Empty dictionary"

    summary_lines = []
    indent_space = f"\t" * indent_level

    for key, value in data_dict.items():
        if value is None:
            summary_lines.append(f" [{key}] (Value: None)")
            continue

        if isinstance(value, (pd.DataFrame, np.ndarray)):
            # This is a LEAF NODE (DF or NDArray)
            shape_str = f"Shape: {value.shape}" if hasattr(value, 'shape') else ""

            # Print the key and shape on the SAME LINE
            summary_lines.append(f"{indent_space} [{key}] ({shape_str})")
                                 ##

            # Get the *content* lines from the helper function
            content_summary = format_data_summary(value)

            if content_summary:
                # Add extra indentation to the content lines
                indented_content = content_summary.replace("\n", f"\n{indent_space}")
                summary_lines.append(f" {indented_content}")

        elif isinstance(value, dict):
            # This is a NESTED DICT (recurse)
            summary_lines.append(f" [{key}] ({len(value)} entries)")
            nested_summary = format_dict_structure(value, indent_level + 1)
            summary_lines.append(nested_summary)

        else:
            # Handle simple metadata values within a dict
            summary_lines.append(f" [{key}] (Value: {str(value)[:50]}...)")

    return "\n".join(summary_lines)

# def format_dict_of_dataframes(data_dict: Dict[str, pd.DataFrame]) -> str:
#     """Formats a dictionary where values are DataFrames."""
#     if not data_dict:
#         return "None"
#
#     summary_lines = []
#     for key, df in data_dict.items():
#         summary_lines.append(f"\t\tKey: '{key}' (Shape: {df.shape})")
#         # Indent the DataFrame summary further
#         df_summary = format_data_summary(df).replace("\t\t- ", "\t\t    - ")
#         summary_lines.append(df_summary)
#
#     return "\n".join(summary_lines)

@dataclass
class SampleData:
    # ----- Identifiers -----
    unique_id: str | None = None
    batch_id: str | None = None
    file: Path | str | None = None
    condition: str | None = None
    description: str | None = None
    replicate: int | None = None
    species: str | None = None

    # ----- Data containers -----
    raw: pd.DataFrame | None = None
    quality_control: pd.DataFrame | None = None
    # baseline_corrected: dict[str, pd.DataFrame] = field(default_factory=dict)   # This can move to quality_control??
    xic: dict[str, pd.DataFrame] = field(default_factory=dict)
    unmixed_chromatograms: dict[str, pd.DataFrame] = field(default_factory=dict)
    peaks_properties: dict[str, pd.DataFrame] = field(default_factory=dict)
    window_df_properties: dict[str, pd.DataFrame] = field(default_factory=dict)

    def summarize(self):
        """Prints a structured summary of the SampleData object contents."""
        print("-" * 60)
        print(f"Sample Summary for sample: {self.unique_id or 'N/A'}")
        print("-" * 60)

        # 1. Print Metadata
        print("\t--- Metadata ---")
        for field_name in ['batch_id', 'condition', 'species', 'replicate']:
            value = getattr(self, field_name)
            if value is not None:
                print(f" {field_name:<12}: {value}")

        # 2. Print Data Containers
        print("\n\t--- Data Containers & Shapes ---")

        # Define the fields that are simple DataFrames
        simple_dfs = ['raw', 'quality_control']
        for field_name in simple_dfs:
            df = getattr(self, field_name)
            print(f"\t .{field_name:<20} (Shape: {df.shape if df is not None else 'None'})")
            if df is not None and not df.empty:
                print(f"\t\t {format_data_summary(df)}")

        # Define the fields that are Dictionaries of DataFrames
        dict_dfs = ['xic',
                    'unmixed_chromatograms',
                    'peaks_properties', 'window_df_properties']

        for field_name in dict_dfs:
            data_dict = getattr(self, field_name)
            count = len(data_dict)
            print(f"\n .{field_name:<20} (Contains {count} entries)")
            if count > 0:
                print(format_dict_structure(data_dict))

    def get_chromatograms(self, chromatogram: str) -> dict[str, pd.DataFrame]:

        if chromatogram == 'xic':
            return self.xic.copy()

        elif chromatogram == 'tic':
            return {"tic": self.quality_control[["tic"]].copy()}

        elif chromatogram == "bpc":
            return {"bpc": self.quality_control[["bpc"]].copy()}

        else:
            raise ValueError(f"Unknown chromatogram type: {chromatogram}")

    def add_chromatograms(self,chromatogram: str, key: str, df: pd.DataFrame):
        if chromatogram == 'xic':
            self.xic[key] = df

        elif chromatogram in ('tic', 'bpc'):
            if chromatogram == 'tic':
                self.quality_control["tic_baseline"] = df["baseline"]
                self.quality_control["tic_corrected"] = df["corrected"]
            else:
                self.quality_control["bpc_baseline"] = df["baseline"]
                self.quality_control["bpc_corrected"] = df["corrected"]

    def qc_df(self):

        quality_control_base = self.raw[['scan_id', 'retention_time']].drop_duplicates()
        tic = self.raw.groupby('scan_id')['intensity'].sum().reset_index(name='tic')
        bpc = self.raw.groupby('scan_id')['intensity'].max().reset_index(name='bpc')
        pps = self.raw.groupby('scan_id').size().rename('peaks_per_scan').reset_index()
        bpc_mz = (self.raw.loc[self.raw.groupby('scan_id')['intensity'].idxmax(), ['scan_id', 'mz']].reset_index(drop=True).rename(columns={'mz': 'bpc_mz'}))

        self.quality_control = (
            quality_control_base
            .merge(tic, on='scan_id')
            .merge(bpc, on='scan_id')
            .merge(pps, on='scan_id')
            .merge(bpc_mz, on='scan_id')
        )
        print(f"\t \033[32m ✓ \033[0m{self.unique_id}")

    def xic_df(self, target_list, tol, tol_type):
        print(f"\t  {self.unique_id} --> {target_list}")
        for metabolite, mz_value in target_list.items():
            # print(f"\t * Extracting metabolite target {metabolite} with mz of {mz_value}")
            if tol_type == "da":
                tol = tol
            elif tol_type == "ppm":
                tol = tol * mz_value / 1e6

            xic_df_mz = self.raw[
                (self.raw["mz"] >= mz_value - tol) &
                (self.raw["mz"] <= mz_value + tol)].copy()

            if len(xic_df_mz) > 0:
                print(f"\t\t \033[32m ✓ \033[0m{metabolite} ({mz_value}): tol={tol:.6f} Da, hits={len(xic_df_mz)}")
            else:
                print(f"\t\t \033[31m x \033[0m{metabolite} ({mz_value}): tol={tol:.6f} Da, hits={len(xic_df_mz)}")

            xic_df_base = self.raw[['scan_id', 'retention_time']].drop_duplicates()

            xic_df = xic_df_base.merge(xic_df_mz[['scan_id', 'intensity']], on='scan_id', how='left')
            xic_df['intensity'] = xic_df['intensity'].fillna(0)
            xic_df = xic_df.sort_values('retention_time').reset_index(drop=True)

            self.xic[metabolite] = xic_df


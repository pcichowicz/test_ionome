#sampleData.py
from dataclasses import dataclass, field
import pandas as pd
from pathlib import Path

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
    baseline_corrected: pd.DataFrame | None = None
    xic: dict[str, pd.DataFrame] = field(default_factory=dict)
    chromatograms: dict[str, pd.DataFrame] = field(default_factory=dict)
    peaks: dict[str, pd.DataFrame] = field(default_factory=dict)

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

    def xic_df(self, target_list, tol):

        for metabolite, mz_value in target_list.items():
            # print(f"\t * Extracting metabolite target {metabolite} with mz of {mz_value}")
            xic_df_mz = self.raw[
                (self.raw["mz"] >= mz_value - tol) &
                (self.raw["mz"] <= mz_value + tol)].copy()
            print(f"\t * {metabolite} ({mz_value}): tol={tol:.6f} Da, hits={len(xic_df_mz)}")
            # print(xic_df_mz[['mz', 'intensity']].head())
            xic_df_base = self.raw[['scan_id', 'retention_time']].drop_duplicates()

            xic_df = xic_df_base.merge(xic_df_mz[['scan_id', 'intensity']], on='scan_id', how='left')
            xic_df['intensity'] = xic_df['intensity'].fillna(0)
            xic_df = xic_df.sort_values('retention_time').reset_index(drop=True)

            self.xic[metabolite] = xic_df


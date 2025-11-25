#sampleData.py
from dataclasses import dataclass, field
import pandas as pd
from pathlib import Path

@dataclass
class SampleData:
    name: str
    path: Path | None = None

    raw: pd.DataFrame | None = None
    baseline_corrected: pd.DataFrame | None = None

    xic: dict[str, pd.DataFrame] = field(default_factory=dict)
    chromatograms: dict[str, pd.DataFrame] = field(default_factory=dict)
    peaks: dict[str, pd.DataFrame] = field(default_factory=dict)

    def has_raw(self): return self.raw is not None
    def has_baseline(self): return self.baseline_corrected is not None
    def has_xic(self): return bool(self.xic)
import yaml
from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd

class SampleMetaData:
    def __init__(self, yaml_file: str | Path):

        with open(f"../config/{yaml_file}", "r") as f:
            data = yaml.safe_load(f)

        self.df = pd.DataFrame(data["samples"])
        self.sample_list = self.df["file"].tolist()

        self.meta_by_unique = self.df.set_index("unique_id").to_dict(orient="index")
        # self.meta_by_sample = (
        #     self.df.groupby("id")
        #     .apply(lambda g: g.to_dict(orient="records"))
        #     .to_dict()
        # )

        if self.df["unique_id"].duplicated().any():
            raise ValueError("Duplicate unique_id values found in sample metadate YAML!")

    def __repr__(self):
        return f"<SampleMetaData>:(sample_list={self.sample_list})"

@dataclass
class SampleData:
    unique_id: str
    path: Path | str | None = None

    #metadata
    batch_id: str | None = None
    condition: str | None = None
    replicate:str | None = None
    description: str | None = None

    meta:  dict = field(default_factory=dict)

    raw: pd.DataFrame | None = None
    baseline_corrected: pd.DataFrame | None = None

    xic: dict[str, pd.DataFrame] = field(default_factory=dict)
    chromatograms: dict[str, pd.DataFrame] = field(default_factory=dict)
    peaks: dict[str, pd.DataFrame] = field(default_factory=dict)

    def has_raw(self): return self.raw is not None
    def has_baseline(self): return self.baseline_corrected is not None
    def has_xic(self): return bool(self.xic)
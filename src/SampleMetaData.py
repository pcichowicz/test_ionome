import yaml
from pathlib import Path
import pandas as pd

class SampleMetaData:
    def __init__(self, yaml_file: str | Path):

        with open(yaml_file, "r") as f:
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
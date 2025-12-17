import pandas as pd
from pathlib import Path

REQUIRED_COLUMNS = {"id", "file", "condition", "replicate", "description", "species"}

def load_samples_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix == ".xlsx":
        df = pd.read_excel(path)
    elif path.suffix in {".csv", ".tsv"}:
        seperator = "\t" if path.suffix == ".tsv" else ","
        df = pd.read_csv(path, sep=seperator)
    else:
        raise ValueError(f"Unsupported file extension: {path.suffix}")

    missing_columns = REQUIRED_COLUMNS - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing columns: {missing_columns}")

    return df

# add in functionallity for other metadata columns
def samples_df_to_yaml(df: pd.DataFrame) -> dict:
    samples = []

    for _, row in df.iterrows():
        unique_id = f"{row['id']}_{row['file'].replace(' ', '-').rsplit('__', 1)[1].removesuffix(".mzML")}_{row['replicate']}"


        sample = {
            "unique_id": str(unique_id),
            "id": str(row["id"]).upper(),
            "file": str(row["file"]),
            "condition": str(row["condition"]).capitalize(),
            "replicate": int(row["replicate"]),
            "description": str(row.get("description", "")).capitalize(),
            "species": str(row.get("species")).capitalize(),
        }

        extra_metadata = {
            col: row[col]
            for col in df.columns
            if col not in REQUIRED_COLUMNS
        }

        # Convert NaNs to None for YAML cleanliness
        extra_metadata = {
            k: (None if pd.isna(v) else v)
            for k, v in extra_metadata.items()
        }

        if extra_metadata:
            sample["metadata"] = extra_metadata

        # sample["unique_id"] = f"{sample['id']}_{sample['file'].replace(' ', '-').rsplit('__', 1)[1].removesuffix(".mzML")}_{sample['replicate']}"

        samples.append(sample)

    return {"samples": samples}
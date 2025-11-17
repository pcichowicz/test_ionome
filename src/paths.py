import yaml
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
with open(ROOT / "config/config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

DATA_DIR = ROOT / cfg["data_dir"]
RESULTS_DIR = ROOT / cfg["results_dir"]
LOGS_DIR = ROOT / cfg["logs_dir"]
CACHED_DIR = ROOT / cfg["cached_dir"]


import yaml
from pathlib import Path

# ROOT = Path(__file__).resolve().parents[1]
# with open(ROOT / "projects/test_ionome/config*.yaml", "r") as f:
#     cfg = yaml.safe_load(f)
#
# DATA_DIR = ROOT / cfg["data_dir"]
# RESULTS_DIR = ROOT / cfg["results_dir"]
# LOGS_DIR = ROOT / cfg["logs_dir"]
# CACHED_DIR = ROOT / cfg["cached_dir"]

def output_path(run_id:str, output_dir) -> Path:
    ROOT = Path(__file__).resolve().parents[1]

    # Regex pattern:
    # ^      -> start of the string
    # config -> matches the literal 'config'
    # .*     -> matches any character (.) zero or more times (*)
    # \.     -> matches a literal dot (dots have special meaning in regex, so we escape it)
    # yaml   -> matches the literal 'yaml'
    # $      -> end of the string

    with open(ROOT / "projects" / run_id / f"config_{run_id}.yaml", "r") as f:
        config = yaml.safe_load(f)

    out = ROOT / "projects" / run_id / config[output_dir]

    return out

import numpy as np
import pandas as pd
import scipy.signal
import warnings
import sys
from colorama import Fore, Style, init
init(autoreset=True)
from datetime import datetime

def normalize_signal(intensity: np.ndarray) -> np.ndarray:
    int_sign = np.sign(intensity)
    norm = (intensity - intensity.min()) / (intensity.max() - intensity.min())
    return int_sign * norm

def detect_peak_indices(signal: np.ndarray, prominence: float) -> np.ndarray:
    peaks, _ = scipy.signal.find_peaks(signal, prominence=prominence)
    return peaks

def calculate_peak_widths(intensity: np.ndarray, peak_indices: np.ndarray, rel_height: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=scipy.signal._peak_finding_utils.PeakPropertyWarning)
        widths, _, left_ips, right_ips = scipy.signal.peak_widths(intensity, peak_indices, rel_height=rel_height)
    return widths, left_ips.astype(int), right_ips.astype(int)

def build_peak_ranges(_left, _right, norm_int_len: int, buffer: int) -> list[np.ndarray]:
    ranges = []
    for l, r in zip(_left, _right):
        rnge = np.arange(int(l - buffer), int(r + buffer))
        rnge = rnge[(rnge >= 0) & (rnge < norm_int_len)]
        ranges.append(rnge)
    return ranges

def remove_subset_ranges(ranges: list[np.ndarray]) -> list[np.ndarray]:
    valid = [True] * len(ranges)
    for i, r1 in enumerate(ranges):
        for j, r2 in enumerate(ranges):
            if i != j and set(r2).issubset(r1):
                valid[j] = False
    return [r for i, r in enumerate(ranges) if valid[i]]

def build_window_df(df: pd.DataFrame, ranges: list[np.ndarray]) -> pd.DataFrame:
    df = df.copy()
    df["time_idx"] = np.arange(len(df))
    df["window_id"] = 0
    df["window_type"] = "peak"

    for i, r in enumerate(ranges):
        df.loc[df["time_idx"].isin(r), "window_id"] = i + 1

    return df

def assign_background_windows(window_df: pd.DataFrame) -> pd.DataFrame:
    bg = window_df[window_df["window_id"] == 0]
    tidx = bg["time_idx"].values

    if not len(bg):
        return window_df

    diff = np.diff(tidx)
    split_inds = np.where(diff > 1)[0]

    if len(split_inds) == 0:
        window_df.loc[bg.index, ["window_id", "window_type"]] = [1, "interpeak"]
        return window_df

    split_inds = np.insert(split_inds + 1, 0, 0)
    split_inds = np.append(split_inds, len(tidx))

    for i, (start, end) in enumerate(zip(split_inds[:-1], split_inds[1:])):
        segment = tidx[start:end]
        if len(segment) >= 10:
            window_df.loc[window_df["time_idx"].isin(segment), "window_id"] = i + 1
            window_df.loc[window_df["time_idx"].isin(segment), "window_type"] = "interpeak"

    return window_df[window_df["window_id"] > 0]

def extract_window_props(window_df: pd.DataFrame,
                         peak_indice,
                         time_col,
                         signal_col,
                         timestep_precision,
                         timestep,
                         widths) -> dict:
    window_dict = {}
    for gid, group in window_df[window_df["window_type"] == "peak"].groupby("window_id"):
        peak_idxs = [i for i in peak_indice if i in group["time_idx"].values]
        peak_inds = [np.where(peak_indice == i)[0][0] for i in peak_idxs]
        window_dict[gid] = {
            "time_range": group[time_col].values,
            "signal": group[signal_col].values,
            "signal_area": group[signal_col].sum(),
            "num_peaks": len(peak_idxs),
            "amplitude": [group[group["time_idx"] == p][signal_col].iloc[0] for p in peak_idxs],
            "location": [np.round(group[group["time_idx"] == p][time_col].iloc[0], timestep_precision) for p in peak_idxs],
            "width": [widths[i] * timestep for i in peak_inds],
        }
    return window_dict

def generate_time_range(df, time_col, integration_window, timestep):
    if not integration_window:
        return df[time_col].values
    if len(integration_window) == 2:
        return np.arange(integration_window[0], integration_window[1], timestep)
    raise RuntimeError("Integration window must be empty or [start, stop].")

def default_param_bounds(amplitude,time_min, time_max):
    return {
        "amplitude": np.sort([0.01 * amplitude, 100 * amplitude]),
        "location": [time_min, time_max],
        "scale": [0, (time_max - time_min) / 2],
        "skew": [-np.inf, np.inf],
    }

def sum_skewnorms(x, *params):
    """
    Sum of skew-normal distributions for curve fitting.
    Each peak is represented by 4 parameters: amplitude, center, width, skew.
    """
    from scipy.stats import skewnorm
    n_params_per_peak = 4
    n_peaks = len(params) // n_params_per_peak
    y = np.zeros_like(x)

    for i in range(n_peaks):
        a, loc, scale, skew = params[i * 4:(i + 1) * 4]
        y += skewnorm.pdf(x, skew, loc, scale) * a

    return y

def compute_skewnorm(x, amplitude, loc, scale, alpha):
    _x = alpha * (x - loc) / scale
    norm = (1 / np.sqrt(2 * np.pi * scale**2)) * np.exp(-((x - loc) ** 2) / (2 * scale**2))
    cdf = 0.5 * (1 + scipy.special.erf(_x / np.sqrt(2)))
    return amplitude * 2 * norm * cdf


##########


# def log_method_entry(color = Fore.GREEN):
#     name = sys._getframe(1).f_code.co_name  # 1 = caller
#     timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     print(f"{color}>[{timestamp}][{name}]{Style.RESET_ALL}", end="\n")

def log_method_entry(color=Fore.GREEN):
    frame = sys._getframe(1)  # caller frame
    method_name = frame.f_code.co_name

    # Try to get class name if 'self' is in scope
    class_name = ""
    if "self" in frame.f_locals:
        class_name = frame.f_locals["self"].__class__.__name__

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    if class_name:
        print(f"{color}>[{timestamp}][{class_name}.{method_name}]{Style.RESET_ALL}")
    else:
        print(f"{color}>[{timestamp}][{method_name}]{Style.RESET_ALL}")
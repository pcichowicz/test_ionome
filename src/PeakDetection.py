import tqdm
from time import sleep
from numpy.ma import reshape
from scipy.signal import find_peaks, peak_widths

from helpers import *



class PeakDetection:
    def __init__(self,df_data, **kwargs):
        self.detect = kwargs
        # self.rt_array = rt_array
        # self.i_array = i_array
        self.dataframe = df_data
        self._timestep = float(np.mean(np.diff(self.dataframe["retention_time"])))
        self._timestep_precision = int(np.abs(np.ceil(np.log10(self._timestep))))

        self.new_df = None

    def detect_peaks(self):

        assign_wind = self._assign_windows()

        peak_properties, unmixed_chromatogram = self._deconvolve_peaks(self.dataframe, assign_wind)

        return assign_wind, peak_properties, unmixed_chromatogram

    def _assign_windows(self):

        prominence = self.detect['prominence']
        rel_height = self.detect['rel_height']
        buffer = self.detect['buffer']

        if not (0 <= rel_height <= 1):
            raise ValueError("rel_height must be between 0 and 1")

        normalized_intensity = normalize_signal(self.dataframe["corrected"])
        peak_indices = detect_peak_indices(normalized_intensity,prominence)
        # print(normalized_intensity)
        # print(peak_indices)

        amps = np.sign(self.dataframe["corrected"][peak_indices])
        _widths = np.zeros_like(amps, dtype=int)
        _left = np.zeros_like(amps, dtype=int)
        _right = np.zeros_like(amps, dtype=int)

        if peak_indices.size == 0:
            raise RuntimeError("No peaks detected, cannot assign windows.")

        for sign_mask, polarity in zip([amps > 0, amps < 0], [1,-1]):
            if np.any(sign_mask):
                sig = np.array(self.dataframe["corrected"] if polarity == 1 else - self.dataframe["corrected"])
                w_half, _, _ = calculate_peak_widths(sig, peak_indices, rel_height=0.5)
                _, l_ips, r_ips = calculate_peak_widths(sig, peak_indices, rel_height=rel_height)

                _widths[sign_mask] = w_half[sign_mask]
                _left[sign_mask] = l_ips[sign_mask]
                _right[sign_mask] = r_ips[sign_mask]

                # print(_widths, _left, _right)

        ranges = build_peak_ranges(_left, _right, len(normalized_intensity), buffer)
        ranges = remove_subset_ranges(ranges)
        # print(ranges)
        # df_vertical = pd.DataFrame(np.vstack((self.dataframe["corrected"], self.dataframe["corrected"])), columns=['retention_time', 'intensity'])

        window_df = build_window_df(self.dataframe, ranges)
        # print(window_df)
        window_df = assign_background_windows(window_df)
        # print(window_df)
        window_props = extract_window_props(window_df, peak_indices, "retention_time", "corrected",
                                            self._timestep_precision, self._timestep, _widths)
        self.new_df = window_props
        # return window_df, window_props
        return window_props

    def _deconvolve_peaks(self,
                          data_frame_df,
                          windows_props,
                          integration_window=None,
                          precision = 9,
                          max_iter=5000):

        if integration_window is None:
            integration_window = []
        if windows_props is None:
            raise RuntimeError("Run _assign_windows() first.")

        iterator = tqdm.tqdm(windows_props.items(),
                             desc="Deconvolving peaks",
                             ncols = 100,
                             leave=True,
                             )
        params_order = ["amplitude", "location", "scale", "skew"]
        t_range = generate_time_range(data_frame_df, "retention_time", integration_window, self._timestep)

        peak_props = {}
        self._param_bounds = []
        self._p0 = []

        for k, v in iterator:
            iterator.set_description(f"\tDeconvolving window {k}/{len(iterator)} with {v["num_peaks"]} peaks ")
            if v["num_peaks"] == 0:
                continue

            window_dict = {}
            p0 = []
            bounds_lower, bounds_upper = [], []

            if v["num_peaks"] >= 10:
                warnings.warn(
                    f"\nToo many peaks ({v['num_peaks']}) are detected between {v['time_range'].min()}-{v['time_range'].max()}.\n"
                    f"This may take some time to finish... "
                )

            for i in range(v["num_peaks"]):
                # Get raw guesses
                amp = max(v["amplitude"][i], 1e-6)  # never allow 0 amplitude
                loc = v["location"][i]
                scale = max(v["width"][i] / 2, self._timestep)  # don't allow 0 scale
                skew = 0  # safe default

                # Clamp location inside window
                t_min, t_max = v["time_range"].min(), v["time_range"].max()
                loc = np.clip(loc, t_min, t_max)

                peak_p0 = [amp, loc, scale, skew]

                # Generate bounds
                default_bounds = default_param_bounds(amp, t_min, t_max)
                lower = [default_bounds[k][0] for k in params_order]
                upper = [default_bounds[k][1] for k in params_order]

                # Ensure 0 skew is allowed
                lower[3] = min(lower[3], -5)
                upper[3] = max(upper[3], 5)

                # Clamp p0 to bounds just in case
                peak_p0 = np.clip(peak_p0, lower, upper)

                p0.extend(peak_p0)
                bounds_lower.extend(lower)
                bounds_upper.extend(upper)

                self._p0.append(p0)
                self._param_bounds.append((bounds_lower, bounds_upper))
                # print(f"p0 {self._p0}")
                # print(f"bounds {self._param_bounds}")
            # fit curves
            popt, _ = scipy.optimize.curve_fit(
                sum_skewnorms,
                v["time_range"],
                v["signal"],
                p0=p0,
                bounds=(bounds_lower, bounds_upper),
                maxfev=max_iter
            )
            # print(f"popt\n{popt}")
            popt = reshape(popt, (v["num_peaks"], 4))
            for i, p in enumerate(popt):
                recon_signal = compute_skewnorm(t_range, *p)
                window_dict[f"peak_{i + 1}"] = {
                    "amplitude": p[0],
                    "retention_time": np.round(p[1], decimals=self._timestep_precision),
                    "scale": p[2],
                    "alpha": p[3],
                    "area": recon_signal.sum(),
                    "reconstructed_signal": recon_signal,
                    "signal_max": np.max(recon_signal),
                }
            # print(f"window_dict\n{window_dict}")
            peak_props[k] = window_dict

            iterator.set_description("\tFinished deconvolution ")

        rows = [
            {
                "retention_time": p["retention_time"],
                "scale": p["scale"],
                "skew": p["alpha"],
                "amplitude": p["amplitude"],
                "area": p["area"],
                "signal_maximum": p["signal_max"]
            }
            for window in peak_props.values()
            for p in window.values()
        ]

        if not rows:
            raise ValueError("No peaks were extracted — rows list is empty")


        peak_prop_df = pd.DataFrame(rows).sort_values("retention_time")
        peak_prop_df["peak_id"] = np.arange(1, len(peak_prop_df) +1).astype(int)

        time = self.dataframe["retention_time"]
        out = np.zeros((len(time), len(peak_prop_df)))

        for i, row in peak_prop_df.iterrows():
            params = [row["amplitude"], row["retention_time"], row["scale"], row["skew"]]
            out[:,i] = compute_skewnorm(time, *params)

        unmixed_chromatogram = np.round(out, decimals=precision)

        return peak_props, unmixed_chromatogram

    def detect_xic_peaks(self, xic_df, height_factor=3,prominence_factor=0.01):
        """
            Detect main peak in a single XIC chromatogram DataFrame.
            Expects columns: 'retention_time', 'intensity' (baseline-corrected already).
            """

        rt = xic_df["retention_time"].values
        y = xic_df["intensity"].values

        # --- 1. noise threshold based on median ---
        noise = np.median(y)
        min_height = height_factor * noise

        # --- 2. find peaks ---
        peaks, props = find_peaks(
            y,
            height=min_height,
            prominence=prominence_factor * max(y)
        )

        if len(peaks) == 0:
            return None  # no peak found

        # --- 3. choose most intense peak (targeted analysis assumption) ---
        main = peaks[np.argmax(props["peak_heights"])]

        # --- 4. compute widths ---
        widths, width_heights, left, right = peak_widths(y, [main], rel_height=0.5)

        # --- 5. compute area under peak ---
        left_i = int(left[0])
        right_i = int(right[0])

        # Enforce a minimum width of 3 points
        if right_i - left_i < 3:
            left_i = max(0, main - 2)
            right_i = min(len(y) - 1, main + 2)

        peak_area = np.trapezoid(y[left_i:right_i + 1], rt[left_i:right_i + 1])
        # left_i, right_i = int(left[0]), int(right[0])
        # peak_area = np.trapz(y[left_i:right_i], rt[left_i:right_i])

        # --- 6. return peak dictionary ---
        return {
            "rt": rt[main],
            "height": y[main],
            "area": peak_area,
            "width_rt": widths[0] * np.median(np.diff(rt)),  # approx width in RT units
            "left_rt": rt[left_i],
            "right_rt": rt[right_i],
        }
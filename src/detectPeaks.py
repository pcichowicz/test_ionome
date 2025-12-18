from numpy.ma.core import reshape
from tqdm import tqdm

# from src.sampleData import SampleData
from src.helpers import *

class DetectPeaks:
    def __init__(self,metabolite:str, sample_data: pd.DataFrame, **kwargs):
        self.metabolite = metabolite
        self.sample_data = sample_data
        self.detect_peak_params = kwargs
        self._timestep = float(np.mean(np.diff(self.sample_data["retention_time"])))
        self._timestep_precision = int(np.abs(np.ceil(np.log10(self._timestep))))

    def detect_peaks(self): # runs all
        log_method_entry()

        window_df_properties = self._assign_windows()
        # print("window_df_properties", window_df_properties)
        peak_properties, unmixed_chromatogram = self._deconvolve_peaks(self.sample_data, window_df_properties)
        return window_df_properties, peak_properties, unmixed_chromatogram
        # return window_df_properties

    def _assign_windows(self):
        """r
        Extracts and assigns sections of chromatogram that has peaks into "windows"
        """
        # print(self.metabolite)
        # window_df_properies = None

        # 1. normalize corrected intensity
        normalized_intensities = normalize_signal(self.sample_data["corrected"])
        # print(f"\t 1. normalized_intensities")
        # print(normalized_intensities)

        # 2. Auto-detect peaks using 'scipy.signal.find_peaks'
        peak_indices = detect_peak_indices(normalized_intensities, self.detect_peak_params['prominence'])
        # print(f"\t 2. Auto-detected peaks (index locations)")
        # print(peak_indices)

        # If no peaks, just add in a blank df???
        if peak_indices.size == 0:
            print("No peaks detected, skipping ...")
            # window_df_properies = {
            #     "time_range": None,
            #     "signal": None,
            #     "signal_area": None,
            #     "num_peaks": None,
            #     "amplitude": None,
            #     "location": None,
            #     "width": None
            # }
            return None

        # 3. Determine peak boundaries,'width'
        widths,l_ips, r_ips = calculate_peak_widths(normalized_intensities,peak_indices,self.detect_peak_params['rel_height'])

        left_ips = np.clip(l_ips, 0, len(normalized_intensities)-1)
        right_ips = np.clip(r_ips, 0, len(normalized_intensities)-1)

        # print(f"\t 3. peak width boundaries (index locations)")
        # print(widths,l_ips, r_ips, sep='\n')
        # print(left_ips, right_ips, sep='\n')

        # 4. Build peak ranges
        peak_ranges = build_peak_ranges(left_ips,right_ips) # remember to refactor the #comment lines (remove them if this works)
        # print(f"\t 4. peak ranges (index locations)")
        # print(peak_ranges)


        # 5. Remove subset ranges
        # print(f"\t 5. Subset ranges")
        peak_ranges = remove_subset_ranges(peak_ranges)
        # print(peak_ranges)

        # 6. Build window dataframe
        window_df = self.sample_data.copy()

        built_window = build_window_df(window_df,peak_ranges)
        # print(f"\t 6. built_window")
        # print(built_window)

        # 7. Assign background windows
        # print(f"\t 7. Assign background windows")
        built_window = assign_background_windows(built_window)
        # print(built_window)

        # 8. Extract window
        # print(f"\t 8. Built window properties")
        window_df_properties = extract_window_props(built_window,
                                                    peak_indices,
                                                    "retention_time",
                                                    "corrected",
                                                    self._timestep_precision,
                                                    self._timestep,
                                                    widths)
        # print(window_df_properties)

        return window_df_properties

    def _deconvolve_peaks(self,
                          dataframe_df: pd.DataFrame = None,
                          window_df_properties: pd.DataFrame = None,
                          integration_window: tuple | None = None, # add later
                          precision:int = 9,
                          max_niter:int = 200000):

        # print(dataframe_df)
        # print(window_df_properties)

        if integration_window is None:
            integration_window = []

        if window_df_properties is None:
            print(f"No windows to deconvolve — skipping metabolite.")
            return None, None

        iterator = tqdm(window_df_properties.items(),
                        desc="Deconvoling peaks",
                        ncols=100,
                        leave=True)

        parameter_order = ["amplitude","location", "scale","skew"]
        time_range = generate_time_range(dataframe_df,"retention_time", integration_window, self._timestep)

        peak_props = {}
        self._param_bounds = []
        self._p0 = []

        for k, v in iterator:
            iterator.set_description(f"\t Deconvolving window {k}/{len(iterator)} with {v['num_peaks']}peaks ")
            if v['num_peaks'] == 0:
                break

            window_dict = {}
            p0= []
            bounds_lower, bounds_upper = [], []

            if v['num_peaks'] >= 10:
                warnings.warn("Alot of peaks found")

            for i in range(v['num_peaks']):
                # raw guess
                amp = max(v['amplitude'][i], 1e-6)
                loc = v['location'][i]
                scale = max(v['width'][i] / 2, self._timestep)
                skew = 0

                # clamp location within window range
                t_min, t_max = v["time_range"].min(), v["time_range"].max()
                loc = np.clip(loc, t_min, t_max)

                peak_p0 = [amp, loc, scale, skew]

                # get bounds
                default_bounds = default_param_bounds(amp,t_min,t_max)
                lower = [default_bounds[k][0] for k in parameter_order]
                upper = [default_bounds[k][1] for k in parameter_order]

                #ensure skew 0 is allowed
                lower[3] = min(lower[3], -5)
                upper[3] = max(upper[3], 5)

                peak_p0 = np.clip(peak_p0, lower, upper)

                p0.extend(peak_p0)
                bounds_lower.extend(lower)
                bounds_upper.extend(upper)

                self._p0.append(p0)
                self._param_bounds.append((bounds_lower, bounds_upper))

            #fit curves

            popt, _ = scipy.optimize.curve_fit(sum_skewnorms,
                                               v['time_range'],
                                               v['signal'],
                                               p0=p0,
                                               bounds=(bounds_lower, bounds_upper),
                                               maxfev=max_niter)

            popt = reshape(popt, (v['num_peaks'], 4))
            for i, p in enumerate(popt):
                reconstructed_signal = compute_skewnorm(time_range, *p)
                window_dict[f"peak_{i+1}"] = {
                    "amplitude": p[0],
                    "retention_time": np.round(p[1], decimals=self._timestep_precision),
                    "scale": p[2],
                    "alpha": p[3],
                    "area": reconstructed_signal.sum(),
                    "reconstructed_signal": reconstructed_signal,
                    "signal_max": np.max(reconstructed_signal),

                }

            peak_props[k] = window_dict

            iterator.set_description("\t Finished deconvolution")

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
            print("No peaks were extracted - rows list is empty")

        peak_prop_df = pd.DataFrame(rows).sort_values("retention_time")
        peak_prop_df["peak_id"] = np.arange(1,len(peak_prop_df) + 1).astype(int)

        time = self.sample_data["retention_time"]
        out = np.zeros((len(time), len(peak_prop_df)))

        for i, row in peak_prop_df.iterrows():
            params = [row["amplitude"], row["retention_time"], row["scale"], row["skew"]]
            out[:, i] = compute_skewnorm(time, *params)

        unmixed_chromatogram = np.round(out, decimals=precision)

        return peak_props, unmixed_chromatogram
"""
Class method to perform baseline correction on Chromatorgam data,
both Total Ion Chromatograms (TIC) and Extracted Ion Chromatograms (XIC)

Two implemented methods, SNIP and AsLS

Perform Asymmetric Least Squares (AsLS) baseline correction.

    This method estimates a smooth baseline `b` for a given signal `y`
    by solving the following minimization problem:

        minimize_b  Σ_i [ w_i * (y_i - b_i)^2 ]  +  λ * Σ_i [ (b''_i)^2 ]

    where:
        - w_i : adaptive weights (smaller for peaks above the baseline)
        - λ (lam): smoothness parameter (controls curvature of the baseline)
        - b''_i : second derivative of baseline, enforcing smoothness

    Reference:
    P. H. C. Eilers & H. F. M. Boelens,
    "Baseline Correction with Asymmetric Least Squares Smoothing",
    Leiden University Medical Centre Report, 2005.
"""
import time
import warnings

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import tqdm
from src.helpers import log_method_entry

class BaselineCorrection:
    def __init__(self,**kwargs):
        self._settings = kwargs


    def asls(self,y,**kwargs):
        # for k,v in kwargs.items():
        #     print(k,v)
        """
        Asymmetric Least Squares (AsLS) baseline correction.

        Parameters
        ----------
        y : array-like
            1D input array of intensities (spectrum or chromatogram).
        lam : float
            Smoothness parameter (lambda). Larger = smoother baseline. 1e4-1e6?
        p : float
            Asymmetry parameter (0 < p < 1). Smaller = baseline forced under peaks.
        niter : int
            Maximum number of iterations.
        tol : float
            Convergence tolerance (relative change in baseline).

        Returns
        -------
        baseline : ndarray
            Estimated baseline vector (same length as y).
        corrected : ndarray
            Baseline-corrected signal (y - baseline, clipped at 0).
        """
        lam = float(kwargs.get('lam', 1e5))
        p = float(kwargs.get('p', 0.01))
        niter = int(kwargs.get('niter', 10))
        tol = float(kwargs.get('tol', 1e-6))

        #Check for numpy array type
        if not isinstance(y, np.ndarray):
            y = np.asarray(y)

        L = len(y)  # number of points in signal

        # Construct second-order difference matrix (D), this enforces "smoothness" of the baseline.
        D = sparse.diags([1, -2, 1], [0, 1, 2], shape=(L - 2, L))
        # D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
        DTD = D.T @ D  # Precompute D^T D (penalty term)

        # Initialize weights (all ones at start, size of y numpy array)
        weights = np.ones(L)
        baseline_asls = np.zeros(L)

        # W = sparse.diags(weights, 0, format="csr")
        # A = (W + lam * DTD).tocsr()
        # b = spsolve(A, B)

        # Iterative fitting
        for i in range(niter):
            # Build diagonal weight matrix
            # W = sparse.diags(weights, 0)
            W = sparse.diags(weights, 0, format="csr")
            # Build system: (W + λ D^T * D), b = W * y
            # A = W + lam * DTD
            A = (W + lam * DTD).tocsr()
            B = W @ y

            # Solve linear system for baseline
            b = spsolve(A, B)

            # Check convergence
            if np.linalg.norm(b - baseline_asls) / (np.linalg.norm(baseline_asls) + 1e-8) < tol:
                baseline_asls = b
                break

            baseline_asls = b

            # Update weights based on residuals
            # Points above baseline (peaks) get weight p (small)
            # Points below baseline get weight (1-p) (large)
            residuals = y - baseline_asls
            weights = p * (residuals > 0) + (1 - p) * (residuals <= 0)

        # Step 3. Subtract baseline
        # corrected_asls = y - baseline_asls
        # corrected_asls[corrected_asls < 0] = 0  # clip negatives to zero

        corrected_asls = y - baseline_asls
        corrected_asls = np.maximum(corrected_asls, 0)
        corrected_asls[y == 0] = 0  # Do not invent signal where none existed

        return baseline_asls, corrected_asls

    def snip(self,raw_df, **kwargs):
        # for k,v in kwargs.items():
            # print(k,v)

        window = kwargs.get('window', 5)
        precision = kwargs.get('precision', 9)
        clip_negatives = kwargs.get('clip_negatives', True)

        # # Arrays
        # retention_time_arr = rt_arr
        # intensity_arr = i_arr
        # Arrays
        retention_time_arr = raw_df['retention_time']
        intensity_arr = raw_df['intensity']

        _timestep = float(np.mean(np.diff(retention_time_arr)))
        # _timestep_precision = int(np.abs(np.ceil(np.log10(_timestep))))


        if (window / _timestep) < 10:
            raise ValueError(
                f"The approximate peak width {window} is too small relative to the time sampling interval ({_timestep})."
                f"Either increase the width or set correct_baseline=False to skip this step."
            )

        min_val, max_val = np.min(intensity_arr), np.max(intensity_arr)

        #Warn for significant negative values.
        if min_val < 0 and (min_val.abs() / max_val.abs()) >= 0.1:
            warnings.warn(
                "The chromatogram has significant negative values. Check results visually to determine if"
                "baseline correction was applied correctly."
            )


        # Sampling step
        retention_timestep = float(np.median(np.diff(retention_time_arr)))

        if retention_timestep <= 0:
            raise ValueError("retention_time must be strictly increasing")

        #Minimum required points within window
        min_points = int(((window / retention_timestep) - 1) // 2)
        if min_points < 1:
            raise ValueError("Window too small for SNIP iterations")

        #Shift to avoid log of negative numbers, shift is returned after correction

        shift = np.abs(np.min(intensity_arr)) + 1 if np.min(intensity_arr) < 0 else 0
        y_shifted = intensity_arr + shift

        # Apply LLS transformation
        # lls_transform = np.log(y_shifted + 1)

        lls_transform = np.log(np.log(np.sqrt(y_shifted + 1) + 1) + 1)
        n_iter = min(min_points,25)
        iterator = tqdm.tqdm(range(n_iter),
                             desc="\tPerforming SNIP baseline correction",
                             ncols=100,
                             leave=True
                             )
        for i in iterator:
            time.sleep(0.01)
            left = np.roll(lls_transform, i)
            right = np.roll(lls_transform, -i)
            avg = (left + right) * 0.5

            avg[:i] = lls_transform[:i]
            avg[-i:] = lls_transform[-i:]
            lls_transform = np.minimum(lls_transform, avg)

        iterator.set_description(f"\tFinished SNIP baseline correction")
        time.sleep(0.01)
        # Inverse baseline
        baseline_snip = (np.exp(np.exp(lls_transform) - 1) - 1) ** 2 - 1

        baseline_snip = np.round(baseline_snip - shift, decimals=precision)

        #Corrected baseline signal
        corrected_snip = intensity_arr - baseline_snip
        if clip_negatives:
            corrected_snip = np.maximum(corrected_snip, 0)

        return baseline_snip,corrected_snip


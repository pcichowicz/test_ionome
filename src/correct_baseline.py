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
import pandas as pd
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
        """r
        Implements an Asymmetric Least Squares Smoothing
        baseline correction algorithm (P. Eilers, H. Boelens 2005)

        Baseline Correction with Asymmetric Least Squares Smoothing
        based on https://web.archive.org/web/20200914144852/https://github.com/vicngtor/BaySpecPlots

        Baseline Correction with Asymmetric Least Squares Smoothing
        Paul H. C. Eilers and Hans F.M. Boelens
        October 21, 2005

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
        lam = float(kwargs.get('lam', 1e6))
        p = float(kwargs.get('p', 0.1))
        niter = int(kwargs.get('niter', 10))
        tol = float(kwargs.get('tol', 1e-6))

        #Check for numpy array type
        if not isinstance(y, np.ndarray):
            y = np.asarray(y)

        L = len(y)  # number of points in signal
        w = np.ones(L)
        D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2),format='csr')
        DTD = D @ D.T

        for i in range(niter):
            W = sparse.spdiags(w,0,L,L)
            Z = W + lam * DTD
            z = spsolve(Z,W @ y)
            w = p * (y > z) + (1-p) * (y < z)




        # Construct second-order difference matrix (D), this enforces "smoothness" of the baseline.
        # D = sparse.diags([1, -2, 1], [0, 1, 2], shape=(L - 2, L))

        # DTD = D.T @ D  # Precompute D^T D (penalty term)
        #
        # # Initialize weights (all ones at start, size of y numpy array)
        # weights = np.ones(L)
        # baseline_asls = np.zeros(L)
        #
        #
        # # Iterative fitting
        # for i in range(niter):
        #     # Build diagonal weight matrix
        #
        #     W = sparse.diags(weights, 0,shape=(L,L), format="csr")
        #
        #     A = (W + lam * DTD).tocsr()
        #     B = W @ y
        #
        #     # Solve linear system for baseline
        #     b = spsolve(A, B)


        corrected_asls = y - z
        corrected_asls = np.maximum(corrected_asls, 0)
        corrected_asls[y == 0] = 0  # Do not invent signal where none existed

        return z, corrected_asls

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


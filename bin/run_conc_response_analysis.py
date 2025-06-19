#!/usr/bin/env python

import os
import argparse
import sys
from pathlib import Path
from typing import Dict, Any, Union, List, Optional, Tuple, Callable
import numpy as np
import pandas as pd
import cmdstanpy
import pickle
import time
import json
from contextlib import contextmanager

from shutil import rmtree
from multiprocessing import Pool
from io import StringIO

from scipy.stats import gmean, beta as beta_dist
from scipy.special import logit, expit, betainc, digamma, polygamma, beta
from scipy.optimize import brentq


@contextmanager
def suppress_stdout_stderr() -> None:
    """
    A context manager for suppressing stdout and stderr in Python.

    This will suppress all print statements, even if they originate in compiled
    C/Fortran sub-functions. It will not suppress raised exceptions.
    """
    # Open a pair of null files
    null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
    # Save the actual stdout (1) and stderr (2) file descriptors
    save_fds = (os.dup(1), os.dup(2))

    try:
        # Assign the null pointers to stdout and stderr
        os.dup2(null_fds[0], 1)
        os.dup2(null_fds[1], 2)
        yield
    finally:
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(save_fds[0], 1)
        os.dup2(save_fds[1], 2)
        # Close the null files
        os.close(null_fds[0])
        os.close(null_fds[1])


class BetaLogistic:
    """
    A class representing a double skew logistic distribution.

    This class provides methods for calculating the PDF, CDF, and quantiles
    of a double skew logistic distribution.
    """

    def __init__(self, mu: float, sigma: float, a: float, b: float) -> None:
        """
        Initialize a BetaLogistic instance.

        Args:
            mu: Mean of the distribution
            sigma: Standard deviation of the distribution
            a: Shape parameter for the left tail
            b: Shape parameter for the right tail
        """
        self.mu, self.sigma, self.a, self.b = mu, sigma, a, b
        self.m = digamma(a) - digamma(b)
        self.s = np.sqrt(polygamma(1, a) + polygamma(1, b))

    def cdf(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate the cumulative distribution function (CDF).

        Args:
            x: Value(s) at which to evaluate the CDF

        Returns:
            CDF evaluated at x
        """
        y = self.m + self.s * (x - self.mu) / self.sigma

        if isinstance(y, (list, np.ndarray)):
            cdf = np.zeros(len(y))
            index = np.where(y <= 0)[0]
            cdf[index] = betainc(self.a, self.b, expit(y[index]))
            index = np.where(y > 0)[0]
            cdf[index] = 1 - betainc(self.b, self.a, expit(-y[index]))
        else:
            cdf = betainc(self.a, self.b, expit(y)) if y <= 0 else 1 - betainc(self.b, self.a, expit(-y))

        return cdf

    def ppf(self, q: float) -> float:
        """
        Calculate the percent point function (inverse of CDF).

        Args:
            q: Quantile at which to evaluate the PPF

        Returns:
            Value x such that P(X <= x) = q
        """
        def func(x: float, *args: Any) -> float:
            return self.cdf(x) - args[0]

        # Find a value of x with func < 0
        lower_bracket = -5
        while True:
            if func(lower_bracket, *(q,)) < 0:
                break
            else:
                lower_bracket -= 5

        # Find a value of x with func > 0
        upper_bracket = 5
        while True:
            if func(upper_bracket, *(q,)) > 0:
                break
            else:
                upper_bracket += 5

        x = brentq(func, a=lower_bracket, b=upper_bracket, args=(q, ))
        return x

    def pdf(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate the probability density function (PDF).

        Args:
            x: Value(s) at which to evaluate the PDF

        Returns:
            PDF evaluated at x
        """
        y = self.m + self.s * (x - self.mu) / self.sigma
        pdf = expit(y) ** self.a * expit(-y) ** self.b / beta(self.a, self.b) * self.s / self.sigma
        return pdf

    def logpdf(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Calculate the log probability density function.

        Args:
            x: Value(s) at which to evaluate the log PDF

        Returns:
            Log PDF evaluated at x
        """
        y = self.m + self.s * (x - self.mu) / self.sigma
        logpdf = (- self.a * np.logaddexp(0, -y) - self.b * np.logaddexp(0, y)
                  - np.log(beta(self.a, self.b)) + np.log(self.s) - np.log(self.sigma))
        return logpdf


def _categorize_by_thresholds(value: float, thresholds: List[float], labels: List[str]) -> str:
    """Helper function to categorize a numeric value using thresholds."""
    bin_index = np.digitize(value, thresholds)
    return labels[min(bin_index, len(labels) - 1)]


def generate_results_summary(
    data: Dict[str, Any],
    samples: Dict[str, Any],
    fit_results: pd.Series,
    diagnostics: str
) -> Dict[str, Any]:
    """
    Generate a JSON-serializable summary of key scientific conclusions.

    This creates consistent categorical classifications for BIFROST results,
    using CDS > 0.5 as the primary threshold for biological significance.
    Designed to be deterministic for reproducible snapshots.

    Args:
        data: Input data dictionary
        samples: MCMC samples dictionary (unused but kept for API consistency)
        fit_results: Fitted response curve results
        diagnostics: Diagnostic string from Stan

    Returns:
        Dictionary containing categorized scientific conclusions
    """
    # Extract core metrics
    pod_samples = fit_results['pod']
    cds_value = float(fit_results['cds'])
    response_curves = fit_results['response']

    # PoD analysis
    finite_pods = pod_samples[np.isfinite(pod_samples)] if len(pod_samples) > 0 else []
    detection_rate = len(finite_pods) / len(pod_samples) if len(pod_samples) > 0 else 0.0

    # Response analysis
    max_response = float(np.max(response_curves))
    min_response = float(np.min(response_curves))
    max_abs_response = max(abs(max_response), abs(min_response))

    # Threshold analysis
    threshold_lower = float(fit_results['response_threshold_lower'])
    threshold_upper = float(fit_results['response_threshold_upper'])
    exceeds_threshold = (max_response > threshold_upper or min_response < threshold_lower)

    # Categorize CDS (biological significance)
    cds_exceeds_threshold = cds_value > 0.5
    cds_category = _categorize_by_thresholds(
        cds_value,
        [0.1, 0.3, 0.5, 0.8],
        ['very_low', 'low', 'moderate', 'high', 'very_high']
    )

    # Categorize detection rate (statistical power)
    signal_strength = _categorize_by_thresholds(
        detection_rate,
        [0.1, 0.5, 0.7, 0.9],
        ['very_low', 'low', 'moderate', 'high', 'very_high']
    )

    # Effect direction
    if abs(max_response) > abs(min_response):
        effect_direction = 'up'
    elif abs(min_response) > abs(max_response):
        effect_direction = 'down'
    else:
        effect_direction = 'none'

    # Effect magnitude (with CDS-informed adjustment)
    base_magnitude = _categorize_by_thresholds(
        max_abs_response,
        [10, 50, 100],
        ['minimal', 'small', 'moderate', 'large']
    )

    # Adjust magnitude based on CDS confidence
    if cds_value < 0.3 and base_magnitude in ['moderate', 'large']:
        effect_magnitude = 'small'
    elif cds_value < 0.1 and base_magnitude != 'minimal':
        effect_magnitude = 'minimal'
    else:
        effect_magnitude = base_magnitude

    # Convergence assessment
    convergence_issues = []
    if 'R-hat greater than 1.01' in diagnostics:
        convergence_issues.append('high_rhat')
    if 'divergent transitions found' in diagnostics and 'No divergent transitions found' not in diagnostics:
        convergence_issues.append('divergent_transitions')
    if 'Treedepth satisfactory' not in diagnostics:
        convergence_issues.append('max_treedepth')
    if 'E-BFMI satisfactory' not in diagnostics:
        convergence_issues.append('low_bfmi')
    if 'effective sample size satisfactory' not in diagnostics:
        convergence_issues.append('low_ess')

    if len(convergence_issues) == 0:
        convergence_status = 'excellent'
    elif len(convergence_issues) == 1 and 'high_rhat' in convergence_issues:
        convergence_status = 'acceptable'
    elif len(convergence_issues) <= 2:
        convergence_status = 'issues_detected'
    else:
        convergence_status = 'poor'

    # Overall interpretation confidence
    if cds_value > 0.8 and convergence_status in ['excellent', 'acceptable']:
        interpretation_confidence = 'high'
    elif cds_value > 0.5 and convergence_status != 'poor':
        interpretation_confidence = 'moderate'
    else:
        interpretation_confidence = _categorize_by_thresholds(
            cds_value,
            [0.3],
            ['very_low', 'low']
        )

    return {
        'scientific_conclusions': {
            'cds_exceeds_threshold': bool(cds_exceeds_threshold),
            'effect_direction': str(effect_direction),
            'effect_magnitude': str(effect_magnitude),
            'effect_exceeds_threshold': bool(exceeds_threshold),
            'cds_category': str(cds_category),
            'interpretation_confidence': str(interpretation_confidence)
        },
        'pod_analysis': {
            'detection_rate': float(round(detection_rate, 6)),
            'statistical_signal_strength': str(signal_strength)
        },
        'response_analysis': {
            'response_shape': list(response_curves.shape) if hasattr(response_curves, 'shape') else None,
            'max_absolute_response': float(round(max_abs_response, 2))
        },
        'model_quality': {
            'convergence_status': str(convergence_status),
            'convergence_issues': list(convergence_issues) if convergence_issues else [],
            'has_regularization_warning': bool('regularizating your model' in diagnostics)
        },
        'sample_metadata': {
            'n_concentrations': int(data['n_conc']),
            'n_batches': int(data['n_batch']),
            'n_samples_total': int(data['n_sample']),
            'mcmc_samples': int(data['n_samp'])
        }
    }


def get_inits(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate initial values for the model parameters.

    Args:
        data: Stan input dictionary

    Returns:
        Dictionary of initial values for model parameters
    """
    log_odds = logit((np.array(data['count']) + 0.5) / (np.array(data['total_count']) + 1))

    mu = np.empty(data['n_batch'])
    for i, idx in enumerate(np.unique(data['batch_index'])):
        mask = np.array(data['batch_index']) == idx
        mu[i] = np.mean(log_odds[mask]) + 10

    return {'log_odds': log_odds, 'mu': mu, 'theta_raw': 0.}


def fit_model(
    path_to_executable: Union[str, Path],
    data: Dict[str, Any],
    seed: Optional[int] = None
) -> Dict[str, Any]:
    """
    Fit the BIFROST model using PyStan.

    Args:
        path_to_executable: Path to compiled Stan model
        data: Data object to be passed to the model
        seed: Optional random seed for reproducibility

    Returns:
        Dictionary containing the posterior samples and diagnostics
    """
    # Attempt using standard settings
    model = cmdstanpy.CmdStanModel(exe_file=path_to_executable)
    fit = model.sample(data=data,
                       chains=4,
                       parallel_chains=1,
                       iter_warmup=500,
                       iter_sampling=250,
                       thin=1,
                       inits=get_inits(data),
                       save_warmup=False,
                       max_treedepth=15,
                       adapt_delta=0.95,
                       seed=seed,
                       show_console=True)

    # Extract diagnostics
    diagnostics = fit.diagnose()

    # Check for multimodality and refit with more chains if detected
    s = 'Split R-hat values satisfactory all parameters.'
    if s not in diagnostics:
        fit = model.sample(data=data,
                           chains=40,
                           parallel_chains=1,  # Use a single core since we're parallelising across input files
                           iter_warmup=500,
                           iter_sampling=250,
                           thin=10,
                           inits=get_inits(data),
                           save_warmup=False,
                           max_treedepth=15,
                           adapt_delta=0.95,
                           seed=seed,
                           show_console=True)

        diagnostics = fit.diagnose()

    # Extract samples
    samples = pd.Series(fit.stan_variables())
    pars = fit.draws_pd()
    for i in ['chain__', 'iter__', 'draw__', 'lp__', 'accept_stat__', 'stepsize__',
              'treedepth__', 'n_leapfrog__', 'divergent__', 'energy__']:
        samples[i] = pars[i]

    return {'samples': samples, 'diagnostics': diagnostics}


def calc_pod_sample(
    conc: np.ndarray,
    response: np.ndarray,
    lower_limit: float,
    upper_limit: float
) -> float:
    """
    Calculate the PoD given a sample of the curve describing the mean response.

    Args:
        conc: Array of concentrations at which the curve has been evaluated
        response: Array containing sample for mean response
        lower_limit: Lower limit of the distribution for the control response
        upper_limit: Upper limit of the distribution for the control response

    Returns:
        Sample for the PoD based on the supplied sample of the concentration-response
    """
    # Determine which direction of largest change
    abs_response_up = np.abs(np.max(response))
    abs_response_down = np.abs(np.min(response))

    if abs_response_up > abs_response_down:
        response_direction = 'up'
    else:
        response_direction = 'down'

    pod = np.inf
    if response_direction == 'up' and np.max(response) > upper_limit:
        index = np.argmax(response)
        for i in range(index):
            if response[index - i] < upper_limit:
                pod = (upper_limit * (conc[index - (i - 1)] - conc[index - i]) -
                       (response[index - i] * conc[index - (i - 1)] - response[index - (i - 1)] * conc[index - i])) \
                      / (response[index - (i - 1)] - response[index - i])
                break

    elif response_direction == 'down' and np.min(response) < lower_limit:
        index = np.argmin(response)
        for i in range(index):
            if response[index - i] > lower_limit:
                pod = (lower_limit * (conc[index - (i - 1)] - conc[index - i]) -
                       (response[index - i] * conc[index - (i - 1)] - response[index - (i - 1)] * conc[index - i])) \
                      / (response[index - (i - 1)] - response[index - i])
                break

    return pod


def get_bifrost_covariance(
    data: Dict[str, Any],
    samples: Dict[str, Any],
    conc: Optional[np.ndarray] = None,
    add_sigma: bool = True
) -> np.ndarray:
    """
    Compute the BIFROST kernel for the supplied concentration arrays.

    Args:
        data: Dictionary of concentration-response data
        samples: Dictionary of parameter estimates
        conc: Optional array of concentrations to extrapolate to
        add_sigma: Whether to add the sigma term to the covariance

    Returns:
        Covariance matrix
    """
    n_samp = len(samples['lp__'])
    if conc is None:
        Sigma = np.zeros((n_samp, data['n_conc'], data['n_conc']))
        for i in range(data['n_conc']):
            ci = data['conc'][i]

            if add_sigma:
                Sigma[:, i, i] += np.square(samples['sigma']) / data['n_treatment_batch']

            theta = samples['theta']
            beta = samples['beta']
            gamma = samples['gamma']
            Sigma[:, i, i] += gamma ** 2 / (1 + np.exp(np.log(19) * (ci - beta) / (theta - beta))) ** 2

            for j in range(i):
                cj = data['conc'][j]
                rho = samples['rho']

                Sigma[:, i, j] += (gamma ** 2 / (1 + np.exp(np.log(19) * (ci - beta) / (theta - beta)))
                                   / (1 + np.exp(np.log(19) * (cj - beta) / (theta - beta)))
                                   * np.exp(- 0.5 * ((ci - cj) / rho) ** 2))

                # Fill opposite diagonal
                Sigma[:, j, i] = Sigma[:, i, j]

    else:
        n = len(conc)
        Sigma = np.zeros((n_samp, n, data['n_conc']))
        for i in range(n):
            ci = conc[i]
            for j in range(data['n_conc']):
                cj = data['conc'][j]
                theta = samples['theta']
                beta = samples['beta']
                gamma = samples['gamma']
                rho = samples['rho']

                Sigma[:, i, j] += (gamma ** 2 / (1 + np.exp(np.log(19) * (ci - beta) / (theta - beta)))
                                   / (1 + np.exp(np.log(19) * (cj - beta) / (theta - beta)))
                                   * np.exp(- 0.5 * ((ci - cj) / rho) ** 2))

    return Sigma


def run_concentration_response_analysis(
    files_to_process: List[Union[str, Path]],
    model_executable: Union[str, Path],
    number_of_cores: int,
    fit_dir: Optional[Union[str, Path]] = None,
    seed: Optional[int] = None
) -> None:
    """
    Fit Pystan model for dataset specified by chemical and cell type.

    Args:
        files_to_process: List of probe .pkl files to process
        model_executable: Path to the compiled Stan model executable
        number_of_cores: Number of cores to use
        fit_dir: Optional directory to contain model fits
        seed: Optional random seed for reproducibility

    Raises:
        ValueError: If fit_dir is not a string
        FileNotFoundError: If any input file does not exist
    """
    # Define path to directory to contain model fits
    if fit_dir is None:
        path_to_fits = Path('Fits')
    elif isinstance(fit_dir, (str, Path)):
        path_to_fits = Path(fit_dir) / 'Fits'
    else:
        raise ValueError('Directory to contain model fits must be specified as a string or Path')

    # Create directory if it does not exist
    path_to_fits.mkdir(parents=True, exist_ok=True)

    # Check all inputs are present
    for f in files_to_process:
        if not Path(f).is_file():
            raise FileNotFoundError(f"Data file '{f}' does not exist")

    # Create list of arguments to pass to standard_analysis function
    fitting_args = [(str(model_executable),
                     i,
                     path_to_fits / f"{Path(i).stem}.pkl",
                     number_of_cores,
                     seed)
                    for i in files_to_process]

    with Pool(number_of_cores) as p:
        p.map(standard_analysis, fitting_args)


def standard_analysis(paths: Tuple[Union[str, Path], ...]) -> None:
    """
    Wrapper for the functions used to fit model and generate plotting data.

    Args:
        paths: Tuple containing paths to:
            - model executable
            - data file
            - fit file
            - number of cores (unused)
            - optional seed for reproducibility
    """
    path_to_executable, path_to_data, path_to_fit, _, seed = paths

    with open(path_to_data, 'rb') as f:
        data = pickle.load(f)

    # Generate posterior samples
    with suppress_stdout_stderr():
        fit_dict = fit_model(path_to_executable, data, seed)

        # Generate model fits
        gen_plotting_data(data,
                          fit_dict['samples'],
                          path_to_fit,
                          fit_dict['diagnostics'])


def get_response_window(samples: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculate response window for calculating PoDs.

    Args:
        samples: Dictionary of parameter samples

    Returns:
        Updated samples dictionary with response window
    """
    rtl = np.array([BetaLogistic(0, s, a, b).ppf(0.05)
                    for s, a, b in zip(samples['sigma'], samples['a'], samples['b'])])
    rtu = np.array([BetaLogistic(0, s, a, b).ppf(0.95)
                    for s, a, b in zip(samples['sigma'], samples['a'], samples['b'])])

    samples['rtl'] = rtl
    samples['rtu'] = rtu

    return samples


def interpolate_treatment_effect(
    data: Dict[str, Any],
    samples: Dict[str, Any]
) -> pd.Series:
    """
    Calculate the posterior predictive mean effect of the treatment.

    Args:
        data: Dictionary of concentration-response data
        samples: Dictionary of parameter estimates

    Returns:
        Series containing concentration-response fit
    """
    n_x = 100
    x = np.linspace(min(np.min(data['conc']), np.min(samples['theta'])), np.max(data['conc']), n_x)
    Sigma = get_bifrost_covariance(data, samples)
    Sigma_inv = np.array([np.linalg.inv(i) for i in Sigma])
    Sigma_extrapolation = get_bifrost_covariance(data, samples, x)
    treatment_response = np.array([
        St.dot(Sinv).dot(tr)
        for St, Sinv, tr in zip(Sigma_extrapolation, Sigma_inv, samples['treatment_response'])])

    # Calculate PoDs
    rtl, rtu = samples['rtl'], samples['rtu']

    pod = np.array([calc_pod_sample(x, i, u, v)
                    for i, u, v in zip(treatment_response, rtl, rtu)]).astype('float')
    cds = np.sum(~np.isinf(pod)) / len(pod)
    pod = pod[~np.isinf(pod)]

    # Convert treatment response to expected count
    median_total_count = np.median(data['total_count'])
    log_odds = np.array([i + j - 10 for i, j in zip(np.mean(samples['mu'], axis=1), treatment_response)])
    prob = expit(log_odds)
    expected_count = prob * median_total_count
    expected_count_percentiles = np.percentile(expected_count, q=(2.5, 50, 97.5), axis=0)

    results_dict = pd.Series({
        'x': x,
        'response': expected_count_percentiles,
        'response_threshold_lower': gmean(expit(np.mean(samples['mu'], axis=1) + rtl - 10) * median_total_count),
        'response_threshold_upper': gmean(expit(np.mean(samples['mu'], axis=1) + rtu - 10) * median_total_count),
        'pod': pod,
        'cds': cds,
    })

    return results_dict


def gen_plotting_data(
    data: Dict[str, Any],
    samples: Dict[str, Any],
    path_to_output: Union[str, Path],
    diagnostics: str
) -> None:
    """
    Generate dose response curves for the BIFROST model.

    Args:
        data: Data used to estimate model parameters
        samples: Posterior samples from the model fit
        path_to_output: Path to which the plotting data will be stored
        diagnostics: Diagnostic string for the fit
    """
    # Calculate response window and add to samples
    samples = get_response_window(samples)

    # Add expected samples values to data file
    data['n_samp'] = len(samples['lp__'])
    data['max_conc'] = np.max(data['conc'])
    data['parameters'] = {}
    pars = ['mu', 'sigma', 'a', 'b',
            'log_odds',
            'treatment_response', 'theta', 'beta', 'gamma', 'rho',
            'rtl', 'rtu'
            ]

    for p in pars:
        if samples[p].ndim == 1:
            data['parameters'][p] = np.mean(samples[p])
        else:
            data['parameters'][p] = np.mean(samples[p], axis=0)

    # Interpolate treatment effects and add to data file
    data['fit'] = interpolate_treatment_effect(data, samples)

    # Add diagnostics
    data['diagnostics'] = diagnostics

    # Save the full pickle file
    with open(path_to_output, 'wb') as f:
        pickle.dump(data, f)

    # Generate and save JSON summary
    json_summary = generate_results_summary(data, samples, data['fit'], diagnostics)
    json_path = Path(path_to_output).with_suffix('.json')
    with open(json_path, 'w') as f:
        json.dump(json_summary, f, indent=2)


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-files', nargs='+', type=str,
                        help='list of probe .pkl files to process separated by spaces')
    parser.add_argument('--model-executable', type=str, help='path to compiled Stan model executable')
    parser.add_argument('--n-cores', type=int, help='number of cores to use in multiprocessing')
    parser.add_argument('--seed', type=int, help='optional random seed for reproducibility')
    args = parser.parse_args()

    run_concentration_response_analysis(args.data_files,
                                        args.model_executable,
                                        args.n_cores,
                                        seed=args.seed)


if __name__ == '__main__':
    main()

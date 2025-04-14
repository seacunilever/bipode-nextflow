#!/usr/bin/env python

import os
import argparse
import sys
import numpy as np
import pandas as pd
import cmdstanpy
import pickle
import time

from shutil import rmtree
from multiprocessing import Pool
from io import StringIO

from scipy.stats import gmean, beta as beta_dist
from scipy.special import logit, expit, betainc, digamma, polygamma
from scipy.optimize import brentq


class suppress_stdout_stderr(object):
    """
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    """

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for _ in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


class BetaLogistic:

    def __init__(self, mu, sigma, a, b):
        """
        Initalises an instance of the class representing a double skew logistic distribution.

        Accepts:
            mu - mean of distribution
            sigma - standard deviation of distribution
            a - shape parameter for left tail
            b - shape parameter for right tail

        Returns:
             None
        """

        self.mu, self.sigma, self.a, self.b = mu, sigma, a, b
        self.m = digamma(a) - digamma(b)
        self.s = np.sqrt(polygamma(1, a) + polygamma(1, b))

    def cdf(self, x):
        """
        Returns the CDF of double skew logistic distribution.

        Accepts:
            x - independent variable

        Returns:
            cdf - cumulative density function evaluated at x
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

    def ppf(self, q):
        """
        Computes quantiles of the double skew logistic distribution

        Accepts:
            q - quantile

        Returns:
            x
        """

        def func(x, *args):
            return self.cdf(x) - args[0]

        # Find a value of x with func < 0
        lower_bracket = -5
        while True:
            if func(lower_bracket, *(q,)) < 0:
                break
            else:
                lower_bracket -= 5

        # Find a value of x with func < 0
        upper_bracket = 5
        while True:
            if func(upper_bracket, *(q,)) > 0:
                break
            else:
                upper_bracket += 5

        x = brentq(func, a=lower_bracket, b=upper_bracket, args=(q, ))

        return x

    def pdf(self, x):
        """
        Returns the PDF of double skew logistic distribution.

        Accepts:
            x - independent variable

        Returns:
            pdf - probability density function evaluated at x
        """

        y = self.m + self.s * (x - self.mu) / self.sigma
        pdf = expit(y) ** self.a * expit(-y) ** self.b / beta(self.a, self.b) * self.s / self.sigma

        return pdf

    def logpdf(self, x):
        """
        Returns the log PDF of double skew logistic distribution.

        Accepts:
            x - independent variable

        Returns:
            logpdf - log probability density function evaluated at x
        """

        y = self.m + self.s * (x - self.mu) / self.sigma
        logpdf = (- self.a * np.logaddexp(0, -y) - self.b * np.logaddexp(0, y)
                  - np.log(beta(self.a, self.b)) + np.log(self.s) - np.log(self.sigma))

        return logpdf


def compile_stan_model(path_to_stan_file):
    """
    Compiles the Stan model

    Args:
        path_to_stan_file (str) - path to Stan code defining the model

    Returns:
        executable_file_path (str) - path to compiled executable
    """

    # Compile PyStan model and pickle
    model = cmdstanpy.CmdStanModel(stan_file=path_to_stan_file)

    return model.exe_file


def get_inits(data):
    """
    Returns point-estimates of the log-odds for all counts in the
    provided data file.

    Accepts:
        data (dict) - Stan input dictionary

    Returns:
         log_odds (dict) - log-odds estimates
    """

    log_odds = logit((np.array(data['count']) + 0.5) / (np.array(data['total_count']) + 1))

    mu = np.empty(data['n_batch'])
    for i, idx in enumerate(np.unique(data['batch_index'])):
        mask = np.array(data['batch_index']) == idx
        mu[i] = np.mean(log_odds[mask]) + 10

    return {'log_odds': log_odds, 'mu': mu, 'theta_raw': 0.}


def fit_model(path_to_executable, data):
    """
    This function fits the BIFROST model using PyStan.
    The stan_utility library is required to export convergence diagnostics

    Accepts:
        model - instance of the BIFROST model
        data (pd.Series/dict) - data object to be passed to the model

    Returns:
        a dictionary containing the posterior samples and the diagnostic string
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
                       )

    # Extract diagnostics if return_diagnostics is true
    diagnostics = fit.diagnose()

    # Check for mulitmodality and refit with more chains if detected
    s = 'Split R-hat values satisfactory all parameters.'
    if s not in diagnostics:
        fit = model.sample(data=data,
                           chains=40,
                           parallel_chains=n_cores,
                           iter_warmup=500,
                           iter_sampling=250,
                           thin=10,
                           inits=get_inits(data),
                           save_warmup=False,
                           max_treedepth=15,
                           adapt_delta=0.95,
                           )

        diagnostics = fit.diagnose()

    # Extract samples
    samples = pd.Series(fit.stan_variables())
    pars = fit.draws_pd()
    for i in ['chain__', 'iter__', 'draw__', 'lp__', 'accept_stat__', 'stepsize__',
              'treedepth__', 'n_leapfrog__', 'divergent__', 'energy__']:
        samples[i] = pars[i]

    return {'samples': samples, 'diagnostics': diagnostics}


def calc_pod_sample(conc, response, lower_limit, upper_limit):
    """
    Calculates the PoD given a sample of the curve describing the mean response

    Accepts:
        conc - 1D array of concentrations at which the curve has been evaluated
        response - 1D array containing sample for mean response
        lower_limit - lower limit of the distribution for the control response
        upper_limit - upper limit of the distribution for the control response

    Returns:
        pod - sample for the PoD based on the supplied sample of the concentration-response
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


def get_bifrost_covariance(data, samples, conc=None, add_sigma=True):
    """
    Computes the BIFROST kernel for the supplied concentration arrays.

    Arguments:
        data (pd.Series) - series of concentration-response
        samples (collections.OrderedDict) - dictionary of parameter estimates
        conc (np.ndarray) (optional) - concentrations to extrapolate to
        add_sigma (bool) - True/False used to decide whether s term is added

    Returns:
        Sigma (np.ndarray) - Covariance matrix
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

                # Fill oppositr diagonal
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


def run_concentration_response_analysis(files_to_process, model_name, number_of_cores, fit_dir=None):
    """
    Uses the multiprocessing module to fit Pystan model for dataset specified by
    chemical and cell type.

    Accepts:    1. analysis_dir - directory containing data to be processed
                2. Path to executable for Stan model
                3. Number of cores to use

    Returns:    None
    """

    # Define path to directory to contain model fits
    if fit_dir is None:
        path_to_fits = f'Fits'
    elif isinstance(fit_dir, str):
        path_to_fits = f'{fit_dir}/Fits'
    else:
        raise ValueError(f'Directory to contain model fits must be specified as a string')

    # Create directory if it does not exist
    if not os.path.exists(path_to_fits):
        os.makedirs(path_to_fits)

    # Check all inputs are present
    for f in files_to_process:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Data file '{f}' does not exist")

    # Compile model
    path_to_model = compile_stan_model(model_name)

    # Create list of arguments to pass to standard_analysis function
    fitting_args = [(path_to_model,
                     i,
                     f'{path_to_fits}/{os.path.splitext(os.path.basename(i))[0]}.pkl')
                    for i in files_to_process]

    with Pool(number_of_cores) as p:
        p.map(standard_analysis, fitting_args)


def standard_analysis(paths):
    """
    Wrapper for the functions used to fit model and generate plotting data

    Also generates concentration-response curves.

    Argument:
        paths (tuple) - tuple of paths to the model instance, data file, and fit file

    Returns:    None
    """

    path_to_executable, path_to_data, path_to_fit = paths

    data = pickle.load(open(path_to_data, 'rb'))

    # Generate posterior samples
    with suppress_stdout_stderr():
        fit_dict = fit_model(path_to_executable, data)

        # Generate model fits
        gen_plotting_data(data,
                          fit_dict['samples'],
                          path_to_fit,
                          fit_dict['diagnostics'])


def get_response_window(samples):
    """
    Updates samples file with response window for calculating PoDs
    """

    rtl = np.array([BetaLogistic(0, s, a, b).ppf(0.05)
                    for s, a, b in zip(samples['sigma'], samples['a'], samples['b'])])
    rtu = np.array([BetaLogistic(0, s, a, b).ppf(0.95)
                    for s, a, b in zip(samples['sigma'], samples['a'], samples['b'])])

    samples['rtl'] = rtl
    samples['rtu'] = rtu

    return samples


def interpolate_treatment_effect(data, samples):
    """
    Calculates the posterior predictive mean effect of the treatment

    Accepts:
        data (pd.Series) - series of concentration-response
        samples (collections.OrderedDict) - dictionary of parameter estimates

    Returns:
        results_dict (pd.Series) - pandas series containing concentration-response fit
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


def gen_plotting_data(data, samples, path_to_output, diagnostics):
    """
    Generate dose response curves for the BIFROST model

    Accepts:
        data (pd.Series/dict) - data used to estimate model parameters
        samples (collections.OrderedDict) - posterior samples from the model fit
        path_to_output (str) - path to which the plotting data will be stored
        diagnostics (str) - diagnostic string for the fit

    Returns:    None
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

    # Interpolate treatment effects and add to data file, indexed by chemical ID
    data['fit'] = interpolate_treatment_effect(data, samples)

    # Add diagnostics
    data['diagnostics'] = diagnostics

    pickle.dump(data, open(path_to_output, 'wb'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-files', nargs='+', type=str,
                        help='list of probe .pkl files to process separated by spaces')
    parser.add_argument('--model-name', type=str, help='model name')
    parser.add_argument('--n-cores', type=int, help='number of cores to use in multiprocessing')
    args = parser.parse_args()

    run_concentration_response_analysis(args.data_files,
                                        args.model_name,
                                        args.n_cores)
#!/usr/bin/env python3

from __future__ import annotations

import os
import time
import logging
from typing import Dict, List, Optional, Tuple, Union, Any
import numpy as np
import pandas as pd
import argparse
import plotly.graph_objects as go
from pylatex import Document, Section, Subsection, Package, Tabular, LongTable, Math, TikZ, \
    Plot, Figure, SubFigure, Matrix, Alignat, NoEscape, MultiColumn, MiniPage, LargeText, LineBreak, \
    NewLine, Itemize, NewPage, PageStyle, Head, Foot, simple_page_number, StandAloneGraphic, Subsubsection
from pylatex.utils import italic, bold

import plotly.io as pio
pio.kaleido.scope.mathjax = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def calculate_summary_statistics(df: pd.Series) -> Dict[str, Union[np.ndarray, float, int]]:
    """Calculates summary statistics for BIFROST analysis.

    This function computes various summary statistics including PoD mean, log2 fold-change
    extrema, and CDS scores for each probe in the dataset.

    Args:
        df: BIFROST summary data as a pandas Series containing probe information.

    Returns:
        A dictionary containing:
            - probe: Array of probe identifiers
            - pod: Array of PoD means for each probe
            - cds: Array of CDS scores for each probe
            - l2fc: Array of log2 fold changes for each probe
            - max_conc: Maximum tested concentration
            - n_samp: Number of samples
            - conc: Array of concentration values

    Raises:
        KeyError: If required keys are missing from the input Series.
        ValueError: If probe data is malformed.
    """
    logger.info("Starting summary statistics calculation...")

    # Extend PoD sample arrays using maximum concentration value to redefine PoD percentile estimates
    probes = np.array(df['probes'])
    logger.info(f"Processing {len(probes)} probes...")

    pod = np.array([np.mean(df[i]['pod']) if len(df[i]['pod']) > 0 else np.nan for i in df['probes']])
    logger.info(f"Calculated PoD means for {np.sum(~np.isnan(pod))} probes")

    cds = np.array([df[i]['cds'] for i in probes])
    logger.info(f"Calculated CDS scores for {len(cds)} probes")

    l2fc = np.empty(probes.shape[0], dtype='float')
    logger.info("Calculating log2 fold changes...")
    for i, probe in enumerate(probes):
        y = np.array(df[probe]['response'][1])
        index = np.argmax(np.abs(np.log2(y / y[0])))
        l2fc[i] = np.log2(y[index] / y[0])
    logger.info("Log2 fold changes calculation complete")

    stats = {
        'probe': probes,
        'pod': pod,
        'cds': cds,
        'l2fc': l2fc,
        'max_conc': df['max_conc'],
        'n_samp': df['n_samp'],
        'conc': df['conc']
    }

    return stats

def filter_summary_statistics(df: Dict[str, np.ndarray], cds_threshold: float) -> Dict[str, np.ndarray]:
    """Filters summary statistics based on CDS threshold.

    Args:
        df: Dictionary containing summary statistics with keys 'probe', 'pod', 'cds', 'l2fc'.
        cds_threshold: Minimum CDS value to keep in filtered results.

    Returns:
        Filtered dictionary with same structure as input, containing only entries where
        CDS >= cds_threshold.

    Raises:
        KeyError: If required keys are missing from input dictionary.
    """
    mask = df['cds'] >= cds_threshold
    filtered_df = df.copy()
    for key in ['probe', 'pod', 'cds', 'l2fc']:
        filtered_df[key] = df[key][mask]
    return filtered_df

def get_confidence_threshold_probability_density(x: np.ndarray) -> np.ndarray:
    """Evaluates probability density for CDS threshold uncertainty.

    This function calculates the probability density for a defined function that describes
    uncertainty in the CDS threshold.

    Args:
        x: Array of values at which to calculate density.

    Returns:
        Array of probability density values corresponding to input x values.

    Note:
        The function uses predefined parameters (tl, tu, a, b, c) that were determined
        empirically for CDS threshold uncertainty modeling.
    """
    dq = np.zeros(len(x))
    index = np.where((x > 0.5) & (x < 1))[0]
    tl, tu, a, b, c = 0.5, 1, 0.38387606, -5.40387609, 2.8775016
    g = ((x[index] - tl) / (tu - tl)) ** (-1 / c) - 1
    dg = - ((x[index] - tl) / (tu - tl)) ** (-1 / c - 1) / (c * (tu - tl))
    h = b - np.log(g) / a
    dh = - dg / (a * g)
    dq[index] = np.exp(-h) / (1 + np.exp(-h)) ** 2 * dh
    return dq

def get_minimum_pod_means(stats: Dict[str, np.ndarray], cds_thresholds: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Computes minimum PoD means for specified CDS thresholds.

    Args:
        stats: Dictionary containing probe statistics with keys 'probe', 'pod', 'cds', 'max_conc'.
        cds_thresholds: Array of CDS threshold values to evaluate.

    Returns:
        Tuple containing:
            - min_means: Array of minimum PoD means for each threshold
            - min_probes: Array of probe IDs corresponding to minimum means
            - min_cds: Array of CDS values for the minimum probes

    Raises:
        KeyError: If required keys are missing from stats dictionary.
    """
    min_means = np.full(len(cds_thresholds), stats['max_conc'])
    min_probes = np.full(len(cds_thresholds), 'Max. conc.', dtype='object')
    min_cds = np.full(len(cds_thresholds), 0, dtype='float')

    for i, threshold in enumerate(cds_thresholds):
        mask = stats['cds'] >= threshold
        if np.sum(mask) > 0:
            pod = stats['pod'][mask]
            index = np.argmin(pod)
            min_means[i] = pod[index]
            min_probes[i] = stats['probe'][mask][index]
            min_cds[i] = stats['cds'][mask][index]

    return min_means, min_probes, min_cds

def fit_pod_histogram(pod_samples: np.ndarray, n_samp: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Creates histogram approximation of PoD distribution.

    Args:
        pod_samples: Array of samples from the PoD distribution.
        n_samp: Maximum number of possible samples.

    Returns:
        Tuple containing:
            - weights: Array of histogram weights (or None if no valid samples)
            - bin_edges: Array of histogram bin edges (or None if no valid samples)

    Note:
        If pod_samples contains infinite values, n_samp is set to the length of pod_samples.
    """
    if np.isinf(pod_samples).any():
        n_samp = len(pod_samples)
    pod_samples = pod_samples[~np.isinf(pod_samples)]
    n_pod_samples = len(pod_samples)

    if n_pod_samples == 0:
        return None, None

    prob_response = n_pod_samples / n_samp
    counts, bin_edges = np.histogram(pod_samples, bins=int(np.sqrt(n_pod_samples)))
    bin_size = bin_edges[1:] - bin_edges[:-1]
    bin_weights = counts * bin_size
    total_weight = np.sum(bin_weights)
    weights = counts / total_weight * prob_response

    return weights, bin_edges

def get_global_pod(stats: Dict[str, np.ndarray]) -> Dict[str, Union[float, np.ndarray, int]]:
    """Calculates global PoD from probe-level statistics.

    This function computes the global point of departure (PoD) by aggregating
    probe-level PoD distributions using a weighted approach based on CDS thresholds.

    Args:
        stats: Dictionary containing probe statistics with keys 'probe', 'pod', 'cds'.

    Returns:
        Dictionary containing:
            - global_pod: Calculated global PoD value
            - num_hits: Expected number of hits
            - means: Array of minimum PoD means for each threshold
            - probes: Array of probe IDs corresponding to minimum means
            - weights: Array of weights used in calculation
            - quantiles: Array of CDS threshold quantiles
            - cds: Array of CDS values for minimum probes

    Raises:
        KeyError: If required keys are missing from stats dictionary.
        ValueError: If stats data is malformed.
    """
    logger.info("Starting global PoD calculation...")

    # Get dictionary of minimum probes
    dq = 0.025
    quantiles = np.arange(0.5, 1 + dq, dq)
    quantiles[-1] -= 1e-6
    logger.info(f"Calculating minimum PoD means for {len(quantiles)} quantiles...")
    min_means, min_probes, min_cds = get_minimum_pod_means(stats, quantiles)

    # Calculate weights and global PoD
    logger.info("Calculating confidence threshold probability density...")
    weights = get_confidence_threshold_probability_density(quantiles)
    weight_sum = np.sum(weights)
    global_pod = 10 ** (np.sum(min_means * weights) / weight_sum)
    logger.info(f"Global PoD calculated: {global_pod}")

    # Calculate number of hits at each confidence threshold
    logger.info("Calculating number of hits at each confidence threshold...")
    num_hits = np.array([np.sum(stats['cds'] >= i) for i in quantiles])
    expected_num_hits = np.round(np.sum(num_hits * weights) / weight_sum)
    logger.info(f"Expected number of hits: {expected_num_hits}")

    results = {
        'global_pod': global_pod,
        'num_hits': expected_num_hits,
        'means': min_means,
        'probes': min_probes,
        'weights': weights,
        'quantiles': quantiles,
        'cds': min_cds
    }

    return results

def get_min_probe_weights(df: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Computes weights for probes contributing to global PoD.

    Args:
        df: Dictionary containing global PoD data with keys 'probes', 'means', 'weights', 'cds'.

    Returns:
        Dictionary containing:
            - probe: Array of probe identifiers
            - weight: Array of weights for each probe
            - min_mean: Array of minimum means for each probe
            - cds: Array of CDS values for each probe

    Raises:
        KeyError: If required keys are missing from input dictionary.
    """
    # Aggregate weights for each unique probe
    probes = np.unique(df['probes'])
    means = np.power(10, [df['means'][df['probes'] == i][0] for i in probes])
    cds = np.array([df['cds'][df['probes'] == i][0] for i in probes])
    weights = np.array([np.sum(df['weights'][df['probes'] == i]) / np.sum(df['weights']) for i in probes])

    rank = np.argsort(weights)[::-1]
    probes, means, cds, weights = probes[rank], means[rank], cds[rank], weights[rank]

    weight_dict = {
        'probe': probes,
        'weight': weights,
        'min_mean': means,
        'cds': cds
    }
    return weight_dict

def get_conc_response_plot(df: pd.Series, probe: str, conc_units: str = '\u03bcM') -> go.Figure:
    """Creates concentration-response plot for a specific probe.

    Args:
        df: BIFROST summary data containing probe information.
        probe: Probe identifier to plot.
        conc_units: Concentration units for plot labels (default: 'μM').

    Returns:
        Plotly figure object containing the concentration-response plot.

    Raises:
        KeyError: If probe data is missing from summary.
        ValueError: If probe data is malformed.
    """
    # Extract plotting data
    conc, conc_index = 10 ** np.array(df['conc']), np.array(df['conc_index'])
    x_range = [np.log10(np.min(conc) / 10), np.log10(np.max(conc) * 2)]
    count, total_count = np.array(df[probe]['count']), np.array(df['total_count'])
    treatment_mask = conc_index > 0
    control_mask = conc_index == 0
    median_total_count = np.median(total_count)
    treatment_x = conc[conc_index[treatment_mask] - 1]
    treatment_y = count[treatment_mask] / total_count[treatment_mask] * median_total_count
    control_y = count[control_mask] / total_count[control_mask] * median_total_count
    response_x, response = 10 ** np.array(df[probe]['x']), np.array(df[probe]['response'])
    max_y = max((np.max(control_y), np.max(treatment_y), np.max(response))) * 1.05

    # Create figure
    fig = go.Figure()

    # Add plot elements
    for i in control_y:
        fig.add_hline(y=i, line_width=1, line_dash='dash',
                     line_color=f'rgba(128, 128, 128, 0.5)', layer='below')

    if df[probe]['cds'] > 0:
        mean_pod = 10 ** np.mean(df[probe]['pod'])
        weights, bin_edges = fit_pod_histogram(np.array(df[probe]['pod']), df['n_samp'])
        if weights is not None and bin_edges is not None:
            bin_edges = 10 ** bin_edges
            weights = weights / np.max(weights) * 0.5 * df[probe]['cds']

            for i, weight in enumerate(weights):
                fig.add_shape(
                    type='rect',
                    line=dict(color='rgba(0,0,0,0)'),
                    fillcolor=f'rgba(102, 51, 153, {weight})',
                    layer='below',
                    x0=bin_edges[i],
                    x1=bin_edges[i + 1],
                    xref='x',
                    y0=0,
                    y1=1,
                    yref='paper',
                )

            fig.add_vline(
                line=dict(color='rgba(102, 51, 153, 1)'),
                layer='below',
                x=mean_pod,
            )

    # Add response curves
    fig.add_trace(go.Scatter(
        x=response_x,
        y=response[0],
        mode='lines',
        marker=dict(color='rgba(255, 0, 0, 0.1)'),
        showlegend=False
    ))

    fig.add_trace(go.Scatter(
        x=response_x,
        y=response[2],
        fill='tonexty',
        fillcolor='rgba(255, 0, 0, 0.1)',
        mode='lines',
        marker=dict(color='rgba(255, 0, 0, 0.1)'),
        showlegend=False
    ))

    fig.add_trace(go.Scatter(
        x=response_x,
        y=response[1],
        mode='lines',
        marker=dict(color='red'),
        showlegend=False
    ))

    # Add treatment points
    fig.add_trace(go.Scatter(
        x=treatment_x,
        y=treatment_y,
        mode='markers',
        marker=dict(
            line=dict(color='black', width=2),
            size=10,
            symbol='x-thin',
        ),
        showlegend=False
    ))

    # Update axes and layout
    fig.update_xaxes(
        title_text=f'Concentration ({conc_units})',
        type='log',
        showline=True,
        linewidth=2,
        linecolor='black'
    )
    fig.update_yaxes(
        title_text='Normalised count',
        showline=True,
        linewidth=2,
        linecolor='black'
    )

    # Add subtitle
    if df[probe]['cds'] > 0:
        subtitle_text = (f'CDS = {df[probe]["cds"]:.3f}, '
                        f'Mean PoD | Response = {10**np.mean(df[probe]["pod"]):.2g} {conc_units}')
    else:
        subtitle_text = f'CDS = {df[probe]["cds"]:.3f}, Mean PoD | Response = N/A'

    fig.update_layout(
        autosize=True,
        margin=dict(l=10, r=5, t=5, b=10),
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        paper_bgcolor='rgba(255,255,255,1)',
        plot_bgcolor='rgba(255,255,255,1)',
        xaxis_range=x_range,
        yaxis_range=[0, max_y],
        title=dict(
            text=f'{probe}',
            font=dict(size=24),
            automargin=True,
            yref='container',
            subtitle=dict(
                text=subtitle_text,
                font=dict(color="gray", size=22)
            )
        )
    )

    return fig

def get_pod_vs_fc_scatter(df: pd.Series) -> go.Scatter:
    """Creates scatter plot of PoD values versus log2 fold changes.

    Args:
        df: Series containing PoD statistics with keys 'pod', 'l2fc', 'probe'.

    Returns:
        Plotly Scatter object for the PoD vs fold-change plot.

    Raises:
        KeyError: If required keys are missing from input Series.
    """
    return go.Scatter(
        x=np.power(10, df['pod']),
        y=df['l2fc'],
        text=df['probe'],
        mode='markers',
        marker=dict(
            symbol='circle',
            size=2,
            color='black',
            line=dict(color='black', width=1)
        ),
        showlegend=False
    )

def get_pod_vs_fc_graph(df: pd.Series, global_pod: float, add_annotations: bool = False,
                       conc_units: str = '\u03bcM') -> go.Figure:
    """Creates complete PoD vs fold-change graph with annotations.

    Args:
        df: Series containing PoD statistics.
        global_pod: Global PoD value to mark on plot.
        add_annotations: Whether to add probe labels (default: False).
        conc_units: Concentration units for plot labels (default: 'μM').

    Returns:
        Plotly Figure object containing the complete graph.

    Raises:
        KeyError: If required keys are missing from input Series.
    """
    # Calculate y-axis limits
    if len(df['l2fc']) > 0:
        ymin, ymax = min(df['l2fc'].min(), -2) - 1, max(df['l2fc'].max(), 2) + 1
    else:
        ymin, ymax = -2, 2

    conc = 10 ** np.array(df['conc'])

    # Create figure
    fig = go.Figure()

    # Add base elements
    fig.add_hline(y=0, line_width=1, line_dash='dash', line_color='black')

    if len(df['pod']) > 0:
        fig.add_trace(get_pod_vs_fc_scatter(df))

    for i in conc:
        fig.add_vline(x=i, line_width=1, line_dash='dash',
                     line_color=f'rgba(128, 128, 128, 0.5)')

    fig.add_vline(x=global_pod, line=dict(color='rgba(255, 0, 0, 1)'))

    # Update axes
    fig.update_xaxes(
        title_text=f'Mean PoD | Response ({conc_units})',
        type='log',
        showline=True,
        linewidth=2,
        linecolor='black',
        range=[np.log10(np.min(conc) / 10), np.log10(np.max(conc) * 2)]
    )
    fig.update_yaxes(
        title_text='Max./min. log2 fold-change',
        showline=True,
        linewidth=2,
        linecolor='black',
        range=(ymin, ymax)
    )

    # Update layout
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        xaxis_showgrid=False,
        yaxis_showgrid=False,
        paper_bgcolor='rgba(255,255,255,1)',
        plot_bgcolor='rgba(255,255,255,1)',
    )

    # Add annotations if requested
    if add_annotations and len(df['pod']) > 0:
        rank = np.argsort(df['pod'])
        pod_masked = df['pod'][rank][:10]
        l2fc_masked = df['l2fc'][rank][:10]
        probe_masked = df['probe'][rank][:10]

        order = np.argsort(l2fc_masked)
        for idx, (pod, fc, probe) in enumerate(zip(pod_masked[order],
                                                 l2fc_masked[order],
                                                 probe_masked[order])):
            fig.add_annotation(
                xref="x",
                yref="y",
                x=pod,
                y=fc,
                axref="x",
                ayref="y",
                ax=pod-0.5,
                ay=ymin + (ymax - ymin) / 11 * (idx + 1),
                text=probe,
                showarrow=True,
                arrowhead=1
            )

    return fig

def format_pod(pod: float, max_dose: float) -> str:
    """Formats PoD value for display.

    Args:
        pod: PoD value to format.
        max_dose: Maximum dose value for comparison.

    Returns:
        Formatted string representation of PoD value.
    """
    if pod < max_dose:
        return format_float(10 ** pod)
    return f'>{format_float(10 ** pod)}'

def format_float(x: float) -> str:
    """Formats float to string with two significant figures.

    Args:
        x: Float value to format.

    Returns:
        String representation of float with two significant figures.
    """
    return f'{float(f"{x:.2g}"):g}'

def get_units(units: str, figure: bool = False) -> str:
    """Converts unit identifiers to LaTeX or figure-friendly format.

    Args:
        units: Unit identifier ('uM', 'ugml-1', or 'mgml-1').
        figure: Whether to use figure-friendly format (default: False).

    Returns:
        Formatted unit string.

    Raises:
        ValueError: If unit identifier is not recognized.
    """
    if not figure:
        unit_map = {
            'uM': r'\unit{\micro M}',
            'ugml-1': r'\unit{\micro g.\milli\litre^{-1}}',
            'mgml-1': r'\unit{mg.\milli\litre^{-1}}'
        }
    else:
        unit_map = {
            'uM': '\u03bcM',
            'ugml-1': '\u03bcg mL$^{-1}$',
            'mgml-1': 'mg mL$^{-1}$'
        }

    if units not in unit_map:
        raise ValueError(f'Unit identifier {units} not recognised')

    return unit_map[units]

def initalise_document(test_substance, cell_type):
    """
    Initialises a summary report.

    Accepts:
        chemical (str) - chemical name
        cell_type (str) - cell type

    Returns:
        doc (pylatex.Document) - base document
    """
    # Create report
    geometry_options = {"tmargin": "3.5cm", "lmargin": "2cm",
                       "rmargin": "2cm", "bmargin": "2cm",
                       "headheight": "3cm", "headsep": "1cm"}
    doc = Document(geometry_options=geometry_options)
    doc.packages.append(Package('hyperref'))
    doc.packages.append(Package('siunitx'))
    doc.packages.append(Package('upgreek'))
    doc.preamble.append(NoEscape(r"\hypersetup{colorlinks,"
                                r"citecolor=black,"
                                r"filecolor=black,"
                                r"linkcolor=black,"
                                r"urlcolor=black}"))

    # Create headers and footers
    header = PageStyle("header")
    with header.create(Head("C")) as header_centre:
        header_centre.append("BIFROST HTTr Analysis Summary")

    # Create left footer
    with header.create(Foot("L")):
        header.append(f'{test_substance} {cell_type}')
    # Create center footer
    with header.create(Foot("C")):
        header.append(simple_page_number())

    doc.preamble.append(header)
    doc.change_document_style("header")

    doc.append(NoEscape(
        rf"""
        \begin{{titlepage}}
        \centering
        \vspace{{1cm}}
        {{\huge BIFROST HTTr Analysis Summary\par}}
        \vspace{{1.5cm}}
        {{\huge\bf {test_substance}\par}}
        {{\huge\bf {cell_type}\par}}
        \vspace{{2cm}}
        \vfill
        % Bottom of the page
        {{\large \today\par}}
        \end{{titlepage}}
        """
    ))
    doc.append(NewPage())

    return doc

def add_introduction_page(doc, test_substance, cell_type, timepoint):
    """
    Returns a page with introductory information

    Accepts:
        doc (pylatex.Document) - document on to which summary page is added.
        test_substance (str) - name of the test substance
        cell_type (str) - type of cell used
        timepoint (str) - exposure duration

    Returns:
        doc (pylatex.Document) - updated document
    """
    doc.append(NewPage())
    with doc.create(Section('Introduction')):
        doc.append(NoEscape(
            rf"""
            This report contains analysis of high-throughput transcriptomics data (HTTr) obtained after
            exposing {cell_type} cells for {timepoint} to {test_substance}. The BIFROST model
            (Bayesian inference for region
            of signal threshold) is a statistical model for analysis of HTTr concentration-response data.
            The model is designed to
            infer a point-of-departure (PoD) from a concentration-response dataset.
            The PoD is an estimate of the minimum effect concentration of the test substance
             for the experimental conditions under which the data were produced.
            PoDs are estimated as probability distributions. \\
            The implementation of the approach used here returns a single PoD for each probe analysed.
            PoD distributions are summarised in
            terms of quantiles of the distribution. The concentration-dependency-score (CDS) is the
            inferred probability that the test substance induces a change in expression below the maximum
             concentration tested.\\
            PoD distributions from individual probes are used to calculate a global PoD, defined as an estimate of a
            minimum effect concentration to induce perturbation in expression of any gene. The global PoD is formally
            an \emph{{expectation}}
             with respect to the nominal concentration of the test substance.\\
            1. A summary section including the global PoD and other overall statistics. Included within this section
            is a plot of the median PoD against the maximum log$_2$ fold-change in expression within the
            concentration-range. \\
            2. Concentration-response plots for probes with the 10 lowest expected PoDs. The entire probe set is
             first filtered for probes with a CDS > 0.5. \\
            3. A table summarising PoD statistics for the lowest 200 probes when ranked by the mean of the
            distribution conditional on there being a response. \\
            4. All PoDs are expressed with respect to the nominal concentration of the test substance.
            """
        ))

    return doc

def add_summary_page(doc, stats, stats_filtered, global_pod, weights, conc_units):
    """
    Returns a page detailing overall summary of the BIFROST analysis.
    """
    logger.info("Starting summary page generation...")

    doc.append(NewPage())
    with doc.create(Section('Summary')):
        logger.info("Creating summary table...")
        doc.append(NoEscape(r'\begin{center}'))
        with doc.create(Tabular('c|c')) as table:
            logger.info("Adding table rows...")
            table.add_row((NoEscape(f'Global PoD ({get_units(conc_units)})'),
                         f'{format_float(global_pod["global_pod"])}'))
            table.add_hline()
            table.add_row((NoEscape(f'Maximum tested concentration ({get_units(conc_units)})'),
                         f'{format_float(10 ** stats["max_conc"])}'))
            table.add_hline()
            table.add_row(('Num. probes analysed', f'{len(stats["probe"])}'))
            table.add_hline()
            table.add_row(('Num. "hits"', f'{int(global_pod["num_hits"])}'))
            table.add_hline()
            table.add_row(('Num. CDS>0.5', f'{np.sum(stats["cds"] > 0.5)}'))
            table.add_hline()
            table.add_row(('Num. CDS=1.0', f'{np.sum(stats["cds"] == 1.0)}'))
            table.add_hline()
            logger.info("Adding probe information to table...")
            order = np.argsort(weights['weight'])[::-1]
            for idx, (pr, we, cds, pod) in enumerate(zip(weights['probe'][order],
                                                       weights['weight'][order],
                                                       weights['cds'][order],
                                                       weights['min_mean'][order],
                                                       )):
                if idx == 0:
                    table.add_row(('Minimum responding probe',
                                 NoEscape(
                                     fr'{pr}, '.replace('_', r'\_') +
                                     f'weight={format_float(we)}, '
                                     f'CDS={format_float(cds)}, '
                                     f'Mean PoD={format_float(pod)} '
                                     f'{get_units(conc_units)}')))
                else:
                    table.add_row(('',
                                 NoEscape(
                                     fr'{pr}, '.replace('_', r'\_') +
                                     f'weight={format_float(we)}, '
                                     f'CDS={format_float(cds)}, '
                                     f'Mean PoD={format_float(pod)} '
                                     f'{get_units(conc_units)}')))
            table.add_hline()
            table.add_row(('Largest fold increase', stats['probe'][np.argmax(stats['l2fc'])]))
            table.add_hline()
            table.add_row(('Largest fold decrease', stats['probe'][np.argmin(stats['l2fc'])]))
        doc.append(NoEscape(r'\end{center}'))

        logger.info("Creating PoD vs fold-change plot...")
        with doc.create(Figure(position='ht!')) as figure:
            logger.info("Generating plot figure...")
            fig_path = f'{os.getcwd()}/figures/pod_fold_change.pdf'
            fig_abs_path = os.path.abspath(fig_path)
            logger.info(f"Plot will be saved to: {fig_path}")

            logger.info("Creating plotly figure...")
            fig = get_pod_vs_fc_graph(stats_filtered, global_pod['global_pod'], add_annotations=True,
                                     conc_units=get_units(conc_units, figure=True))

            logger.info("Writing plot to PDF...")
            fig.write_image(fig_path)
            logger.info("Plot PDF created successfully")

            logger.info("Adding plot to LaTeX document...")
            figure.add_image(fig_abs_path, width=NoEscape(r'0.85\linewidth'))
            figure.add_caption(NoEscape(f"""
            Maximum fold-change in expression over the tested concentration-range plotted against the probe-level PoD
            (mean given response). The red vertical dashed line is plotted at the global PoD (nominal concentration).
            Vertical grey lines are
            placed at the experimental test substance concentrations.
            """))
            logger.info("Plot added to document successfully")

    logger.info("Summary page generation complete")
    return doc

def add_concentration_response_plots(doc, probes_to_plot, df, conc_units):
    """
    Adds concentration-response plots to the document.
    """
    for i in range(int(np.ceil(len(probes_to_plot) / 2))):
        probe_set = probes_to_plot[2 * i: 2 * i + 2]

        with doc.create(Figure(position='h!')):
            with doc.create(SubFigure(position='b', width=NoEscape(r'0.45\linewidth'))) as left_fig:
                fig_path = f'{os.getcwd()}/figures/{probe_set[0]}.pdf'
                fig_abs_path = os.path.abspath(fig_path)
                fig = get_conc_response_plot(df, probe_set[0], get_units(conc_units, figure=True))
                fig.write_image(fig_abs_path)
                left_fig.add_image(fig_abs_path, width=NoEscape(r'\linewidth'))

            if len(probe_set) > 1:
                with doc.create(SubFigure(position='b', width=NoEscape(r'0.45\linewidth'))) as right_fig:
                    fig_path = f'{os.getcwd()}/figures/{probe_set[1]}.pdf'
                    fig_abs_path = os.path.abspath(fig_path)
                    fig = get_conc_response_plot(df, probe_set[1], get_units(conc_units, figure=True))
                    fig.write_image(fig_abs_path)
                    right_fig.add_image(fig_abs_path, width=NoEscape(r'\linewidth'))

    return doc

def add_concentration_response_plot_section(doc, df, stats, weights, conc_units):
    """
    Adds pages with concentration-response plots

    Accepts:
        doc (pylatex.Document) - document on to which summary page is added.
        df (pd.Series) - BIFROST summary
        stats (dict) - dictionary of summary statistics
        conc_units (str) - test substance concentration units to display in tables/plots

    Returns:
        doc (pylatex.Document) - updated document
    """
    doc.append(NewPage())
    with doc.create(Section('Concentration-response plots')):
        doc.append(NoEscape(
            r"""
            Concentration-response plots for the lowest 10 responding probes (ranked by mean of distribution
            less than maximum concentration tested given CDS>0.5). Quantities indicated with each plot item are
            as follows: \\
            \textbf{Black markers}: normalised count for treated sample \\
            \textbf{Horizontal grey line}: normalised count for solvent control sample \\
            \textbf{Red bands}: 90\% centred interval indicating plausible range for the expected count \\
            \textbf{Red line (middle)}: median of the distribution of the expected count \\
            \textbf{Purple bands}: histogram-like representation of PoD distribution where band intensity represents
             bin height \\
            \textbf{Purple line}: estimate of mean of PoD distribution, conditional on there being a response \\
            \textbf{CDS}: concentration-dependency-score. Proportion of the PoD distribution less than the maximum
            concentration tested \\
            """
        ))

        with doc.create(Subsection('Probes with non-zero global PoD weight')):
            # Filter out 'Max. conc.' placeholder values
            valid_probes = weights['probe'][weights['probe'] != 'Max. conc.']
            probes_to_plot = valid_probes[np.argsort(weights['weight'][weights['probe'] != 'Max. conc.'])]
            if len(probes_to_plot) > 0:
                doc = add_concentration_response_plots(doc, probes_to_plot, df, conc_units)
            else:
                doc.append(NoEscape(r"\textbf{No probes with non-zero global PoD weight found.}"))

        doc.append(NewPage())
        with doc.create(Subsection('Probes with largest fold changes')):
            index = np.argsort(stats['l2fc'])
            probes_to_plot = np.concatenate((stats['probe'][index[:2]], stats['probe'][index[-2:]]))
            if len(probes_to_plot) > 0:
                doc = add_concentration_response_plots(doc, probes_to_plot, df, conc_units)
            else:
                doc.append(NoEscape(r"\textbf{No probes with fold changes found.}"))

        doc.append(NewPage())
        with doc.create(Subsection('10 lowest means with CDS > 0.5')):
            n_probe = len(stats['probe'])
            probes_to_plot = stats['probe'][np.argsort(stats['pod'])][:min(n_probe, 10)]
            if len(probes_to_plot) > 0:
                doc = add_concentration_response_plots(doc, probes_to_plot, df, conc_units)
            else:
                doc.append(NoEscape(r"\textbf{No probes with CDS > 0.5 found.}"))

    return doc

def add_pod_summary_table(doc, df, stats, conc_units):
    """
    Adds pages with concentration-response plots

    Accepts:
        doc (pylatex.Document) - document on to which summary page is added.
        df (pd.Series) - BIFROST summary
        stats (dict) - dictionary of summary statistics
        conc_units (str) - test substance concentration units to display in tables/plots

    Returns:
        doc (pylatex.Document) - updated document
    """
    doc.append(NewPage())
    with doc.create(Section('Probe-level PoD statistics (CDS > 0.5)')):
        n_probe = len(stats['probe'])
        if n_probe > 0:
            with doc.create(LongTable('|c|c|ccccc|')) as data_table:
                data_table.add_hline()
                data_table.add_row(('Probe', 'CDS',
                                  MultiColumn(5, align='c|',
                                            data=NoEscape(fr'PoD percentiles ({get_units(conc_units)})'))))
                data_table.add_row(('', '', '5th', '25th', '50th', '75th', '95th'))
                data_table.add_hline()
                data_table.end_table_header()
                data_table.add_hline()
                data_table.add_row((MultiColumn(7, align='|c|',
                                              data='Continued on next page'),))
                data_table.add_hline()
                data_table.end_table_footer()
                data_table.add_hline()
                data_table.end_table_last_footer()

                probes = stats['probe'][np.argsort(stats['pod'])][:(min(100, n_probe))]
                cds = stats['cds'][np.argsort(stats['pod'])][:(min(100, n_probe))]
                for i, (j, k) in enumerate(zip(probes, cds)):
                    n_pod_samples = len(df[j]['pod'])
                    extended_pod_samples = np.concatenate((df[j]['pod'],
                                        [df['max_conc'] for _ in range(df['n_samp'] - n_pod_samples)]))
                    pod_percentiles = np.percentile(extended_pod_samples, q=(5, 25, 50, 75, 95))
                    data_table.add_row(np.concatenate(([j, format_float(k)],
                                                     [format_pod(u, df['max_conc'])
                                                      for u in pod_percentiles])))

    return doc

def create_bifrost_httr_report(df, report_name,
                              test_substance, cell_type,
                              timepoint='24 hours', conc_units='uM'):
    """
    Creates a PDF summarising HTTr inferences obtained with the BIFROST method

    Accepts:
        df (pd.Series) - BIFROST summary file
        report_name (str) - filename for compiled PDF
        test_substance (str) - test substance name
        cell_type (str) - cell type name
        timepoint (str) - exposure duration within experiment
        conc_units (str) - test substance concentration units

    Returns:
        None
    """
    logger.info(f"Starting report generation for {test_substance} on {cell_type}")

    # Make directory for storing figures, if not already present
    figures_dir = f'{os.getcwd()}/figures'
    logger.info(f"Creating figures directory at {figures_dir}")
    if not os.path.exists(figures_dir):
        os.mkdir(figures_dir)

    # Compute summary statistics and global PoD
    logger.info("Calculating summary statistics...")
    stats = calculate_summary_statistics(df)
    logger.info(f"Found {len(stats['probe'])} probes in summary statistics")

    logger.info("Filtering summary statistics...")
    stats_filtered = filter_summary_statistics(stats, cds_threshold=0.5)
    logger.info(f"After filtering: {len(stats_filtered['probe'])} probes")

    logger.info("Calculating global PoD...")
    global_pod = get_global_pod(stats)
    logger.info(f"Global PoD calculated: {global_pod['global_pod']}")

    logger.info("Calculating probe weights...")
    weights = get_min_probe_weights(global_pod)
    logger.info(f"Calculated weights for {len(weights['probe'])} probes")

    # Initialize document
    logger.info("Initializing LaTeX document...")
    doc = initalise_document(test_substance, cell_type)

    # Contents page
    logger.info("Adding table of contents...")
    doc.append(NoEscape(r"\tableofcontents"))

    # Introduction page
    logger.info("Adding introduction page...")
    doc = add_introduction_page(doc, test_substance, cell_type, timepoint)

    # Overall summary statistics
    logger.info("Adding summary page...")
    doc = add_summary_page(doc, stats, stats_filtered, global_pod, weights, conc_units)

    # Concentration-response plots for minimum-responding genes
    logger.info("Adding concentration-response plots...")
    doc = add_concentration_response_plot_section(doc, df, stats, weights, conc_units)

    # Table of probe-level PoD quantiles
    logger.info("Adding PoD summary table...")
    doc = add_pod_summary_table(doc, df, stats, conc_units)

    # Generate PDF
    logger.info(f"Generating PDF report: {report_name}.pdf")
    doc.generate_pdf(f'{report_name}', clean_tex=True)
    logger.info("Report generation complete")

def parse_args() -> argparse.Namespace:
    """Parse command line arguments for report generation.

    Returns:
        Namespace containing parsed arguments:
            - summary_file: Path to summary JSON file
            - test_substance: Name of test substance (default: 'MyChemical')
            - cell_type: Type of cell used (default: 'MyCell')
            - output_name: Name for output report (default: 'test_report')
            - timepoint: Exposure duration (default: '24 hours')
            - conc_units: Concentration units (default: 'uM', choices: ['uM', 'ugml-1', 'mgml-1'])
    """
    parser = argparse.ArgumentParser(description='Create Bifrost HTTR reports from summary data')
    parser.add_argument('--summary-file',
                      help='Path to the summary JSON file (default: summary.json.zip)')
    parser.add_argument('--test-substance',
                      default='MyChemical',
                      help='Name of the test substance (default: MyChemical)')
    parser.add_argument('--cell-type',
                      default='MyCell',
                      help='Type of cell used in the test (default: MyCell)')
    parser.add_argument('--output-name',
                      default='test_report',
                      help='Name for the output report (default: test_report)')
    parser.add_argument('--timepoint',
                      default='24 hours',
                      help='Exposure duration within experiment (default: 24 hours)')
    parser.add_argument('--conc-units',
                      default='uM',
                      choices=['uM', 'ugml-1', 'mgml-1'],
                      help='Concentration units (default: uM)')
    return parser.parse_args()

def main() -> None:
    """Main entry point for report generation script.

    This function:
    1. Parses command line arguments
    2. Loads the summary data from JSON
    3. Generates the BIFROST HTTR report

    Raises:
        FileNotFoundError: If summary file doesn't exist
        ValueError: If summary file is malformed
        RuntimeError: If report generation fails
    """
    args = parse_args()

    try:
        # Load summary file
        df = pd.read_json(args.summary_file,
                         typ='series', orient='index', compression='zip')

        create_bifrost_httr_report(df, args.output_name, args.test_substance, args.cell_type,
                                  timepoint=args.timepoint, conc_units=args.conc_units)
    except FileNotFoundError as e:
        logger.error(f"Summary file not found: {e}")
        raise
    except ValueError as e:
        logger.error(f"Invalid summary file format: {e}")
        raise
    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        raise RuntimeError(f"Failed to generate report: {e}")

if __name__ == '__main__':
    main()

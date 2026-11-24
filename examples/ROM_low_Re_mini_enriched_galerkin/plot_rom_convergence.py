#!/usr/bin/env python3
"""Plot full-order and inner ROM convergence from Neko-TOP optimisation CSVs.

Usage:
    python plot_rom_convergence.py

Place this script alongside optimization_data.csv. The input and output
paths are fixed relative to the script, regardless of the working directory. The largest recorded inner iteration is used to
space the ROM trajectories. Individual objective components are unweighted.

The CSV convention is that the ROM reference (inner=0) at outer k belongs
with positive inner iterations labelled outer k+1. A cycle with N planned
inner steps occupies [k,k+1], with inner j plotted at k+j/N.
"""
from pathlib import Path
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

# Okabe–Ito palette, in the same order for every objective.
PALETTE = ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9')
OBJECTIVES = {
    'total_objective': r'Total objective, $\mathcal{J}$',
    'scalar_mixing': 'Scalar mixing',
    'viscous_dissipation': 'Viscous dissipation',
    'brinkman_dissipation': 'Brinkman dissipation',
}


def read_data(path):
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    required = {'source', 'outer_iteration', 'inner_iteration', *OBJECTIVES}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'Missing CSV columns: {sorted(missing)}')
    df['source'] = df['source'].astype(str).str.strip().str.lower()
    for col in ('outer_iteration', 'inner_iteration'):
        df[col] = pd.to_numeric(df[col], errors='raise')
        if not np.isfinite(df[col]).all() or (df[col] < 0).any() or (df[col] % 1 != 0).any():
            raise ValueError(f'{col} must contain non-negative integers')
        df[col] = df[col].astype(int)
    for col in (*OBJECTIVES, *(f'{c}.weight' for c in OBJECTIVES if c != 'total_objective')):
        if col in df:
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(r'[dD]', 'E', regex=True), errors='raise')
    return df


def collect_cycles(df):
    fom = df[df.source == 'full_order'].sort_values('outer_iteration')
    rom = df[df.source == 'rom']
    if fom.empty:
        raise ValueError('No full_order rows found')
    if fom.outer_iteration.duplicated().any():
        raise ValueError('Duplicate full-order outer iterations; select a single run')
    if rom.duplicated(['outer_iteration', 'inner_iteration']).any():
        raise ValueError('Duplicate ROM iteration indices; select a single run')
    refs = rom[rom.inner_iteration == 0].sort_values('outer_iteration')
    updates = rom[rom.inner_iteration > 0]
    max_inner = int(updates.inner_iteration.max()) if not updates.empty else 1
    nsteps = max_inner
    cycles = []
    # Include legacy logs without reference rows, but do not invent an inner-0 value.
    starts = sorted(set(refs.outer_iteration) | {int(k)-1 for k in updates.outer_iteration})
    for k in starts:
        if k < 0:
            raise ValueError('Unexpected positive inner iteration at outer 0')
        ref = refs[refs.outer_iteration == k]
        steps = updates[updates.outer_iteration == k+1].sort_values('inner_iteration')
        if ref.empty:
            warnings.warn(f'ROM cycle {k+1}: no reference at outer {k}; plotting recorded updates only')
        rows = pd.concat([ref, steps], ignore_index=True)
        x = k + rows.inner_iteration.to_numpy()/nsteps
        cycles.append((k+1, x, rows))
    return fom, cycles


def plot_objective(fom, cycles, column, outdir):
    label = OBJECTIVES[column]

    def values(rows):
        return rows[column].to_numpy()

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(fom.outer_iteration, values(fom), '-o', color='#222222',
            linewidth=2.4, markersize=6, label='Full-order objective', zorder=5)
    any_ref = False
    for cycle_id, x, rows in cycles:
        colour = PALETTE[(cycle_id-1) % len(PALETTE)]
        y = values(rows)
        ax.plot(x, y, '--s', color=colour, linewidth=1.8, markersize=5,
                label=f'ROM cycle {cycle_id}', zorder=3)
        if rows.iloc[0].inner_iteration == 0:
            any_ref = True
            ax.plot(x[0], y[0], marker='s', markerfacecolor='white',
                    markeredgecolor=colour, markeredgewidth=1.5, markersize=8.5,
                    linestyle='none', zorder=4)
    handles, labels = ax.get_legend_handles_labels()
    if any_ref:
        handles.append(Line2D([], [], marker='s', markerfacecolor='white',
                              markeredgecolor='#555555', markeredgewidth=1.5,
                              markersize=7.5, linestyle='none'))
        labels.append('ROM reference (inner 0)')
    last = max([int(fom.outer_iteration.max())] + [cid for cid, _, _ in cycles])
    ax.set(xlabel='Outer optimisation iteration', ylabel=label,
           xticks=np.arange(int(fom.outer_iteration.min()), last+1))
    ax.grid(alpha=0.20)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles, labels, frameon=False, loc='best', fontsize=9)
    fig.tight_layout()
    stem = column
    for ext in ('png', 'svg'):
        fig.savefig(outdir / f'{stem}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    directory = Path(__file__).resolve().parent
    source = directory / 'optimization_data.csv'
    outdir = directory / 'rom_convergence_plots'
    df = read_data(source)
    fom, cycles = collect_cycles(df)
    outdir.mkdir(parents=True, exist_ok=True)
    for column in OBJECTIVES:
        plot_objective(fom, cycles, column, outdir)
    print(f'Saved four PNG and SVG figures to {outdir}')


if __name__ == '__main__':
    main()


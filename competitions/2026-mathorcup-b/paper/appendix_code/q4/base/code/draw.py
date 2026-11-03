from pathlib import Path
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = ['#777777', '#4477AA', '#66AABB', '#AA3377']
METHODS = ['A_single_fill', 'B_MF_fill', 'Bm_MF_count_matched', 'C_MF_active']
LABELS = ['Single GP', 'MF passive', 'MF matched', 'MF active']


def fingerprint(path):
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size}


def sol():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    spec_path = args.root / 'results/PLOT_INPUTS.json'
    spec = json.loads(spec_path.read_text())
    inputs = [fingerprint(spec_path)]
    figures = []
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 9,
                         'axes.labelsize': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.frameon': False, 'axes.linewidth': .8})

    def asset(key):
        path = args.root / spec[key]
        if not path.resolve().is_relative_to(args.root.resolve()):
            raise ValueError('plot input must be inside delivery root')
        inputs.append(fingerprint(path))
        return path

    def canvas():
        return plt.subplots(figsize=(7.2, 4.2), layout='constrained')

    def save(fig, name, note):
        fig.canvas.draw()
        (args.out / (name + '.alignment.json')).write_text(json.dumps({
            'verdict': 'NOT APPLICABLE', 'reason': 'one scientific plot-area', 'plot_areas': 1}))
        fig.savefig(args.out / (name + '.pdf'), bbox_inches='tight')
        fig.savefig(args.out / (name + '.svg'), bbox_inches='tight')
        fig.savefig(args.out / (name + '.png'), dpi=600, bbox_inches='tight')
        figures.append({'name': name, 'interpretation': note,
                        'files': [fingerprint(args.out / (name + '.' + e)) for e in ('pdf', 'svg', 'png')]})
        plt.close(fig)

    if 'morris' in spec:
        d = json.loads(asset('morris').read_text())
        fig, ax = canvas()
        ax.scatter(d['mu_star'], d['sigma'], color=COLORS[1], s=40)
        for x, y, label in zip(d['mu_star'], d['sigma'], ['V', 'T', 'μ', 'ξ']):
            ax.annotate(label, (x, y), xytext=(6, 5), textcoords='offset points')
        ax.set(xlabel='Mean absolute elementary effect (m)', ylabel='Elementary-effect SD (m)',
               xlim=(0, max(d['mu_star']) * 1.2), ylim=(0, max(d['sigma']) * 1.2))
        save(fig, '01_morris', 'Eight complete trajectories; domain-dependent screening of raw depth margin (m).')

    if 'hydrograph' in spec:
        d = np.load(asset('hydrograph'))
        series = d['series']
        support_rows = []
        for coordinate in (348660., 348720., 348780.):
            x, y = d['x'], d['y']
            j = int(np.searchsorted(x, coordinate))
            a = (coordinate-x[j-1])/(x[j]-x[j-1])
            selected = (y >= 3134880) & (y <= 3135600)
            h = (1-a)*d['snapshot_U'][:, :, j-1, 0] + a*d['snapshot_U'][:, :, j, 0]
            peak = h[:, selected].max(axis=1)
            hits = np.flatnonzero(peak >= .1)
            support_rows.append({'target_x_m': coordinate,
                'snapshot_peak_h_lower_bound_m': float(peak.max()),
                'first_recorded_exceedance_time_s': float(d['snapshot_t'][hits[0]]) if len(hits) else None})
        (args.out/'TARGET_SUPPORT_SENSITIVITY.json').write_text(json.dumps({
            'results': support_rows, 'nominal_snapshot_interval_s': 30,
            'interpretation': 'Saved-field spatial support comparison; peak lower bounds and first observed snapshot exceedances, not stepwise first-arrival brackets.'}, indent=2))
        for column, name, ylabel in [(1, '02_target_depth', 'Target peak depth (m)'),
                                     (4, '03_target_discharge', 'Signed downstream discharge (m³/s)')]:
            fig, ax = canvas()
            ax.plot(series[:, 0], series[:, column], color=COLORS[1], lw=1.5)
            if column == 1:
                ax.axhline(.1, color=COLORS[0], ls='--', lw=.8, label='Arrival threshold: 0.1 m')
                ax.legend(loc='upper right')
            ax.set(xlabel='Time after prescribed entry (s)', ylabel=ylabel)
            save(fig, name, 'Stored same-time cross-section series; signed discharge and peak depth have distinct meanings.')

    if 'metrics' in spec:
        rows = list(csv.DictReader(asset('metrics').open()))
        for stratum, suffix in [('global', 'global'), ('development_predicted_boundary', 'band')]:
            fig, ax = canvas()
            for seed in sorted({int(r['seed']) for r in rows}):
                selected = [next(r for r in rows if r['method'] == method and int(r['seed']) == seed and r['stratum'] == stratum) for method in METHODS]
                values = [float(r['disagreement_rate']) if r['disagreement_rate'] else np.nan for r in selected]
                ax.plot(range(4), values, color='#BBBBBB', lw=.8, zorder=1)
                ax.scatter(range(4), values, c=COLORS, s=35, zorder=2, marker={137:'o',271:'s',503:'^'}[seed], label=f'Seed {seed}')
            ax.set_xticks(range(4), LABELS)
            ax.set(ylabel='Classification disagreement / valid HF cases', ylim=(-.025, 1.025))
            ax.legend(loc='upper right')
            save(fig, '04_comparison_' + suffix, 'Single matched seed; no cross-seed stability estimate; all valid high-layer cases retained, unknowns excluded only from valid-case denominator.')

    if 'costs' in spec:
        rows = list(csv.DictReader(asset('costs').open()))
        fig, ax = canvas()
        for j, method in enumerate(METHODS):
            v = [float(r['physical_child_wall_sum_s']) / 60 for r in rows if r['method'] == method]
            ax.scatter(np.full(len(v), j), v, color=COLORS[j], s=35)
        ax.set_xticks(range(4), LABELS)
        ax.set(ylabel='Required physical child wall-time sum (min)')
        save(fig, '05_cost', 'Sum of measured required physical child wall times only; GP/acquisition timings missing, so this is not full cost, CPU time, or parallel elapsed time.')

    if 'surface' in spec:
        d = np.load(asset('surface'))
        a = d['axis']; v = 2e5 * 4**a; mu = .03 + .12*a
        fig, ax = canvas()
        img = ax.pcolormesh(v, mu, d['mean_g'], shading='auto', cmap='RdBu_r')
        if d['mean_g'].min() <= 0 <= d['mean_g'].max():
            ax.contour(v, mu, d['mean_g'], levels=[0], colors='black', linewidths=1)
        fig.colorbar(img, ax=ax, label='Frozen GP mean log margin g')
        ax.set(xlabel='Prescribed volume (m³)', ylabel='Coulomb coefficient mu')
        save(fig, '06_surface', 'One frozen active model at T=960 s and xi=500 m/s²; zero contour is a surrogate boundary, not a natural probability.')

    if 'full_field' in spec:
        d = np.load(asset('full_field'))
        fig, ax = canvas()
        xx, yy = np.meshgrid(d['x'], d['y'])
        wet = np.ma.masked_less(d['max_h'], .1)
        img = ax.pcolormesh(xx/1000, yy/1000, wet, shading='auto', cmap='Blues', rasterized=True)
        if 'full_route' in spec:
            route = np.loadtxt(asset('full_route'), delimiter=',', skiprows=1)
            ax.plot(route[:,0]/1000, route[:,1]/1000, color='#666666', lw=.6, label='Reference axis')
        ax.scatter([350.640], [3135.040], marker='>', s=35, color='#AA3377', label='Prescribed entry', zorder=5, clip_on=False)
        ax.plot([340.500, 341.100], [3129.103736725446]*2, color='#AA3377', lw=2, label='Port target')
        ax.legend(loc='lower right', fontsize=8)
        ax.set(xlabel='UTM easting (km)', ylabel='UTM northing (km)', aspect='equal',
               xlim=((d['x'].min()-20)/1000, (d['x'].max()+20)/1000),
               ylim=((d['y'].min()-20)/1000, (d['y'].max()+20)/1000))
        fig.colorbar(img, ax=ax, label='Space-time maximum depth ≥ 0.1 m (m)')
        save(fig, '07_full_footprint', 'Ever-threshold depth footprint differs from final instantaneous inundation; axis line is a geometric reference.')

    if 'inverse' in spec:
        d = np.load(asset('inverse'))
        fig, ax = canvas()
        ax.scatter(d['normalized_arrival_error'], d['normalized_depth_error'], s=4, color='#BBBBBB', rasterized=True)
        indices = d['pareto_indices'].astype(int)
        ax.plot(d['normalized_arrival_error'][indices], d['normalized_depth_error'][indices], color=COLORS[3], lw=1.2)
        ax.axvline(1, ls='--', color=COLORS[0], lw=.8)
        ax.axhline(1, ls='--', color=COLORS[0], lw=.8)
        ax.set(xlabel='Absolute arrival error / 60 s', ylabel='Absolute peak-depth error / 0.75 m')
        save(fig, '08_inverse', 'Synthetic truth at independent 15 m discretization; candidate surface from 30 m library, with separate physical representative rechecks.')

    (args.out / 'PLOT_RECEIPT.json').write_text(json.dumps({'source': fingerprint(Path(__file__)), 'inputs': inputs,
          'figures': figures, 'solver_executed': False, 'GP_fitted': False,
          'visual_and_collision_review': 'pending'}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    sol()

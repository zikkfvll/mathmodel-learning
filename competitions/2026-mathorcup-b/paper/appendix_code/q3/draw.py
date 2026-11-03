from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import argparse, datetime, hashlib, json, math, os, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
STYLE = {
    'font_family': 'sans-serif',
    'fonts': ['Arial', 'Noto Sans CJK SC', 'Microsoft YaHei', 'PingFang SC', 'DejaVu Sans'],
    'font_size': 8,
    'tick_size': 7,
    'title_size': 9,
    'legend_size': 7,
    'panel_label_size': 10,
    'status_size': 7,
    'supertitle_size': 9,
    'panel_label_offset_pt': (-22, 8),
    'panel_label_anchor': (0, 1),
    'panel_label_weight': 'bold',
    'title_pad_pt': 10,
    'panel_xycoords': 'axes fraction',
    'panel_textcoords': 'offset points',
    'panel_annotation_clip': False,
    'status_horizontal_alignment': 'center',
    'status_vertical_alignment': 'center',
    'pending_text_anchor': (0.5, 0.5),
    'q4_status_x_fraction': 0.5,
    'supertitle_y': 0.98,
    'sizes': {
        'grid': (7.2, 5.4),
        'row': (7.2, 3.2),
        'validation': (7.2, 3.4),
        'sensitivity_row': (7.6, 3.6),
        'sensitivity_grid': (7.6, 6.0),
    },
    'figures': {
        'response': {
            'shape': (2, 2),
            'size': 'grid',
            'legend_columns': 3,
            'layout': {'left': 0.1, 'right': 0.98, 'bottom': 0.1, 'top': 0.86, 'wspace': 0.29, 'hspace': 0.45},
        },
        'learning': {
            'shape': (2, 2),
            'size': 'grid',
            'legend_columns': 4,
            'layout': {'left': 0.1, 'right': 0.98, 'bottom': 0.1, 'top': 0.86, 'wspace': 0.29, 'hspace': 0.4},
        },
        'validation': {
            'shape': (1, 2),
            'size': 'validation',
            'legend_columns': None,
            'layout': {'left': 0.1, 'right': 0.98, 'bottom': 0.2, 'top': 0.78, 'wspace': 0.3},
        },
        'distribution': {
            'shape': (2, 2),
            'size': 'grid',
            'legend_columns': 4,
            'layout': {'left': 0.17, 'right': 0.98, 'bottom': 0.1, 'top': 0.85, 'wspace': 0.48, 'hspace': 0.43},
        },
        'lead': {
            'shape': (1, 2),
            'size': 'row',
            'legend_max_columns': 3,
            'layout': {'left': 0.1, 'right': 0.98, 'bottom': 0.21, 'top': 0.76, 'wspace': 0.3},
        },
        'sensitivity_phase': {
            'shape': (1, 2),
            'size': 'sensitivity_row',
            'legend_columns': 2,
            'layout': {'left': 0.09, 'right': 0.91, 'bottom': 0.2, 'top': 0.77, 'wspace': 0.4},
        },
        'sensitivity_limits': {
            'shape': (2, 2),
            'size': 'sensitivity_grid',
            'legend_columns': 2,
            'layout': {'left': 0.12, 'right': 0.97, 'bottom': 0.1, 'top': 0.86, 'wspace': 0.4, 'hspace': 0.55},
        },
        'q4': {
            'shape': (1, 2),
            'size': 'row',
            'layout': {'left': 0.13, 'right': 0.98, 'bottom': 0.2, 'top': 0.76, 'wspace': 0.3},
        },
    },
    'colors': {
        'truth': '#252525',
        'baseline': '#B57A34',
        'passive': '#777777',
        'active': '#246A98',
        'failure': '#A94C47',
        'class_threshold': '#aaaaaa',
        'fallback_legend': '#333333',
        'bound_connector': '#999999',
    },
    'line_styles': {'truth': '-', 'baseline': ':', 'passive': '--', 'active': '-'},
    'mode_styles': {'raw': {'line': '-', 'marker': 'o'}, 'fallback': {'line': '--', 'marker': 's'}},
    'query_markers': {'passive': 'o', 'active': 'x'},
    'lead_markers': ['o', 's', '^', 'D', 'v'],
    'lead_line_styles': ['-', '--', ':', '-.'],
    'lead_color_order': ['truth', 'active', 'baseline', 'passive', 'failure'],
    'line_width': 1.2,
    'marker_size': 3,
    'line_widths': {
        'score': 1.1,
        'class_threshold': 0.7,
        'baseline_reference': 1,
        'validation_target': 1,
        'bound_connector': 0.8,
        'distribution_truth': 1.3,
        'distribution_baseline': 1.4,
        'distribution_gp': 1,
        'interval': 2,
        'lead': 1,
        'local_zero': 1,
        'q4_zero': 0.9,
    },
    'mark_sizes': {'validation': 4, 'interval_endpoint': 7, 'lead': 3},
    'scatter': {
        'query_area': 15,
        'query_edge_width': 0.8,
        'query_zorder': 3,
        'class_display_jitter': 0.025,
        'nohit_area': 30,
    },
    'interval_marker': '|',
    'nohit_marker': 's',
    'errorbar_capsize': 2,
    'hollow_marker_facecolor': 'none',
    'legend_marker_line': '',
    'default_rc': {
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'axes.unicode_minus': True,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'text.usetex': False,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
    },
    'views': {
        'offset_xlim': (-1.55, 1.55),
        'class_ylim': (-0.13, 1.13),
        'class_yticks': [0, 0.5, 1],
        'miss_bottom': -0.06,
        'miss_minimum_top': 1,
        'miss_headroom': 0.15,
        'validation_xlim': (-0.25, 2.25),
        'validation_lower_floor': 0.9,
        'validation_lower_padding': 0.01,
        'validation_upper': 1.005,
        'cdf_xlim': {'strong': (0.765, 0.825), 'near': (71.3, 73)},
        'cdf_ylim': (0, 1.02),
    },
    'sensitivity': {
        'phase_cmap': 'cividis',
        'no_hit_color': '#D9D9D9',
        'phase_shading': 'flat',
        'phase_rasterized': False,
        'colorbar_fraction': 0.05,
        'colorbar_pad': 0.025,
        'nominal_marker': 'o',
        'nominal_area': 35,
        'nominal_face': 'white',
        'nominal_edge_width': 0.8,
        'nominal_zorder': 4,
        'unknown_marker': 'x',
        'initial_marker': '^',
        'category_area': 15,
        'category_edge_width': 0.8,
        'bar_width': 0.55,
        'amplification_area': 35,
        'amplification_marker': 'o',
        'amplification_bottom_factor': 1.6,
        'amplification_top_factor': 1.8,
        'amplification_ratio_anchor': (0.06, 0.92),
        'amplification_ratio_ha': 'left',
        'regime_tick_rotation': 0,
        'chi_xlim_padding': 0.01,
        'log_tick_size': 8,
        'nohit_strip_y': 1.02,
        'nohit_strip_area': 23,
        'chi_legend_loc': 'upper right',
        'chi_legend_anchor': (1, 0.96),
        'condition_markers': {'strong': 'o', 'near': 's'},
        'condition_lines': {'strong': '--', 'near': '-'},
        'width_legend_loc': 'upper left',
    },
    'grid': False,
    'dpi': 300,
    'formats': ['.png', '.pdf', '.svg'],
    'alignment_tolerance_pt': 1.5,
    'legend_anchor': (0.54, 0.99),
    'legend': {'frameon': False, 'loc': 'upper center'},
}
REGIMES = ['strong', 'near']
NAMES = {'strong': 'Strong shear', 'near': 'Near boundary'}
BASE_KEYS = {'strong': 'Q2_strong_design', 'near': 'Q2_near_boundary_comparison'}
FAMILIES = ['nominal_triangular', 'uniform_support', 'lower_mode_triangular']
FAMILY_LABELS = ['Nominal', 'Uniform', 'Lower mode']
REQUIRED = ['summary.json', 'baseline_results.json', 'gp_metrics.csv', 'gp_predictions.csv.gz', 'baseline_test_metrics.json', 'validation.json', 'lead_margins.csv', 'gp_training.json', 'test_truth.csv', 'cost_comparison.csv']

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, obj):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def load_json(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError) as e:
        raise ValueError(f'Cannot read JSON input {Path(p).name}: {e}') from e

def need_columns(frame, fields, filename):
    missing = sorted(set(fields) - set(frame.columns))
    if missing:
        raise ValueError(f'{filename}: missing required columns {missing}')

def apply_style():
    plt.rcParams.update({
        'font.family': STYLE['font_family'],
        'font.sans-serif': STYLE['fonts'],
        'font.size': STYLE['font_size'],
        'axes.labelsize': STYLE['font_size'],
        'axes.titlesize': STYLE['title_size'],
        'xtick.labelsize': STYLE['tick_size'],
        'ytick.labelsize': STYLE['tick_size'],
        'legend.fontsize': STYLE['legend_size'],
        **STYLE['default_rc'],
        'axes.grid': STYLE['grid'],
        'savefig.dpi': STYLE['dpi'],
    })

def panel(ax, label, title=None):
    ax.annotate(label, STYLE['panel_label_anchor'], xycoords=STYLE['panel_xycoords'], xytext=STYLE['panel_label_offset_pt'], textcoords=STYLE['panel_textcoords'], weight=STYLE['panel_label_weight'], fontsize=STYLE['panel_label_size'], annotation_clip=STYLE['panel_annotation_clip'])
    if title:
        ax.set_title(title, pad=STYLE['title_pad_pt'])

def new_figure(name):
    config = STYLE['figures'][name]
    fig, axes = plt.subplots(*config['shape'], figsize=STYLE['sizes'][config['size']])
    fig.subplots_adjust(**config['layout'])
    return (fig, axes)

def require_matplotlib_panel_alignment(fig, path, tolerance_pt=1.5):
    fig.canvas.draw()
    w, h = fig.get_size_inches() * 72
    rects = []
    for i, ax in enumerate(fig.axes):
        if getattr(ax, '_q3_auxiliary', False):
            continue
        p = ax.get_position()
        rects.append(dict(id=str(i), left=p.x0 * w, bottom=p.y0 * h, right=p.x1 * w, top=p.y1 * h, width=p.width * w, height=p.height * h))
    checks = []
    for i, a in enumerate(rects):
        for b in rects[i + 1:]:
            same_row = abs(a['bottom'] - b['bottom']) < tolerance_pt and abs(a['top'] - b['top']) < tolerance_pt
            same_col = abs(a['left'] - b['left']) < tolerance_pt and abs(a['right'] - b['right']) < tolerance_pt
            if same_row:
                checks.append(abs(a['height'] - b['height']) <= tolerance_pt and abs(a['width'] - b['width']) <= tolerance_pt)
            if same_col:
                checks.append(abs(a['width'] - b['width']) <= tolerance_pt)
            overlap_x = min(a['right'], b['right']) - max(a['left'], b['left'])
            overlap_y = min(a['top'], b['top']) - max(a['bottom'], b['bottom'])
            checks.append(not (overlap_x > tolerance_pt and overlap_y > tolerance_pt))
    answer = dict(status='pass' if all(checks) else 'fail', tolerance_pt=tolerance_pt, scope='rectangular comparable plot-area edges/size/nonoverlap only', width_pt=float(w), height_pt=float(h), panels=rects, checks=len(checks), exemptions=[{
        'axis': i,
        'reason': 'auxiliary colorbar encodes the same field; excluded from comparable scientific-panel width/edge groups',
    } for i, ax in enumerate(fig.axes) if getattr(ax, '_q3_auxiliary', False)])
    save(path, answer)
    if not all(checks):
        raise ValueError('Final subplot plot-area alignment failed')

def export(fig, name, out, audit):
    require_matplotlib_panel_alignment(fig, audit / f'{name}.alignment.json', STYLE['alignment_tolerance_pt'])
    for ext in STYLE['formats']:
        fig.savefig(out / f'{name}{ext}', dpi=STYLE['dpi'])
    plt.close(fig)

def read_inputs(results, sensitivity_results=None):
    for f in REQUIRED:
        if not (results / f).is_file():
            raise ValueError(f'Missing required saved analysis result: {f}')
    d = {k: load_json(results / f'{k}.json') for k in ['summary', 'baseline_results', 'baseline_test_metrics', 'validation', 'gp_training']}
    for key, file in [('metrics', 'gp_metrics.csv'), ('pred', 'gp_predictions.csv.gz'), ('lead', 'lead_margins.csv'), ('truth', 'test_truth.csv'), ('cost', 'cost_comparison.csv')]:
        try:
            d[key] = pd.read_csv(results / file)
        except (OSError, ValueError) as e:
            raise ValueError(f'Cannot read {file}: {e}') from e
    need_columns(d['metrics'], ['regime', 'method', 'budget', 'mode', 'subset', 'mae_s', 'false_negative', 'unknown_predictions', 'total_nominal_q2_calls'], 'gp_metrics.csv')
    need_columns(d['pred'], ['regime', 'method', 'budget', 'mode', 'd', 'status', 'tau_s', 'p_hit', 'true_status', 'true_tau_s'], 'gp_predictions.csv.gz')
    need_columns(d['truth'], ['regime', 'd', 'status', 'tau_s'], 'test_truth.csv')
    if not {'edges_s', 'conditions'} <= set(d['baseline_results']):
        raise ValueError('baseline_results.json: edges_s and conditions required')
    edges = np.asarray(d['baseline_results']['edges_s'], float)
    if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
        raise ValueError('baseline time edges must be finite and strictly increasing')
    for key in ['truth', 'pred']:
        frame = d[key]
        if not frame.status.isin(['hit', 'no_hit', 'initial', 'unknown']).all():
            raise ValueError(f'{key}: invalid outcome status')
        if frame.loc[frame.status != 'hit', 'tau_s'].notna().any():
            raise ValueError(f'{key}: non-hit outcome cannot carry future first-passage time')
        times = frame.loc[frame.status == 'hit', 'tau_s']
        if times.isna().any() or ((times <= 0) | (times > 73)).any():
            raise ValueError(f'{key}: invalid hit time')
    d['lead'] = normalize_lead(d['lead'])
    d['validation_rows'] = normalize_validation(d['validation'])
    for regime in REGIMES:
        sp = spec_of(base_condition(d, regime), 'robust_shortest')
        rows = d['lead'][d['lead'].regime == regime]
        if not len(rows):
            raise ValueError(f'lead_margins.csv: missing {regime} delay scenarios')
        interval = sp.get('finite_interval_s')
        if interval is not None:
            shift = rows.publication_delay_s + rows.action_delay_s
            if not np.allclose(rows.finite_lower_s + shift, interval[0], atol=1e-09, rtol=0) or not np.allclose(rows.finite_upper_s + shift, interval[1], atol=1e-09, rtol=0):
                raise ValueError(f'lead_margins.csv: {regime} margins are not the formal-set translation')
        if not (rows.includes_no_hit.astype(bool) == bool(sp.get('includes_no_hit_by_horizon'))).all():
            raise ValueError(f'lead_margins.csv: {regime} no-hit category changed during translation')
    cost_alias = {
        'construction_wall_s': 'build_s',
        'query_wall_s': 'query_s',
        'probability_set_wall_s': 'probability_s',
        'first_use_total_wall_s': 'first_use_total_s',
        'repeat_query_total_wall_s': 'repeat_use_total_s',
    }
    d['cost'] = d['cost'].rename(columns={a: b for a, b in cost_alias.items() if a in d['cost'] and b not in d['cost']})
    cost = d['cost']
    cost_fields = ['build_s', 'query_s', 'probability_s', 'first_use_total_s', 'repeat_use_total_s']
    need_columns(cost, ['regime', 'method', 'budget', 'mode', 'construction_q2_calls', 'fallback_q2_calls'] + cost_fields, 'cost_comparison.csv')
    if not np.isfinite(cost[cost_fields].to_numpy(float)).all() or (cost[cost_fields] < 0).any().any():
        raise ValueError('cost_comparison.csv: required workload costs must be finite and nonnegative')
    if not np.allclose(cost.first_use_total_s, cost.build_s + cost.query_s + cost.probability_s, atol=1e-08, rtol=1e-09):
        raise ValueError('cost_comparison.csv: first-use total does not equal construction plus same queries plus probability task')
    if not np.allclose(cost.repeat_use_total_s, cost.query_s + cost.probability_s, atol=1e-08, rtol=1e-09):
        raise ValueError('cost_comparison.csv: repeated workload total is inconsistent')
    if (results / 'q4_spatial_demo.csv').is_file():
        d['q4_demo'] = normalize_q4(pd.read_csv(results / 'q4_spatial_demo.csv'))
        for regime in REGIMES:
            rows = d['q4_demo'][(d['q4_demo'].regime == regime) & (d['q4_demo'].status == 'arrived')]
            sp = spec_of(base_condition(d, regime), 'robust_shortest')
            interval = sp['finite_interval_s']
            if not np.allclose(rows.finite_arrival_lower_s - rows.travel_time_s, interval[0], atol=1e-08, rtol=0) or not np.allclose(rows.finite_arrival_upper_s - rows.travel_time_s, interval[1], atol=1e-08, rtol=0):
                raise ValueError('q4_spatial_demo.csv: assumed zero bridge must translate the formal local candidate exactly')
            if not (rows.includes_no_hit.astype(bool) == bool(sp.get('includes_no_hit_by_horizon'))).all():
                raise ValueError('q4_spatial_demo.csv: local no-hit category was lost')
    else:
        d['q4_demo'] = None
    d['sensitivity'] = read_sensitivity(results / 'sensitivity' if sensitivity_results is None else sensitivity_results)
    return d

def normalize_lead(frame):
    alias = {
        'd_pub_s': 'publication_delay_s',
        'd_act_s': 'action_delay_s',
        'lower_lead_s': 'finite_lower_s',
        'upper_lead_s': 'finite_upper_s',
        'lead_lower_s': 'finite_lower_s',
        'lead_upper_s': 'finite_upper_s',
        'includes_no_hit_by_horizon': 'includes_no_hit',
    }
    frame = frame.rename(columns={a: b for a, b in alias.items() if a in frame and b not in frame}).copy()
    need_columns(frame, ['regime', 'publication_delay_s', 'action_delay_s', 'finite_lower_s', 'finite_upper_s', 'includes_no_hit'], 'lead_margins.csv')
    if frame[['publication_delay_s', 'action_delay_s']].isna().any().any():
        raise ValueError('lead_margins.csv: missing delays')
    if (frame[['publication_delay_s', 'action_delay_s']] < 0).any().any():
        raise ValueError('lead_margins.csv: negative delays')
    if (frame.finite_lower_s > frame.finite_upper_s).any():
        raise ValueError('lead_margins.csv: inverted finite margin interval')
    return frame

def normalize_q4(frame):
    alias = {
        'target_arrival_status': 'status',
        'arrival_lower_s': 'finite_arrival_lower_s',
        'arrival_upper_s': 'finite_arrival_upper_s',
        'lead_lower_s': 'finite_lead_lower_s',
        'lead_upper_s': 'finite_lead_upper_s',
        'd_pub_s': 'publication_delay_s',
        'd_act_s': 'action_delay_s',
        'lower_lead_s': 'finite_lead_lower_s',
        'upper_lead_s': 'finite_lead_upper_s',
        'effective_lead_lower_s': 'finite_lead_lower_s',
        'effective_lead_upper_s': 'finite_lead_upper_s',
    }
    frame = frame.rename(columns={a: b for a, b in alias.items() if a in frame and b not in frame}).copy()
    need_columns(frame, ['regime', 'scenario_id', 'status', 'publication_delay_s', 'action_delay_s', 'finite_lead_lower_s', 'finite_lead_upper_s', 'travel_time_s', 'horizon_s', 'position_m', 'finite_arrival_lower_s', 'finite_arrival_upper_s', 'terrain_identity', 'stage_bridge', 'includes_no_hit', 'natural_warning_radius_m'], 'q4_spatial_demo.csv')
    known = ['arrived', 'valid_no_arrival_by_horizon', 'no_arrival_by_horizon', 'unknown', 'failed']
    if not frame.status.isin(known).all():
        raise ValueError('q4_spatial_demo.csv: unsupported propagation outcome')
    if (frame.horizon_s <= 0).any() or (frame.position_m <= 0).any():
        raise ValueError('q4_spatial_demo.csv: positive horizon and target distance required')
    if not np.allclose(frame.position_m, 1000, rtol=0, atol=1e-09) or not np.allclose(frame.horizon_s, 180, rtol=0, atol=1e-09):
        raise ValueError('q4_spatial_demo.csv: this fixed synthetic figure requires the declared 1000m target and 180s propagation horizon')
    if frame.natural_warning_radius_m.notna().any():
        raise ValueError('q4_spatial_demo.csv: synthetic target must not carry a natural warning radius')
    if frame[['terrain_identity', 'stage_bridge']].isna().any().any():
        raise ValueError('q4_spatial_demo.csv: explicit terrain and bridge identities required')
    if frame.loc[frame.status != 'arrived', ['finite_lead_lower_s', 'finite_lead_upper_s']].notna().any().any():
        raise ValueError('q4_spatial_demo.csv: no-arrival/unknown must not carry a finite target lead')
    if frame.loc[frame.status != 'arrived', ['travel_time_s', 'finite_arrival_lower_s', 'finite_arrival_upper_s']].notna().any().any():
        raise ValueError('q4_spatial_demo.csv: nonarrival must retain empty travel/arrival times')
    arrived = frame[frame.status == 'arrived']
    if ((arrived.travel_time_s <= 0) | (arrived.travel_time_s > arrived.horizon_s)).any():
        raise ValueError('q4_spatial_demo.csv: arrived travel time outside declared horizon')
    shift = arrived.publication_delay_s + arrived.action_delay_s
    if not np.allclose(arrived.finite_arrival_lower_s - shift, arrived.finite_lead_lower_s, atol=1e-08, rtol=0) or not np.allclose(arrived.finite_arrival_upper_s - shift, arrived.finite_lead_upper_s, atol=1e-08, rtol=0):
        raise ValueError('q4_spatial_demo.csv: target lead does not equal same-scenario arrival minus delays')
    if arrived[['finite_lead_lower_s', 'finite_lead_upper_s']].isna().any().any() or (arrived.finite_lead_lower_s > arrived.finite_lead_upper_s).any():
        raise ValueError('q4_spatial_demo.csv: invalid arrived finite interval')
    return frame

def normalize_validation(value):
    rows = value.get('claims', value.get('results', value.get('families', [])))
    if isinstance(rows, dict):
        flat = []
        for regime, v in rows.items():
            seq = v.get('families', v.get('claims', [])) if isinstance(v, dict) else v
            for row in seq:
                flat.append(dict(regime=regime, **row))
        rows = flat
    if not isinstance(rows, list):
        raise ValueError('validation.json: claims/results must contain rows')
    frame = pd.DataFrame(rows)
    if not len(frame):
        return frame
    aliases = {
        'condition': 'regime',
        'n_samples': 'n',
        'included': 'k',
        'empirical_coverage': 'empirical_mass',
        'one_sided_lower_mass': 'lower_bound',
        'one_sided_lower_bound': 'lower_bound',
        'one_sided_lower': 'lower_bound',
        'certified': 'passed',
        'all_target_pass': 'passed',
    }
    frame = frame.rename(columns={a: b for a, b in aliases.items() if a in frame and b not in frame})
    need_columns(frame, ['regime', 'family', 'n', 'k', 'empirical_mass', 'lower_bound', 'passed'], 'validation.json claims')
    if not all((isinstance(x, (bool, np.bool_)) for x in frame.passed)):
        raise ValueError('validation.json: passed must contain Boolean values, not strings')
    if not np.isfinite(frame[['n', 'k', 'empirical_mass', 'lower_bound']].to_numpy(float)).all():
        raise ValueError('validation.json: finite counts, proportions and bounds required')
    if ((frame.k < 0) | (frame.n <= 0) | (frame.k > frame.n)).any() or not np.equal(frame[['n', 'k']], np.floor(frame[['n', 'k']])).all().all():
        raise ValueError('validation.json: invalid integer binomial counts')
    target = 1 - float(value.get('alpha', 0.05))
    if not (frame.passed == (frame.lower_bound >= target - 1e-12)).all():
        raise ValueError('validation.json: pass status contradicts the exact lower bound')
    if not np.allclose(frame.k / frame.n, frame.empirical_mass, rtol=0, atol=1e-10):
        raise ValueError('validation empirical proportions do not match k/n')
    if ((frame.lower_bound < 0) | (frame.lower_bound > frame.empirical_mass + 1e-12)).any():
        raise ValueError('validation bounds outside valid range')
    return frame

def training_rows(data, regime, method):
    v = data['gp_training']
    if regime in v and method in v[regime]:
        return pd.DataFrame(v[regime][method])
    if f'{regime}_{method}' in v:
        return pd.DataFrame(v[f'{regime}_{method}'])
    if isinstance(v.get('rows'), list):
        f = pd.DataFrame(v['rows'])
        return f[(f.regime == regime) & (f.method == method)]
    raise ValueError(f'gp_training.json: missing {regime}/{method} training rows')

def base_condition(data, regime):
    conditions = data['baseline_results']['conditions']
    return conditions[BASE_KEYS[regime]] if BASE_KEYS[regime] in conditions else conditions[regime]

def family_cdf(x, family):
    u = np.clip((np.asarray(x) + 1.5) / 3, 0, 1)
    if family == 'uniform_support':
        return u
    if family == 'lower_mode_triangular':
        return 2 * u - u * u
    if family == 'nominal_triangular':
        return np.where(u <= 0.5, 2 * u * u, 1 - 2 * (1 - u) ** 2)
    raise ValueError(f'Unsupported declared family {family}')

def source_cdf(frame, edges, family='nominal_triangular'):
    f = frame.sort_values('d')
    x = f.d.to_numpy()
    expected = -1.5 + (np.arange(len(f)) + 0.5) * 3 / len(f)
    if not np.allclose(x, expected, rtol=0, atol=1e-09):
        raise ValueError('CDF diagnostic requires the saved equal-cell midpoint input grid')
    weights = np.diff(family_cdf(np.linspace(-1.5, 1.5, len(f) + 1), family))
    hit = f.status.to_numpy() == 'hit'
    t = f.tau_s.to_numpy()
    c = np.array([weights[hit & (t <= e)].sum() for e in edges])
    return (c, float(weights[f.status.to_numpy() == 'no_hit'].sum()), float(weights[f.status.to_numpy() == 'unknown'].sum()))

def spec_of(c, name):
    s = c[name]
    return s.get('set', s) if isinstance(s, dict) else s

def fmt(x, digits=6):
    if x is None or (isinstance(x, (float, np.floating)) and (not np.isfinite(x))):
        return '—'
    if isinstance(x, (bool, np.bool_)):
        return '是' if x else '否'
    if digits == 0:
        return f'{float(x):.0f}'
    if abs(float(x)) < 0.001 and float(x) != 0:
        return f'{float(x):.3g}'
    return f'{float(x):.{digits}f}'.rstrip('0').rstrip('.')

def interval_text(spec, digits=3):
    if spec is None:
        return '无可行集合'
    interval = spec.get('finite_interval_s')
    text = '∅' if interval is None else f'[{fmt(interval[0], digits)}, {fmt(interval[1], digits)}]'
    if spec.get('includes_no_hit_by_horizon', False):
        text += '∪{⊥}'
    return text

def response_figure(data, out, audit):
    fig, axes = new_figure('response')
    C = STYLE['colors']
    for col, regime in enumerate(REGIMES):
        truth = data['truth'][data['truth'].regime == regime].sort_values('d')
        h = truth.status == 'hit'
        ax = axes[0, col]
        ax.plot(truth.d[h], truth.tau_s[h], color=C['truth'], lw=STYLE['line_width'])
        for method, mark in STYLE['query_markers'].items():
            tr = training_rows(data, regime, method)
            h = tr.status == 'hit'
            ax.scatter(tr.d[h], tr.tau_s[h], s=STYLE['scatter']['query_area'], color=C[method], marker=mark, linewidths=STYLE['scatter']['query_edge_width'], zorder=STYLE['scatter']['query_zorder'])
            known = tr.status.isin(['hit', 'no_hit'])
            offset = STYLE['scatter']['class_display_jitter'] if method == 'active' else -STYLE['scatter']['class_display_jitter']
            axes[1, col].scatter(tr.d[known], (tr.status[known] == 'hit').astype(float) + offset, s=STYLE['scatter']['query_area'], color=C[method], marker=mark, linewidths=STYLE['scatter']['query_edge_width'])
        ax.set(xlabel='Coherent velocity offset (m/s)', ylabel='Elapsed first-hit time (s)', xlim=STYLE['views']['offset_xlim'])
        panel(ax, chr(97 + col), NAMES[regime])
        f = data['pred']
        p = f[(f.regime == regime) & (f.method == 'active') & (f.budget == f.budget.max()) & (f['mode'] == 'raw')].sort_values('d')
        axes[1, col].plot(p.d, p.p_hit, color=C['active'], lw=STYLE['line_widths']['score'])
        axes[1, col].axhline(0.5, color=C['class_threshold'], ls=STYLE['line_styles']['baseline'], lw=STYLE['line_widths']['class_threshold'])
        axes[1, col].set(xlabel='Coherent velocity offset (m/s)', ylabel='Hit class / GP score', xlim=STYLE['views']['offset_xlim'], ylim=STYLE['views']['class_ylim'], yticks=STYLE['views']['class_yticks'])
        panel(axes[1, col], chr(99 + col))
    fig.legend(handles=[Line2D([], [], color=C['truth'], label='Direct Q2'), Line2D([], [], color=C['passive'], marker=STYLE['query_markers']['passive'], ls=STYLE['legend_marker_line'], label='Passive queries'), Line2D([], [], color=C['active'], marker=STYLE['query_markers']['active'], ls=STYLE['legend_marker_line'], label='Active queries')], ncol=STYLE['figures']['response']['legend_columns'], bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig01_mapping_queries', out, audit)

def learning_figure(data, out, audit):
    fig, axes = new_figure('learning')
    C = STYLE['colors']
    m = data['metrics']
    for col, regime in enumerate(REGIMES):
        for method in ['passive', 'active']:
            for mode, appearance in STYLE['mode_styles'].items():
                ls = appearance['line']
                mark = appearance['marker']
                f = m[(m.regime == regime) & (m.method == method) & (m['mode'] == mode) & (m.subset == 'all')].sort_values('budget')
                axes[0, col].plot(f.budget, f.mae_s, color=C[method], ls=ls, marker=mark, ms=STYLE['marker_size'])
                axes[1, col].plot(f.total_nominal_q2_calls, f.false_negative, color=C[method], ls=ls, marker=mark, ms=STYLE['marker_size'])
        b = data['baseline_test_metrics'].get(regime, {})
        if b.get('separated_002', {}).get('mae_s') is not None:
            axes[0, col].axhline(b['separated_002']['mae_s'], color=C['baseline'], ls=STYLE['line_styles']['baseline'], lw=STYLE['line_widths']['baseline_reference'])
        axes[0, col].set(xlabel='Training Q2 calls', ylabel='Elapsed-time MAE (s)', xticks=sorted(m.budget.unique()))
        axes[1, col].set(xlabel='Training + fallback Q2 calls', ylabel='Missed true hits (count)', ylim=(STYLE['views']['miss_bottom'], max(STYLE['views']['miss_minimum_top'], m[(m.regime == regime) & (m.subset == 'all')].false_negative.max() + STYLE['views']['miss_headroom'])))
        axes[1, col].yaxis.set_major_locator(MaxNLocator(integer=True))
        for row in [0, 1]:
            panel(axes[row, col], chr(97 + row * 2 + col), NAMES[regime] if row == 0 else None)
    fig.legend(handles=[Line2D([], [], color=C['passive'], label='Passive GP'), Line2D([], [], color=C['active'], label='Active GP'), Line2D([], [], color=C['fallback_legend'], ls=STYLE['mode_styles']['fallback']['line'], label='With fallback'), Line2D([], [], color=C['baseline'], ls=STYLE['line_styles']['baseline'], label='Separated baseline MAE')], ncol=STYLE['figures']['learning']['legend_columns'], bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig02_gp_accuracy_cost', out, audit)

def validation_figure(data, out, audit):
    fig, axes = new_figure('validation')
    v = data['validation_rows']
    C = STYLE['colors']
    for col, regime in enumerate(REGIMES):
        ax = axes[col]
        ax.axhline(0.95, color=C['failure'], lw=STYLE['line_widths']['validation_target'], ls=STYLE['mode_styles']['fallback']['line'])
        f = v[v.regime == regime].copy() if len(v) else v
        if not len(f):
            ax.text(*STYLE['pending_text_anchor'], 'Random validation pending', ha=STYLE['status_horizontal_alignment'], va=STYLE['status_vertical_alignment'], transform=ax.transAxes)
        else:
            f = f.set_index('family').reindex(FAMILIES)
            if f.empirical_mass.isna().any():
                raise ValueError(f'validation.json: three declared families required for {regime}')
            x = np.arange(3)
            ax.plot(x, f.empirical_mass, color=C['truth'], marker=STYLE['mode_styles']['raw']['marker'], ls=STYLE['mode_styles']['raw']['line'], ms=STYLE['mark_sizes']['validation'])
            ax.plot(x, f.lower_bound, color=C['baseline'], marker=STYLE['mode_styles']['fallback']['marker'], ls=STYLE['mode_styles']['fallback']['line'], ms=STYLE['mark_sizes']['validation'])
            for j, (_, row) in enumerate(f.iterrows()):
                ax.plot([j, j], [row.lower_bound, row.empirical_mass], color=C['bound_connector'], lw=STYLE['line_widths']['bound_connector'])
        ax.set(xticks=range(3), xticklabels=FAMILY_LABELS, xlim=STYLE['views']['validation_xlim'], ylim=(min(STYLE['views']['validation_lower_floor'], float(f.lower_bound.min()) - STYLE['views']['validation_lower_padding']) if len(f) else STYLE['views']['validation_lower_floor'], STYLE['views']['validation_upper']), ylabel='Frozen-set probability mass')
        panel(ax, chr(97 + col), NAMES[regime])
    handles = [Line2D([], [], color=C['failure'], ls=STYLE['mode_styles']['fallback']['line'], label='Target 0.95')]
    if len(v):
        handles = [Line2D([], [], color=C['truth'], marker=STYLE['mode_styles']['raw']['marker'], label='New random sample mass'), Line2D([], [], color=C['baseline'], marker=STYLE['mode_styles']['fallback']['marker'], ls=STYLE['mode_styles']['fallback']['line'], label='Simultaneous exact lower bound')] + handles
    fig.legend(handles=handles, ncol=STYLE['figures']['validation']['legend_columns'] or len(handles), bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig03_frozen_validation', out, audit)

def distribution_figure(data, out, audit):
    fig, axes = new_figure('distribution')
    C = STYLE['colors']
    edges = np.asarray(data['baseline_results']['edges_s'])
    for col, regime in enumerate(REGIMES):
        ax = axes[0, col]
        truth = data['truth'][data['truth'].regime == regime]
        c, _, _ = source_cdf(truth, edges)
        ax.plot(edges, c, color=C['truth'], lw=STYLE['line_widths']['distribution_truth'])
        bc = base_condition(data, regime)
        pr = bc['propagation'][0]
        cdfv = np.asarray(pr.get('cdf', np.r_[0, np.cumsum(pr['finite_mass'])]))
        ax.plot(edges, cdfv, color=C['baseline'], ls=STYLE['line_styles']['baseline'], lw=STYLE['line_widths']['distribution_baseline'])
        p = data['pred']
        maxbudget = p.budget.max()
        for method in ['passive', 'active']:
            f = p[(p.regime == regime) & (p.method == method) & (p.budget == maxbudget) & (p['mode'] == 'raw')]
            c, _, _ = source_cdf(f, edges)
            ax.plot(edges, c, color=C[method], ls=STYLE['line_styles'][method], lw=STYLE['line_widths']['distribution_gp'])
        ax.set(xlim=STYLE['views']['cdf_xlim'][regime], ylim=STYLE['views']['cdf_ylim'], xlabel='Elapsed first-hit time (s)', ylabel='Unconditional finite CDF')
        panel(ax, chr(97 + col), NAMES[regime])
        ax = axes[1, col]
        for i, (key, label) in enumerate([('nominal_central', 'Nominal central'), ('nominal_shortest', 'Nominal shortest'), ('robust_shortest', 'Robust shortest')]):
            s = spec_of(bc, key)
            finite = s.get('finite_interval_s') if s else None
            if finite is not None:
                ax.plot(finite, [i, i], color=C['baseline'], lw=STYLE['line_widths']['interval'], marker=STYLE['interval_marker'], ms=STYLE['mark_sizes']['interval_endpoint'])
                if s.get('includes_no_hit_by_horizon'):
                    ax.scatter(finite[1], i, s=STYLE['scatter']['nohit_area'], marker=STYLE['nohit_marker'], facecolors=STYLE['hollow_marker_facecolor'], edgecolors=C['baseline'])
        ax.set(yticks=[0, 1, 2], yticklabels=['Nominal central', 'Nominal shortest', 'Robust shortest'], xlabel='Elapsed interval endpoints (s)')
        ax.invert_yaxis()
        panel(ax, chr(99 + col))
    fig.legend(handles=[Line2D([], [], color=C['truth'], label='Q2 midpoint diagnostic'), Line2D([], [], color=C['baseline'], ls=STYLE['line_styles']['baseline'], label='Integrated 1D mapping'), Line2D([], [], color=C['passive'], ls=STYLE['line_styles']['passive'], label='Passive GP raw'), Line2D([], [], color=C['active'], label='Active GP raw')], ncol=STYLE['figures']['distribution']['legend_columns'], bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig04_probability_sets', out, audit)

def lead_figure(data, out, audit):
    lead = data['lead']
    actions = sorted(lead.action_delay_s.unique())
    C = STYLE['colors']
    fig, axes = new_figure('lead')
    marks = STYLE['lead_markers']
    styles = STYLE['lead_line_styles']
    for col, regime in enumerate(REGIMES):
        ax = axes[col]
        f = lead[lead.regime == regime]
        for j, action in enumerate(actions):
            z = f[f.action_delay_s == action].sort_values('publication_delay_s')
            if not len(z):
                continue
            color = C[STYLE['lead_color_order'][j % len(STYLE['lead_color_order'])]]
            center = (z.finite_lower_s + z.finite_upper_s) / 2
            err = (z.finite_upper_s - z.finite_lower_s) / 2
            ax.errorbar(z.publication_delay_s, center, yerr=err, color=color, ls=styles[j % len(styles)], marker=marks[j % len(marks)], ms=STYLE['mark_sizes']['lead'], capsize=STYLE['errorbar_capsize'], lw=STYLE['line_widths']['lead'], label=f'Action delay {fmt(action, 1)} s')
        ax.axhline(0, color=C['failure'], ls=STYLE['line_styles']['baseline'], lw=STYLE['line_widths']['local_zero'])
        ax.set(xlabel='Publication elapsed time (s)', ylabel='Finite local lead interval (s)')
        panel(ax, chr(97 + col), NAMES[regime] + ('; + no-hit' if f.includes_no_hit.astype(bool).any() else ''))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=min(STYLE['figures']['lead']['legend_max_columns'], len(actions)), bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig05_local_lead', out, audit)

def q4_figure(data, out, audit):
    demo = data['q4_demo']
    if demo is None:
        return False
    pairs = demo[['publication_delay_s', 'action_delay_s']].drop_duplicates().sort_values(['publication_delay_s', 'action_delay_s'])
    pub = float(pairs.iloc[-1].publication_delay_s)
    act = float(pairs.iloc[-1].action_delay_s)
    fig, axes = new_figure('q4')
    C = STYLE['colors']
    for col, regime in enumerate(REGIMES):
        f = demo[(demo.regime == regime) & (demo.publication_delay_s == pub) & (demo.action_delay_s == act)].sort_values('scenario_id')
        if not len(f):
            raise ValueError('q4_spatial_demo.csv: selected saved delay pair missing in a regime')
        ax = axes[col]
        labels = []
        for j, (_, r) in enumerate(f.iterrows()):
            labels.append(f'S{j + 1}')
            if r.status == 'arrived':
                ax.plot([r.finite_lead_lower_s, r.finite_lead_upper_s], [j, j], color=C['active'], lw=STYLE['line_widths']['interval'], marker=STYLE['interval_marker'], ms=STYLE['mark_sizes']['interval_endpoint'])
            else:
                label = f'No arrival by {fmt(r.horizon_s, 1)} s' if 'no_arrival' in r.status else str(r.status)
                ax.text(STYLE['q4_status_x_fraction'], j, label, transform=ax.get_yaxis_transform(), va=STYLE['status_vertical_alignment'], ha=STYLE['status_horizontal_alignment'], fontsize=STYLE['status_size'], color=C['passive'])
        ax.axvline(0, color=C['failure'], ls=STYLE['line_styles']['baseline'], lw=STYLE['line_widths']['q4_zero'])
        ax.set(yticks=range(len(f)), yticklabels=labels, xlabel='Conditional target lead interval (s)', ylabel='Fixed synthetic Q4 scenario')
        ax.invert_yaxis()
        panel(ax, chr(97 + col), NAMES[regime])
    fig.suptitle(f'Synthetic 8-degree slope; 1000 m; elapsed publication {fmt(pub, 1)} s, action {fmt(act, 1)} s', y=STYLE['supertitle_y'], fontsize=STYLE['supertitle_size'])
    export(fig, 'fig06_synthetic_target_lead', out, audit)
    return True
SENSITIVITY_FILES = ['phase_map.csv', 'chi_perturbation.csv', 'crossing_amplification.csv', 'local_derivatives.csv', 'derivative_step_stability.csv', 'distribution_sensitivity.csv', 'set_sensitivity.csv', 'summary.json']

def read_sensitivity(folder):
    for name in SENSITIVITY_FILES:
        if not (folder / name).is_file():
            raise ValueError('Missing required saved sensitivity result: sensitivity/' + name)
    d = {name[:-4]: pd.read_csv(folder / name) for name in SENSITIVITY_FILES if name.endswith('.csv')}
    d['summary'] = load_json(folder / 'summary.json')
    d['legacy_numerical_convergence'] = pd.read_csv(folder / 'legacy_numerical_convergence.csv') if (folder / 'legacy_numerical_convergence.csv').is_file() else None
    fields = {
        'phase_map': ['regime', 'q0', 'offset_m_s', 'status', 'tau_s'],
        'chi_perturbation': ['regime', 'chi_multiplier', 'offset_m_s', 'status', 'tau_s'],
        'crossing_amplification': ['regime', 'status', 'tau_s', 'qdot_cross_s_inv', 'seconds_per_0_001_q', 'horizon_slack_s', 'shear_qdot_s_inv', 'drainage_qdot_s_inv', 'drainage_to_abs_shear_ratio'],
        'local_derivatives': ['regime', 'parameter', 'step', 'absolute_parameter_step', 'minus_status', 'plus_status', 'minus_tau_s', 'plus_tau_s', 'derivative', 'derivative_unit', 'elasticity', 'parameter_reference_scale', 'scaled_sensitivity'],
        'derivative_step_stability': ['regime', 'parameter', 'finite_step_count', 'total_step_count'],
        'distribution_sensitivity': ['regime', 'support_halfwidth_m_s', 'family', 'target_mass', 'hit_mass', 'no_hit_mass', 'unknown_mass', 'total_mass'],
        'set_sensitivity': ['regime', 'support_halfwidth_m_s', 'target_mass', 'method', 'finite_lower_s', 'finite_upper_s', 'finite_width_s', 'includes_no_hit', 'worst_model_mass'],
    }
    for name, cols in fields.items():
        need_columns(d[name], cols, 'sensitivity/' + name + '.csv')
    for name in ['phase_map', 'chi_perturbation']:
        frame = d[name]
        if not frame.regime.isin(REGIMES).all() or not frame.status.isin(['hit', 'no_hit', 'unknown', 'initial']).all():
            raise ValueError('invalid sensitivity condition/status: ' + name)
        if frame.loc[frame.status != 'hit', 'tau_s'].notna().any():
            raise ValueError('non-hit sensitivity members must not carry finite time')
        hit_times = frame.loc[frame.status == 'hit', 'tau_s'].to_numpy(float)
        if not np.isfinite(hit_times).all() or not ((hit_times > 0) & (hit_times <= 73)).all():
            raise ValueError('hit sensitivity members require finite time in(0,73]')
    phase = d['phase_map']
    if phase.duplicated(['regime', 'q0', 'offset_m_s']).any():
        raise ValueError('duplicate sensitivity phase design point')
    for regime in REGIMES:
        frame = phase[phase.regime == regime]
        if len(frame) != frame.q0.nunique() * frame.offset_m_s.nunique():
            raise ValueError('incomplete phase design grid')
    amp = d['crossing_amplification']
    if len(amp) != len(REGIMES) or set(amp.regime) != set(REGIMES):
        raise ValueError('one nominal amplification row per condition required')
    if not np.allclose(amp.seconds_per_0_001_q, 0.001 / abs(amp.qdot_cross_s_inv), rtol=1e-09, atol=1e-12):
        raise ValueError('local amplification does not match saved crossing slope')
    if not np.allclose(amp.shear_qdot_s_inv + amp.drainage_qdot_s_inv, amp.qdot_cross_s_inv, rtol=1e-08, atol=1e-12):
        raise ValueError('saved mechanical terms do not sum to qdot')
    dist = d['distribution_sensitivity']
    if not np.isfinite(dist[['hit_mass', 'no_hit_mass', 'unknown_mass', 'total_mass']]).all().all() or (dist[['hit_mass', 'no_hit_mass', 'unknown_mass']] < -1e-10).any().any():
        raise ValueError('invalid sensitivity probability mass')
    if not np.allclose(dist.total_mass, dist.hit_mass + dist.no_hit_mass + dist.unknown_mass, atol=1e-08, rtol=0):
        raise ValueError('saved sensitivity total_mass contradicts component masses')
    if not np.allclose(dist.hit_mass + dist.no_hit_mass + dist.unknown_mass, 1, atol=1e-08, rtol=0):
        raise ValueError('sensitivity probability mass does not close')
    local = d['local_derivatives']
    if not local.parameter.isin(['offset_m_s', 'q0', 'chi_relative']).all():
        raise ValueError('unsupported sensitivity derivative parameter')
    if not np.isfinite(local.absolute_parameter_step.to_numpy(float)).all() or not (local.absolute_parameter_step > 0).all():
        raise ValueError('finite-difference absolute step must be positive finite')
    for side in ['minus', 'plus']:
        status = local[side + '_status']
        tau = local[side + '_tau_s']
        if not status.isin(['hit', 'no_hit', 'initial', 'unknown']).all():
            raise ValueError('invalid finite-difference member status')
        times = tau[status == 'hit'].to_numpy(float)
        if not np.isfinite(times).all() or not ((times > 0) & (times <= 73)).all():
            raise ValueError('finite-difference hit side needs finite time in(0,73]')
        if tau[status != 'hit'].notna().any():
            raise ValueError('finite-difference non-hit side cannot carry time')
    both = local.minus_status.eq('hit') & local.plus_status.eq('hit')
    if not np.isfinite(local.loc[both, 'derivative'].to_numpy(float)).all():
        raise ValueError('both-hit finite-difference pair requires its finite derivative, not category-transition NaN')
    expected = (local.loc[both, 'plus_tau_s'] - local.loc[both, 'minus_tau_s']) / (2 * local.loc[both, 'absolute_parameter_step'])
    if not np.allclose(local.loc[both, 'derivative'], expected, rtol=1e-07, atol=1e-08):
        raise ValueError('finite difference disagrees with saved two-side times and step')
    if local.loc[~both, ['derivative', 'elasticity', 'scaled_sensitivity']].notna().any().any():
        raise ValueError('category-transition pair must not have a finite derivative or normalized derivative')
    nominal = d['crossing_amplification'].set_index('regime')
    for index, row in local[both].iterrows():
        scale = 1.5 if row.parameter == 'offset_m_s' else 0.5 if row.parameter == 'q0' else float(d['summary']['conditions'][row.regime])
        if not math.isfinite(row.parameter_reference_scale) or not math.isclose(row.parameter_reference_scale, scale, rel_tol=1e-10, abs_tol=1e-14):
            raise ValueError('finite-difference reference scale differs from declared design')
        scaled = row.derivative * scale / float(nominal.loc[row.regime, 'tau_s'])
        if not math.isfinite(row.scaled_sensitivity) or not math.isclose(row.scaled_sensitivity, scaled, rel_tol=1e-07, abs_tol=1e-08):
            raise ValueError('scaled sensitivity disagrees with derivative, reference scale and nominal time')
        if row.parameter == 'offset_m_s':
            if pd.notna(row.elasticity):
                raise ValueError('ordinary elasticity at nominal d=0 must remain undefined')
        elif not math.isfinite(row.elasticity) or not math.isclose(row.elasticity, scaled, rel_tol=1e-07, abs_tol=1e-08):
            raise ValueError('ordinary elasticity disagrees with finite difference')
    sets = d['set_sensitivity']

    def strict_bool(value):
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, str) and value in ['True', 'False']:
            return value == 'True'
        raise ValueError('sensitivity set includes_no_hit must contain Boolean True/False, not arbitrary strings')
    sets['includes_no_hit'] = sets.includes_no_hit.map(strict_bool)
    if not np.isfinite(sets[['finite_lower_s', 'finite_upper_s', 'finite_width_s', 'worst_model_mass', 'target_mass']].to_numpy(float)).all():
        raise ValueError('finite sensitivity candidate requires finite endpoints and masses')
    if (sets.finite_width_s < 0).any() or not np.allclose(sets.finite_upper_s - sets.finite_lower_s, sets.finite_width_s, atol=1e-09, rtol=0):
        raise ValueError('sensitivity set endpoint/width mismatch')
    if (sets.worst_model_mass < sets.target_mass - 1e-10).any():
        raise ValueError('sensitivity set violates its declared numerical mass target')
    legacy = d['legacy_numerical_convergence']
    if legacy is not None:
        need_columns(legacy, ['base_status', 'tight_status', 'base_tau_s', 'tight_tau_s', 'tau_abs_difference_s', 'status_agrees', 'base_rtol', 'tight_rtol', 'base_max_step_s', 'tight_max_step_s'], 'sensitivity/legacy_numerical_convergence.csv')
        if len(legacy) != 6:
            raise ValueError('historical numerical convergence requires exactly six saved cases')
        if not legacy.status_agrees.map(lambda value: isinstance(value, (bool, np.bool_)) and bool(value)).all():
            raise ValueError('historical convergence status_agrees must be Boolean True for every case')
        if not legacy.base_status.isin(['hit', 'no_hit']).all() or not legacy.base_status.eq(legacy.tight_status).all():
            raise ValueError('historical convergence base/tight statuses must agree')
        for base, tight, ratio in [('base_rtol', 'tight_rtol', 10.0), ('base_max_step_s', 'tight_max_step_s', 2.0)]:
            values = legacy[[base, tight]].to_numpy(float)
            if not np.isfinite(values).all() or not (values > 0).all() or (not np.allclose(legacy[base] / legacy[tight], ratio, rtol=1e-12, atol=0)):
                raise ValueError('historical numerical tightening ratios contradict the declared 10x tolerance / 2x step design')
        hit = legacy.base_status.eq('hit')
        columns = ['base_tau_s', 'tight_tau_s', 'tau_abs_difference_s']
        if legacy.loc[~hit, columns].notna().any().any():
            raise ValueError('historical no-hit convergence case must leave all time fields empty')
        values = legacy.loc[hit, columns].to_numpy(float)
        if not np.isfinite(values).all() or not ((values[:, :2] > 0) & (values[:, :2] <= 73)).all() or (not (values[:, 2] >= 0).all()):
            raise ValueError('historical hit convergence cases require valid finite times and nonnegative differences')
        expected = abs(legacy.loc[hit, 'base_tau_s'] - legacy.loc[hit, 'tight_tau_s'])
        if not np.allclose(legacy.loc[hit, 'tau_abs_difference_s'], expected, rtol=1e-07, atol=1e-12):
            raise ValueError('historical convergence difference contradicts saved finite times')
    return d

def grid_cell_edges(x):
    x = np.asarray(x, float)
    if len(x) < 2 or not np.all(np.diff(x) > 0):
        raise ValueError('ordered design coordinates required')
    return np.r_[x[0] - (x[1] - x[0]) / 2, (x[:-1] + x[1:]) / 2, x[-1] + (x[-1] - x[-2]) / 2]

def sensitivity_phase_figure(data, out, audit):
    from matplotlib.colors import Normalize
    from matplotlib.patches import Patch
    d = data['sensitivity']['phase_map']
    st = STYLE['sensitivity']
    fig, axes = new_figure('sensitivity_phase')
    for col, regime in enumerate(REGIMES):
        f = d[d.regime == regime]
        q = np.sort(f.q0.unique())
        x = np.sort(f.offset_m_s.unique())
        pivot = f.pivot(index='q0', columns='offset_m_s', values='tau_s').reindex(index=q, columns=x)
        vals = pivot.to_numpy(float)
        mask = np.ma.masked_invalid(vals)
        cmap = plt.get_cmap(st['phase_cmap']).copy()
        cmap.set_bad(st['no_hit_color'])
        maximum = float(np.nanmax(vals))
        mesh = axes[col].pcolormesh(grid_cell_edges(x), grid_cell_edges(q), mask, cmap=cmap, norm=Normalize(0, maximum), shading=st['phase_shading'], rasterized=st['phase_rasterized'])
        cb = fig.colorbar(mesh, ax=axes[col], fraction=st['colorbar_fraction'], pad=st['colorbar_pad'])
        cb.ax._q3_auxiliary = True
        cb.set_label('Elapsed first-hit time (s)')
        axes[col].scatter([0], [0.5], marker=st['nominal_marker'], s=st['nominal_area'], facecolors=st['nominal_face'], edgecolors=STYLE['colors']['truth'], linewidths=st['nominal_edge_width'], zorder=st['nominal_zorder'])
        for status, marker in [('unknown', st['unknown_marker']), ('initial', st['initial_marker'])]:
            group = f[f.status == status]
            if len(group):
                axes[col].scatter(group.offset_m_s, group.q0, marker=marker, s=st['category_area'], color=STYLE['colors']['failure'], linewidths=st['category_edge_width'])
        axes[col].set(xlabel='Whole-path velocity offset (m/s)', ylabel='Initial stress ratio q0', xlim=(x[0], x[-1]), ylim=(q[0], q[-1]))
        panel(axes[col], chr(97 + col), NAMES[regime])
    handles = [Patch(facecolor=st['no_hit_color'], label='No hit within 73 s'), Line2D([], [], marker=st['nominal_marker'], ls=STYLE['legend_marker_line'], markerfacecolor=st['nominal_face'], markeredgecolor=STYLE['colors']['truth'], label='Nominal q0=0.5, offset=0')]
    fig.legend(handles=handles, ncol=STYLE['figures']['sensitivity_phase']['legend_columns'], bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    export(fig, 'fig07_initial_velocity_sensitivity', out, audit)

def sensitivity_limits_figure(data, out, audit):
    sens = data['sensitivity']
    amp = sens['crossing_amplification'].set_index('regime')
    sets = sens['set_sensitivity']
    chi = sens['chi_perturbation']
    st = STYLE['sensitivity']
    C = STYLE['colors']
    fig, axes = new_figure('sensitivity_limits')
    values = [float(amp.loc[r, 'seconds_per_0_001_q']) for r in REGIMES]
    axes[0, 0].scatter(range(len(REGIMES)), values, color=[C['truth'], C['active']], s=st['amplification_area'], marker=st['amplification_marker'])
    axes[0, 0].set_ylim(min(values) / st['amplification_bottom_factor'], max(values) * st['amplification_top_factor'])
    axes[0, 0].text(*st['amplification_ratio_anchor'], 'Near / strong = ' + fmt(values[1] / values[0], 1), transform=axes[0, 0].transAxes, fontsize=STYLE['status_size'], va=STYLE['status_vertical_alignment'], ha=st['amplification_ratio_ha'])
    axes[0, 0].set(xticks=range(len(REGIMES)), xticklabels=[NAMES[r] for r in REGIMES], ylabel='Time / 0.001 state margin (s)', yscale='log')
    axes[0, 0].tick_params(axis='x', labelrotation=st['regime_tick_rotation'])
    panel(axes[0, 0], 'a', 'Local transversal amplification')
    ax = axes[0, 1]
    f = chi[(chi.regime == 'near') & np.isclose(chi.offset_m_s, 0)]
    order = f.sort_values('chi_multiplier')
    hit = order[order.status == 'hit']
    no = order[order.status == 'no_hit']
    ax.plot(order.chi_multiplier, order.tau_s, color=C['active'], ls=STYLE['line_styles']['active'], marker=STYLE['mode_styles']['raw']['marker'], lw=STYLE['line_widths']['distribution_gp'], ms=STYLE['marker_size'])
    if len(no):
        ax.scatter(no.chi_multiplier, np.full(len(no), st['nohit_strip_y']), transform=ax.get_xaxis_transform(), marker=STYLE['nohit_marker'], s=st['nohit_strip_area'], facecolors=STYLE['hollow_marker_facecolor'], edgecolors=C['passive'], clip_on=False, label='No hit (category)')
    ax.set(xlabel='χ / nominal χ (near-boundary condition)', ylabel='Elapsed first-hit time (s)', xlim=(float(order.chi_multiplier.min()) - st['chi_xlim_padding'], float(order.chi_multiplier.max()) + st['chi_xlim_padding']))
    panel(ax, 'b', 'Conditional χ perturbation')
    if len(no):
        ax.legend(loc=st['chi_legend_loc'], frameon=STYLE['legend']['frameon'], bbox_to_anchor=st['chi_legend_anchor'])
    for regime in REGIMES:
        a = sets[(sets.regime == regime) & (sets.method == 'robust_shortest') & np.isclose(sets.target_mass, 0.95)].sort_values('support_halfwidth_m_s')
        a = a.drop_duplicates('support_halfwidth_m_s')
        color = C['truth'] if regime == 'strong' else C['active']
        marker = st['condition_markers'][regime]
        line = st['condition_lines'][regime]
        axes[1, 0].plot(a.support_halfwidth_m_s, a.finite_width_s, color=color, ls=line, marker=marker, lw=STYLE['line_width'], ms=STYLE['marker_size'], label=NAMES[regime])
        b = sets[(sets.regime == regime) & (sets.method == 'robust_shortest') & np.isclose(sets.support_halfwidth_m_s, 1.5)].sort_values('target_mass')
        axes[1, 1].plot(b.target_mass, b.finite_width_s, color=color, ls=line, marker=marker, lw=STYLE['line_width'], ms=STYLE['marker_size'], label=NAMES[regime])
    axes[1, 0].set(xlabel='Offset support half-width (m/s)', ylabel='Robust finite width (s)', yscale='log')
    panel(axes[1, 0], 'c', 'Input-support sensitivity')
    axes[1, 1].set(xlabel='Target model mass', ylabel='Robust finite width (s)', yscale='log')
    panel(axes[1, 1], 'd', 'Coverage-target sensitivity')
    handles, labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=STYLE['figures']['sensitivity_limits']['legend_columns'], bbox_to_anchor=STYLE['legend_anchor'], **STYLE['legend'])
    for ax in [axes[0, 0], axes[1, 0], axes[1, 1]]:
        ax.tick_params(axis='y', labelsize=st['log_tick_size'])
    export(fig, 'fig08_boundary_probability_sensitivity', out, audit)

def sol():
    ap = argparse.ArgumentParser(description='Q3 plots')
    base = Path(__file__).resolve().parent
    ap.add_argument('--results', type=Path, default=base / 'results')
    ap.add_argument('--output', type=Path, default=base / 'figures_reproduced')
    ap.add_argument('--sensitivity-results', type=Path, default=None)
    args = ap.parse_args()
    results = args.results.resolve()
    out = args.output.resolve()
    if results == out or out.is_relative_to(results):
        raise ValueError('Figure output must not be inside saved analysis results')
    if out.exists():
        raise FileExistsError('output directory exists: use a new output directory')
    data = read_inputs(results, None if args.sensitivity_results is None else args.sensitivity_results.resolve())
    out.mkdir(parents=True, exist_ok=True)
    audit = out / 'audit'
    audit.mkdir(exist_ok=True)
    apply_style()
    response_figure(data, out, audit)
    learning_figure(data, out, audit)
    validation_figure(data, out, audit)
    distribution_figure(data, out, audit)
    lead_figure(data, out, audit)
    has_q4 = q4_figure(data, out, audit)
    sensitivity_phase_figure(data, out, audit)
    sensitivity_limits_figure(data, out, audit)
    save(audit / 'PLOT_RECEIPT.json', dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), source_sha256=sha(__file__), result_version=data['summary'].get('version', data['summary'].get('run_id')), input_sha256={f: sha(results / f) for f in REQUIRED + (['q4_spatial_demo.csv'] if has_q4 else [])}, sensitivity_input_sha256={name: sha((results / 'sensitivity' if args.sensitivity_results is None else args.sensitivity_results.resolve()) / name) for name in SENSITIVITY_FILES + (['legacy_numerical_convergence.csv'] if ((results / 'sensitivity' if args.sensitivity_results is None else args.sensitivity_results.resolve()) / 'legacy_numerical_convergence.csv').is_file() else [])}, output_sha256={p.name: sha(p) for p in out.iterdir() if p.is_file()}, q2_calls=0, training_calls=0, random_draws=0, styles=STYLE, audit_scope='self-contained rectangular axes geometry only; external PDF and manual visual audits recorded separately', independent_scientific_acceptance='pending'))
    print(json.dumps(dict(figure_count=8 if has_q4 else 7, q2_calls=0), ensure_ascii=False))
if __name__ == '__main__':
    sol()

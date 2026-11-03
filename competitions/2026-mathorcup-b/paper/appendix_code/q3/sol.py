from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_name, '1')
import sys
sys.dont_write_bytecode = True
import argparse, hashlib, json, math, platform, time, warnings
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import dataclass, asdict
from decimal import Decimal
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
import scipy, sklearn
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import beta, binom
from sklearn.gaussian_process import GaussianProcessClassifier, GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern
from sklearn.exceptions import ConvergenceWarning
PACKAGE_ROOT = Path(__file__).resolve().parent
ROOT = PACKAGE_ROOT
DATA = PACKAGE_ROOT / 'data'
OUT = PACKAGE_ROOT / 'results'
ATOL = (1e-11, 1e-06, 1e-10)
G = 9.81
_DENSE_NODES = np.array([0.0, 1 / 3, 2 / 3, 1.0])
_CUBIC_INVERSE = np.linalg.inv(np.vander(_DENSE_NODES, 4, increasing=True))
FAMILIES = ['nominal_triangular', 'uniform_support', 'lower_mode_triangular']
CONDITIONS = {'Q2_strong_design': 0.02, 'Q2_near_boundary_comparison': 0.0002265378784906}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def fingerprint(path):
    return sha(path)

def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def stamp():
    return datetime.now(timezone.utc).isoformat()

@dataclass(frozen=True)
class Case:
    m0: float
    H0: float
    mcrit: float
    rho_s: float
    rho_f: float
    mu: float
    a: float
    sigma0: float
    delta: float
    Kfs_transfer_m_s: float
    water_conversion_mu_Pa_s: float
    theta_deg: float
    chi_m_inv: float
    curve_shift: float = 0.0
    start: float = 27.0
    end: float = 100.0
    mode: str = 'coupled'
    shape_n: float = 2.0
    initial_stress_ratio: float = 0.5

    @property
    def s_init(self):
        return self.initial_stress_ratio * self.E

    @property
    def qcrit(self):
        return (self.shape_n - 1) / (2 * (self.shape_n + 1))

    @property
    def ki(self):
        return self.Kfs_transfer_m_s * self.water_conversion_mu_Pa_s / (self.rho_f * G)

    @property
    def Hs(self):
        return self.m0 * self.H0

    @property
    def E(self):
        return (self.rho_s - self.rho_f) * G * math.cos(math.radians(self.theta_deg)) * self.Hs

class Forcing:

    def __init__(self, path):
        d = pd.read_csv(path)
        self.t = d.time_since_refined_onset_s.to_numpy(float)
        self.u = d.velocity_m_s.to_numpy(float)
        if not (np.all(np.diff(self.t) > 0) and np.all(self.u >= 0)):
            raise ValueError('Forcing requires ordered times and nonnegative velocities')

    def speed(self, case, t):
        if case.start < self.t[0] or case.end > self.t[-1]:
            raise ValueError('Requested forcing window exceeds digitized data; no extrapolation')
        return max(0.0, float(np.interp(t + case.start, self.t, self.u)) + case.curve_shift)

def constitutive(case, m, s, gamma):
    if not 0 < m < 1 or s <= -case.sigma0:
        raise ValueError('Constitutive state outside domain')
    alpha = case.a / (m * (s + case.sigma0))
    N = 0.0 if gamma == 0 else case.mu * gamma / (case.rho_s * gamma ** 2 * case.delta ** 2 + max(s, 0.0))
    meq = case.mcrit / (1 + math.sqrt(N))
    k = case.ki * math.exp((case.m0 - m) / 0.04)
    return (alpha, N, meq, k)

def closure(case, m, s):
    H = case.Hs / m
    zeta = s / case.E
    if zeta >= case.qcrit:
        e = (case.shape_n + 1) / case.shape_n * (case.E / 2 - s)
        r = 1.0
        Jfactor = case.shape_n * e / H
    else:
        r = math.sqrt(max(0.0, zeta / case.qcrit))
        e = case.E * (1 - (case.shape_n - 1) * r / case.shape_n)
        Jfactor = case.E / H
    return (H, e, r, Jfactor)

def rhs(case, forcing, t, y):
    m, s, _ = y
    U = forcing.speed(case, t)
    gamma = 0.0 if case.mode == 'zero_shear' else case.chi_m_inv * U
    alpha, N, meq, k = constitutive(case, m, s, gamma)
    H, e, r, Jfactor = closure(case, m, s)
    if case.mode == 'drained':
        dm = -m * gamma * (m - meq)
        J = H * dm / m
        ds = 0.0
    else:
        J = k * Jfactor / case.mu
        dm = m * J / H
        ds = (gamma * (m - meq) + J / H) / alpha
    return (dm, ds, J)

def profile(case, m, s, x):
    H, e, r, _ = closure(case, m, s)
    if s / case.E >= case.qcrit:
        pe = e * (1 - x ** case.shape_n)
        effective = case.E * (1 - x) - pe
    else:
        rn = max(r, 1e-100)
        pe = np.where(x <= r, e - case.E * x ** case.shape_n / (case.shape_n * rn ** (case.shape_n - 1)), case.E * (1 - x))
        effective = np.where(x <= r, case.E * ((case.shape_n - 1) * r / case.shape_n - x + x ** case.shape_n / (case.shape_n * rn ** (case.shape_n - 1))), 0.0)
    hydro = case.rho_f * G * math.cos(math.radians(case.theta_deg)) * H * (1 - x)
    return (pe, effective, hydro + pe)

def initial_classification(case):
    q0 = case.initial_stress_ratio
    B = case.rho_f * G * math.cos(math.radians(case.theta_deg)) * case.H0
    qmax = 0.5 + B / ((case.shape_n + 1) * case.E)
    if not math.isfinite(q0) or q0 < 0 or q0 > qmax + 1e-12:
        status = 'invalid_initial_state'
    elif math.isclose(q0, 0, rel_tol=0, abs_tol=1e-12):
        status = 'initial_closure_limit'
    elif math.isclose(q0, qmax, rel_tol=0, abs_tol=1e-12):
        status = 'initial_saturation_boundary'
    elif q0 < case.qcrit - 1e-12:
        status = 'initial_positive_thickness_layer'
    elif math.isclose(q0, case.qcrit, rel_tol=0, abs_tol=1e-12):
        status = 'initial_first_layer_boundary'
    else:
        status = 'admissible_unreached_initial_state'
    f0 = 1 - math.sqrt(q0 / case.qcrit) if 0 <= q0 < case.qcrit else 0.0 if q0 >= case.qcrit else None
    if status == 'invalid_initial_state':
        f0 = None
    return {
        'q0': q0,
        's_init_Pa': case.s_init,
        'qcrit': case.qcrit,
        'qmax': qmax,
        'initial_status': status,
        'initial_fluid_supported_fraction': f0,
        'initial_condition_satisfied': None if status in ('initial_closure_limit', 'invalid_initial_state') else status in ('initial_positive_thickness_layer', 'initial_first_layer_boundary'),
        'dynamic_first_hit_elapsed_s': None,
        'identity': 'mathematical closure admission only; not local observed or credible range',
    }

def run_case(case, forcing, numerics, dense=True, keep_solution=False):
    if not (0 < case.m0 < 1 and 0 < case.mcrit < 1 and (case.H0 > 0) and (case.ki > 0) and (case.mu > 0) and (case.shape_n > 1)):
        raise ValueError('Invalid physical input')
    admission = initial_classification(case)
    if admission['initial_status'] != 'admissible_unreached_initial_state':
        raise ValueError('Initial state is classified separately: ' + admission['initial_status'])
    if case.mode == 'drained' and case.initial_stress_ratio != 0.5:
        raise ValueError('Historical ideal-drained control is restricted to hydrostatic initial state')

    def zero_stress(t, y):
        return y[1]
    zero_stress.terminal = True
    zero_stress.direction = -1

    def desaturation(t, y):
        H, e, r, _ = closure(case, y[0], y[1])
        return e + case.rho_f * G * math.cos(math.radians(case.theta_deg)) * H / case.shape_n
    desaturation.terminal = True
    desaturation.direction = -1

    def fluid_layer(t, y):
        return y[1] / case.E - case.qcrit
    fluid_layer.terminal = False
    fluid_layer.direction = -1
    duration = case.end - case.start
    sol = solve_ivp(lambda t, y: rhs(case, forcing, t, y), (0, duration), (case.m0, case.s_init, 0.0), method='Radau', rtol=numerics['rtol'], atol=(1e-11, 1e-06, 1e-10), max_step=numerics['max_step_s'], dense_output=True, events=(zero_stress, desaturation, fluid_layer))
    if not sol.success:
        raise RuntimeError(sol.message)
    state = 'completed_conditional_window'
    if len(sol.t_events[0]):
        state = 'mean_stress_zero_closure_limit'
    if len(sol.t_events[1]):
        state = 'saturation_hypothesis_boundary'
    if dense:
        ts = np.unique(np.r_[np.arange(0, sol.t[-1], numerics['output_step_s']), sol.t[-1], sol.t_events[2]])
    else:
        ts = np.unique(np.r_[sol.t, sol.t_events[2]])
    ys = sol.sol(ts)
    rows = []
    residual = []
    for t, (m, s, budget) in zip(ts, ys.T):
        H, e, r, _ = closure(case, m, s)
        U = forcing.speed(case, t)
        gamma = 0.0 if case.mode == 'zero_shear' else case.chi_m_inv * U
        alpha, N, meq, k = constitutive(case, m, s, gamma)
        dm, ds, J = rhs(case, forcing, t, (m, s, budget))
        if case.mode == 'drained':
            e = 0.0
            r = 1.0
        resid = gamma * (m - meq) - alpha * ds + dm / m
        residual.append(resid)
        rows.append({
            'elapsed_s': t,
            'source_relative_s': t + case.start,
            'U_input_m_s': U,
            'gamma_input_s_inv': gamma,
            'm': m,
            'H_m': H,
            'mean_effective_stress_Pa': s,
            'mean_effective_ratio': s / case.E,
            'basal_excess_Pa': e,
            'basal_excess_ratio': e / case.E,
            'fluid_supported_fraction': 1 - r,
            'drainage_flux_m_s': J,
            'cumulative_outflow_m': budget,
            'skeleton_thickness_m': H * r,
            'skeleton_to_delta_ratio': H * r / case.delta,
            'alpha_Pa_inv': alpha,
            'k_m2': k,
            'N': N,
            'meq': meq,
            'solid_volume_residual_m': m * H - case.Hs,
            'water_budget_residual_m': H - case.H0 + budget,
            'constitutive_residual_s_inv': resid,
        })
    frame = pd.DataFrame(rows)
    x = np.linspace(0, 1, 1001)
    profile_checks = []
    for j in np.unique(np.linspace(0, len(frame) - 1, min(21, len(frame))).astype(int)):
        row = frame.iloc[j]
        pe, eff, total = profile(case, row.m, row.mean_effective_stress_Pa, x)
        if case.mode == 'drained':
            eff = case.E * (1 - x)
            total = case.rho_f * G * math.cos(math.radians(case.theta_deg)) * row.H_m * (1 - x)
        profile_checks.append((np.min(eff), np.min(total), abs(np.trapezoid(eff, x) - row.mean_effective_stress_Pa) / case.E))
    onset = float(sol.t_events[2][0]) if len(sol.t_events[2]) else None
    summary = {
        'q0': case.initial_stress_ratio,
        's_init_Pa': case.s_init,
        'initial_status': admission['initial_status'],
        'status': state,
        'onset_elapsed_s': onset,
        'onset_source_relative_s': None if onset is None else onset + case.start,
        'max_fluid_fraction': float(frame.fluid_supported_fraction.max()),
        'min_mean_effective_ratio': float(frame.mean_effective_ratio.min()),
        'final_m': float(frame.m.iloc[-1]),
        'final_H_m': float(frame.H_m.iloc[-1]),
        'end_elapsed_s': float(sol.t[-1]),
        'ki_m2': case.ki,
        'E_Pa': case.E,
        'max_solid_residual_m': float(frame.solid_volume_residual_m.abs().max()),
        'max_water_residual_m': float(frame.water_budget_residual_m.abs().max()),
        'max_constitutive_residual_s_inv': float(max((abs(x) for x in residual))),
        'min_profile_effective_stress_Pa': float(min((x[0] for x in profile_checks))),
        'min_profile_total_pore_Pa': float(min((x[1] for x in profile_checks))),
        'max_profile_mean_relative_error': float(max((x[2] for x in profile_checks))),
        'min_skeleton_to_delta_ratio': float(frame.skeleton_to_delta_ratio.min()),
        'continuum_scale_warning': bool(frame.skeleton_to_delta_ratio.min() < 10),
        'continuum_warning_identity': '10 delta is a disclosed resolution diagnostic, not a calibrated formation boundary',
    }
    if keep_solution:
        return (frame, summary, sol)
    return (frame, summary)
Q2_NAMESPACE = SimpleNamespace(Case=Case, Forcing=Forcing, constitutive=constitutive, closure=closure, rhs=rhs, initial_classification=initial_classification, run_case=run_case, G=G)
Q2_CFG = None
FINAL = None

def configure(data_dir):
    global DATA, Q2_CFG, FINAL
    DATA = Path(data_dir)
    for name in ('q2_config.json', 'forcing.csv', 'frozen_validation.json'):
        if not (DATA / name).is_file():
            raise FileNotFoundError('required delivered input missing: data/' + name)
    Q2_CFG = json.loads((DATA / 'q2_config.json').read_text())
    FINAL = json.loads((DATA / 'frozen_validation.json').read_text())
    if Q2_CFG['version'] != 'Q2-ACTIVE-DATA-20261029-v3':
        raise ValueError('wrong frozen Q2 version')
    expected = json.loads((DATA / 'INPUT_HASHES.json').read_text())
    for name, value in expected.items():
        if sha(DATA / name) != value:
            raise ValueError('input fingerprint mismatch: ' + name)
    forcing = pd.read_csv(DATA / 'forcing.csv')
    required = ['time_since_refined_onset_s', 'velocity_m_s']
    if any((c not in forcing for c in required)):
        raise ValueError('forcing requires ' + str(required))
    a = forcing[required].to_numpy(float)
    if not np.isfinite(a).all() or not np.all(np.diff(a[:, 0]) > 0) or np.any(a[:, 1] < 0):
        raise ValueError('invalid forcing values/order/units')
    if a[0, 0] > 27 or a[-1, 0] < 100:
        raise ValueError('forcing does not cover frozen source-relative27-100s')

def load_frozen_q2():
    return Q2_NAMESPACE

def frozen_config():
    if Q2_CFG is None:
        configure(DATA)
    return Q2_CFG

def load_default_forcing():
    f = Forcing(DATA / 'forcing.csv')
    f.identity = 'same_event_retrospective_digitized_Q2_v3'
    f.metadata = {
        'source_file': 'data/forcing.csv',
        'source_sha256': sha(DATA / 'forcing.csv'),
        'not_available_for_pre_event_forecast': True,
    }
    return f

def audit_dense_first_hit(sol, *, stress_scale_E, qcrit, rtol, atol_stress=ATOL[1], state_index=1, known_hit_time=None):
    answer = {
        'method': 'accepted_step_Radau_cubic_stationary_point_audit',
        'audit_status': 'separated_no_contact',
        'recovered_crossing_tau_s': None,
        'contact_candidates': [],
        'min_interpolant_margin_q': None,
        'max_reconstruction_residual_q': 0.0,
        'max_contact_tolerance_q': 0.0,
        'accepted_steps_checked': 0,
        'tolerance_identity': '10*(atol_s/E + rtol*max|q|) heuristic; no certified global error',
        'continuous_physical_first_hit_certified': False,
    }
    if sol.sol is None or len(sol.t) < 2 or stress_scale_E <= 0:
        answer.update(audit_status='audit_unavailable_unknown')
        return answer
    interpolants = getattr(sol.sol, 'interpolants', None)
    if interpolants is None or any((type(p).__name__ != 'RadauDenseOutput' for p in interpolants)):
        answer.update(audit_status='unsupported_dense_interpolant_unknown')
        return answer
    cutoff = float(sol.t[-1]) if known_hit_time is None else min(float(sol.t[-1]), float(known_hit_time))
    crossings = []
    minima = []
    bad_reconstruction = False
    for a, b0 in zip(sol.t[:-1], sol.t[1:]):
        a = float(a)
        b = min(float(b0), cutoff)
        if a >= cutoff or b <= a:
            continue
        width = b - a
        times = a + width * _DENSE_NODES
        q = np.asarray(sol.sol(times)[state_index], dtype=float) / stress_scale_E
        if not np.isfinite(q).all():
            answer['audit_status'] = 'nonfinite_dense_path_unknown'
            return answer
        values = q - qcrit
        coeff = _CUBIC_INVERSE @ values
        mid_margin = float(sol.sol(a + 0.5 * width)[state_index] / stress_scale_E - qcrit)
        reconstruction = abs(float(np.polynomial.polynomial.polyval(0.5, coeff)) - mid_margin)
        tol = max(128 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(q)))), 10 * (float(atol_stress) / stress_scale_E + float(rtol) * float(np.max(np.abs(q)))))
        answer['max_reconstruction_residual_q'] = max(answer['max_reconstruction_residual_q'], reconstruction)
        answer['max_contact_tolerance_q'] = max(answer['max_contact_tolerance_q'], tol)
        answer['accepted_steps_checked'] += 1
        if reconstruction > tol:
            bad_reconstruction = True
        extrema = []
        for root in np.polynomial.polynomial.polyroots([coeff[1], 2 * coeff[2], 3 * coeff[3]]):
            if abs(root.imag) < 1e-09 and 0 < float(root.real) < 1:
                extrema.append(float(root.real))
        cuts = sorted(set([0.0, *extrema, 1.0]))
        cut_times = [a + width * s for s in cuts]
        margins = [float(sol.sol(t)[state_index] / stress_scale_E - qcrit) for t in cut_times]
        local_min = min(margins)
        minima.append(local_min)
        for j, (s, t, g) in enumerate(zip(cuts, cut_times, margins)):
            if abs(g) <= tol:
                is_interior = 0 < s < 1
                if known_hit_time is None or is_interior:
                    answer['contact_candidates'].append({
                        'elapsed_s': t,
                        'margin_q': g,
                        'tolerance_q': tol,
                        'location': 'interior_stationary' if is_interior else 'step_boundary',
                        'step_start_s': a,
                        'step_end_s': b,
                    })
        for left, right, gl, gr in zip(cut_times[:-1], cut_times[1:], margins[:-1], margins[1:]):
            if gl >= 0 and gr < -tol:
                root = left if gl == 0 else brentq(lambda t: float(sol.sol(t)[state_index] / stress_scale_E - qcrit), left, right, xtol=1e-12, rtol=1e-14)
                crossings.append(float(root))
    answer['min_interpolant_margin_q'] = min(minima) if minima else None
    if bad_reconstruction:
        answer['audit_status'] = 'cubic_reconstruction_unresolved_unknown'
    elif crossings:
        answer.update(audit_status='recovered_crossing', recovered_crossing_tau_s=min(crossings))
    elif answer['contact_candidates']:
        answer['audit_status'] = 'contact_candidate'
    return answer

def build_case(q0=None, chi_m_inv=None, start_s=None, horizon_s=None):
    q2 = load_frozen_q2()
    cfg = frozen_config()
    values = dict(cfg['reference'])
    values['initial_stress_ratio'] = float(values['initial_stress_ratio'] if q0 is None else q0)
    values['chi_m_inv'] = float(values['chi_m_inv'] if chi_m_inv is None else chi_m_inv)
    if not math.isfinite(values['initial_stress_ratio']):
        raise ValueError('q0 must be finite; NaN/Inf cannot form a serializable conditional member')
    if not math.isfinite(values['chi_m_inv']) or values['chi_m_inv'] < 0:
        raise ValueError('chi must be a finite nonnegative conditional input')
    start = float(cfg['forcing']['start_relative_s'] if start_s is None else start_s)
    horizon = float(cfg['forcing']['end_relative_s'] - start if horizon_s is None else horizon_s)
    if not math.isfinite(start) or not math.isfinite(horizon) or horizon <= 0:
        raise ValueError('A finite positive supported horizon is required')
    values.update(start=start, end=start + horizon, curve_shift=0.0, mode='coupled', shape_n=cfg['initial_state_enhancement']['shape_n'])
    return q2.Case(**values)

class VelocityForcing:

    def __init__(self, times_s, velocity_m_s, *, identity='declared_Q3_input_path', metadata=None):
        self.t = np.asarray(times_s, dtype=float).copy()
        self.u = np.asarray(velocity_m_s, dtype=float).copy()
        if self.t.ndim != 1 or self.t.shape != self.u.shape or self.t.size < 2:
            raise ValueError('Forcing requires two or more paired 1-D time/velocity values')
        if not (np.all(np.isfinite(self.t)) and np.all(np.isfinite(self.u)) and np.all(np.diff(self.t) > 0) and np.all(self.u >= 0)):
            raise ValueError('Forcing requires ordered finite times and nonnegative finite speed')
        self.identity = str(identity)
        self.metadata = dict(metadata or {})
        self.path_sha256 = hashlib.sha256(self.t.astype('<f8').tobytes() + self.u.astype('<f8').tobytes()).hexdigest()

    def speed(self, case, t):
        if case.start < self.t[0] or case.end > self.t[-1]:
            raise ValueError('Requested path window exceeds input support; no extrapolation')
        at = float(t) + case.start
        if at < self.t[0] - 1e-10 or at > self.t[-1] + 1e-10:
            raise ValueError('Solver evaluation exceeds declared forcing support')
        return max(0.0, float(np.interp(at, self.t, self.u)) + case.curve_shift)

def _cut_diagnostics(q2, case, forcing, at, y):
    m, s, water_out = (float(v) for v in y)
    H, _, _, _ = q2.closure(case, m, s)
    dm, ds, flux = q2.rhs(case, forcing, float(at), (m, s, water_out))
    return {
        'elapsed_s': float(at),
        'source_relative_s': float(at + case.start),
        'm': m,
        's_Pa': s,
        'q': s / case.E,
        'H_m': H,
        'qdot_s_inv': float(ds / case.E),
        'drainage_flux_m_s': float(flux),
        'cumulative_outflow_m': water_out,
        'solid_budget_residual_m': float(m * H - case.Hs),
        'water_budget_residual_m': float(H - case.H0 + water_out),
    }

def _terminal_first_hit(q2, case, forcing, numerics):

    def zero(t, y):
        return y[1]
    zero.terminal = True
    zero.direction = -1

    def desaturation(t, y):
        H, e, _, _ = q2.closure(case, y[0], y[1])
        return e + case.rho_f * q2.G * math.cos(math.radians(case.theta_deg)) * H / case.shape_n
    desaturation.terminal = True
    desaturation.direction = -1

    def first_layer(t, y):
        return y[1] / case.E - case.qcrit
    first_layer.terminal = True
    first_layer.direction = -1
    return solve_ivp(lambda t, y: q2.rhs(case, forcing, t, y), (0.0, case.end - case.start), (case.m0, case.s_init, 0.0), method='Radau', rtol=numerics['rtol'], atol=ATOL, max_step=numerics['max_step_s'], dense_output=True, events=(zero, desaturation, first_layer))

def run_member(q0=None, chi_m_inv=None, forcing=None, *, start_s=None, horizon_s=None, mode='terminal_first_hit', numerics=None, trajectory_points=5, case=None):
    started = time.perf_counter()
    q2 = load_frozen_q2()
    cfg = frozen_config()
    if case is None:
        case = build_case(q0, chi_m_inv, start_s, horizon_s)
    else:
        if any((x is not None for x in (q0, chi_m_inv, start_s, horizon_s))):
            raise ValueError('Pass a case or conditional input arguments, not both')
        expected = build_case(case.initial_stress_ratio, case.chi_m_inv, case.start, case.end - case.start)
        if asdict(case) != asdict(expected):
            raise ValueError('Case contains an unfrozen material/closure/mode change')
    if mode not in ('terminal_first_hit', 'frozen_full_window'):
        raise ValueError('Unknown numerical calling mode')
    setting = dict(cfg['numerics'])
    if numerics is not None:
        if set(numerics) - {'rtol', 'max_step_s', 'output_step_s'}:
            raise ValueError('Only documented numerical settings may change')
        setting.update(numerics)
    if any((not math.isfinite(float(setting[k])) or float(setting[k]) <= 0 for k in ('rtol', 'max_step_s', 'output_step_s'))):
        raise ValueError('Numerical settings must be finite and positive')
    forcing = load_default_forcing() if forcing is None else forcing
    admission = q2.initial_classification(case)
    result = {
        'status': 'unknown',
        'tau_s': None,
        'qdot_cross_s_inv': None,
        'first_hit_valid': False,
        'whole_window_valid': None,
        'initial_classification': admission,
        'case': asdict(case),
        'material_and_closure_identity': 'exact_pinned_Q2_reference_no_fit',
        'chi_identity': 'unidentified_conditional_input_transfer_design_not_local_estimate',
        'forcing_identity': getattr(forcing, 'identity', 'caller_supplied_speed_protocol'),
        'forcing_metadata': getattr(forcing, 'metadata', {}),
        'mode': mode,
        'numerics': {**setting, 'method': 'Radau', 'atol': list(ATOL)},
        'solver_status': None,
        'cut_diagnostics': None,
        'hit_diagnostics': None,
        'samples': [],
        'failure': None,
        'dense_event_audit': None,
        'recorded_solver_hit_tau_s': None,
        'horizon_s': case.end - case.start,
        'probability_identity': 'one conditional member; no probability assigned here',
    }
    if admission['initial_status'] in ('initial_first_layer_boundary', 'initial_positive_thickness_layer'):
        result.update(status='initial', solver_status=admission['initial_status'], initial_boundary_time_s=0.0, initial_is_not_future_prediction=True)
    elif admission['initial_status'] != 'admissible_unreached_initial_state':
        result.update(solver_status=admission['initial_status'], failure='Initial closure/saturation limit or illegal state: not a no-hit outcome')
    else:
        try:
            if mode == 'frozen_full_window':
                _, summary, sol = q2.run_case(case, forcing, setting, dense=False, keep_solution=True)
                solver_status = summary['status']
                result['frozen_summary'] = summary
            else:
                sol = _terminal_first_hit(q2, case, forcing, setting)
                solver_status = 'completed_conditional_window'
                if len(sol.t_events[0]):
                    solver_status = 'mean_stress_zero_closure_limit'
                elif len(sol.t_events[1]):
                    solver_status = 'saturation_hypothesis_boundary'
                elif len(sol.t_events[2]):
                    solver_status = 'first_layer_boundary_terminated'
                if not sol.success:
                    solver_status = 'numerical_failure'
            result['solver_status'] = solver_status
            result['solver_evaluations'] = {'nfev': int(sol.nfev), 'njev': int(sol.njev), 'nlu': int(sol.nlu)}
            end = float(sol.t[-1])
            hits = sol.t_events[2]
            exits = [float(t) for events in sol.t_events[:2] for t in events]
            earliest_exit = min(exits) if exits else math.inf
            hit = float(hits[0]) if len(hits) else None
            result['recorded_solver_hit_tau_s'] = hit
            whole_valid = bool(sol.success and (not exits) and math.isclose(end, case.end - case.start, rel_tol=0, abs_tol=1e-09))
            result['whole_window_valid'] = whole_valid
            audit = audit_dense_first_hit(sol, stress_scale_E=case.E, qcrit=case.qcrit, rtol=setting['rtol'], known_hit_time=hit)
            result['dense_event_audit'] = audit
            recovered = audit['recovered_crossing_tau_s']
            if recovered is not None and (hit is None or recovered < hit):
                hit = recovered
                result['event_detection_identity'] = 'Q3_Radau_dense_interpolant_recovered_crossing'
            elif hit is not None:
                result['event_detection_identity'] = 'original_Q2_directional_event_with_dense_prefix_audit'
            unresolved = 'unknown' in audit['audit_status']
            if hit is None and audit['audit_status'] == 'contact_candidate':
                unresolved = True
            if hit is not None:
                early_contacts = [c for c in audit['contact_candidates'] if c['location'] == 'interior_stationary' and c['elapsed_s'] < hit - 1e-09]
                if early_contacts:
                    unresolved = True
                    result['earlier_contact_candidate'] = True
                if recovered is None:
                    at_hit = _cut_diagnostics(q2, case, forcing, hit, sol.sol(hit))
                    last_width = max(1e-12, float(sol.t[-1] - sol.t[-2]))
                    slope_scale = audit['max_contact_tolerance_q'] / last_width
                    if abs(at_hit['qdot_s_inv']) <= slope_scale:
                        unresolved = True
                        result['flat_recorded_event_contact_candidate'] = True
            if unresolved:
                result.update(status='unknown', first_hit_valid=False, failure='Dense numerical path has unresolved near-threshold contact or audit gap; not no_hit')
            elif sol.success and hit is not None and (hit <= earliest_exit + 1e-10):
                result.update(status='hit', tau_s=hit, first_hit_valid=True)
                hit_diag = _cut_diagnostics(q2, case, forcing, hit, sol.sol(hit))
                result['hit_diagnostics'] = hit_diag
                result['qdot_cross_s_inv'] = hit_diag['qdot_s_inv']
                if mode == 'terminal_first_hit' and end < case.end - case.start - 1e-09:
                    result['whole_window_valid'] = None
                    result['post_hit_window'] = 'not_integrated_for_stopped_first_hit_object'
                elif not whole_valid:
                    result['post_hit_window'] = 'later_exit_does_not_erase_earlier_valid_first_hit'
            elif whole_valid:
                result.update(status='no_hit', first_hit_valid=False)
            else:
                result.update(status='unknown', failure='Exited/failed before a valid recorded first hit')
            result['cut_diagnostics'] = _cut_diagnostics(q2, case, forcing, end, sol.sol(end))
            if trajectory_points:
                count = max(2, int(trajectory_points))
                ts = np.unique(np.r_[np.linspace(0, end, count), [] if hit is None else [hit]])
                result['samples'] = [_cut_diagnostics(q2, case, forcing, t, sol.sol(t)) for t in ts]
        except (ValueError, RuntimeError, FloatingPointError, OverflowError) as error:
            result.update(status='unknown', solver_status='input_numerical_or_physical_exception', failure=f'{type(error).__name__}: {error}', whole_window_valid=False)
    result['elapsed_wall_s'] = time.perf_counter() - started
    return result

def convergence_check(q0=0.5, chi_m_inv=0.02, forcing=None, *, start_s=None, horizon_s=None):
    cfg = frozen_config()
    base = dict(cfg['numerics'])
    tight = {**base, 'rtol': base['rtol'] / 10, 'max_step_s': base['max_step_s'] / 2}
    a = run_member(q0, chi_m_inv, forcing, start_s=start_s, horizon_s=horizon_s, numerics=base, trajectory_points=0)
    b = run_member(q0, chi_m_inv, forcing, start_s=start_s, horizon_s=horizon_s, numerics=tight, trajectory_points=0)
    delta = abs(a['tau_s'] - b['tau_s']) if a['tau_s'] is not None and b['tau_s'] is not None else None
    return {
        'base': a,
        'tight': b,
        'status_agrees': a['status'] == b['status'],
        'absolute_tau_difference_s': delta,
        'identity': 'same_frozen_rhs_and_atol; rtol_tightened_10x_and_max_step_halved',
    }

def cdf(x, family):
    u = np.clip((np.asarray(x) + 1.5) / 3, 0, 1)
    if family == 'uniform_support':
        return u
    if family == 'lower_mode_triangular':
        return 2 * u - u * u
    if family == 'nominal_triangular':
        return np.where(u <= 0.5, 2 * u * u, 1 - 2 * (1 - u) ** 2)
    raise ValueError(family)

class BaselineModel:

    def __init__(self, rows):
        self.rows = sorted(rows, key=lambda r: r['d'])
        self.x = np.array([r['d'] for r in self.rows])
        self.status = np.array([r['status'] for r in self.rows])
        self.tau = np.array([np.nan if r['tau_s'] is None else r['tau_s'] for r in self.rows])
        hits = self.status == 'hit'
        self.hx = self.x[hits]
        self.ht = self.tau[hits]
        if np.any(~np.isin(self.status, ['hit', 'no_hit'])):
            raise ValueError('unknown/initial blocks baseline')
        if not len(self.hx):
            raise ValueError('no hit interpolation support')
        self.maximum_sample_reversal_s = float(max(0, np.max(np.diff(self.ht))))
        if np.any(self.status[np.flatnonzero(hits)[0]:] != 'hit'):
            raise ValueError('nonmonotone classes')
        self.interp = PchipInterpolator(self.hx, self.ht, extrapolate=False)
        self.boundary = None if np.all(hits) else [float(self.x[self.x < self.hx[0]][-1]), float(self.hx[0])]

    @classmethod
    def from_json(cls, path):
        return cls(json.loads(Path(path).read_text())['rows'])

    def predict(self, d):
        d = np.asarray(d, dtype=float)
        s = np.full(d.shape, 'unknown', dtype='<U8')
        t = np.full(d.shape, np.nan)
        valid = np.isfinite(d) & (d >= max(-1.5, self.x[0])) & (d <= min(1.5, self.x[-1]))
        if self.boundary is None:
            hit = valid
        else:
            s[valid & (d <= self.boundary[0])] = 'no_hit'
            hit = valid & (d >= self.boundary[1])
        s[hit] = 'hit'
        t[hit] = self.interp(d[hit])
        bad = (s == 'hit') & ~np.isfinite(t)
        s[bad] = 'unknown'
        t[bad] = np.nan
        return (s, t)

    def propagated(self, edges):
        if self.x[0] > -1.5 or self.x[-1] < 1.5:
            raise ValueError('probability propagation requires complete offset support [-1.5,1.5]; partial supplied map is not a probability model')
        edges = np.asarray(edges, dtype=float)
        if edges.ndim != 1 or len(edges) < 2 or (not np.isfinite(edges).all()) or (edges[0] != 0) or (not np.all(np.diff(edges) > 0)) or (edges[-1] < np.max(self.ht)):
            raise ValueError('probability time grid must start at0, increase, and include every finite mapped first-hit time')
        out = []
        for f in FAMILIES:
            values = []
            lower = []
            upper = []
            for t in edges:
                total = 0.0
                lo_mass = 0.0
                hi_mass = 0.0
                for left, right, yl, yr in zip(self.hx[:-1], self.hx[1:], self.ht[:-1], self.ht[1:]):
                    mass = float(cdf(right, f) - cdf(left, f))
                    if max(yl, yr) <= t:
                        total += mass
                        lo_mass += mass
                        hi_mass += mass
                    elif min(yl, yr) <= t:
                        root = brentq(lambda z: float(self.interp(z)) - t, left, right, xtol=1e-12)
                        total += float(cdf(right, f) - cdf(root, f)) if yr < yl else float(cdf(root, f) - cdf(left, f))
                        hi_mass += mass
                values.append(total)
                lower.append(lo_mass)
                upper.append(hi_mass)
            values = np.array(values)
            mass = np.diff(values)
            none = 0 if self.boundary is None else float(cdf(self.boundary[0], f))
            unresolved = 0 if self.boundary is None else float(cdf(self.boundary[1], f) - cdf(self.boundary[0], f))
            out.append({
                'family': f,
                'finite_mass': mass.tolist(),
                'no_hit_mass': none,
                'boundary_unresolved_mass': unresolved,
                'cdf': values.tolist(),
                'sample_segment_cdf_lower': lower,
                'sample_segment_cdf_upper': upper,
            })
        for row in out:
            total = sum(row['finite_mass']) + row['no_hit_mass'] + row['boundary_unresolved_mass']
            if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-08):
                raise ValueError('probability mass does not cover complete declared input support')
        return out

def shortest_set(edges, masses, no_hit, target=0.95):
    edges = np.asarray(edges)
    a = np.atleast_2d(masses)
    n = np.asarray(no_hit)
    K = a.shape[1]
    p = np.c_[np.zeros(len(a)), np.cumsum(a, axis=1)]
    best = None
    operations = 0

    def offer(i, j, inc):
        nonlocal best
        value = (np.zeros(len(a)) if i is None else p[:, j] - p[:, i]) + (n if inc else 0)
        if value.min() + 1e-12 < target:
            return
        spec = {
            'finite_interval_s': None if i is None else [float(edges[i]), float(edges[j])],
            'finite_width_s': 0.0 if i is None else float(edges[j] - edges[i]),
            'includes_no_hit_by_horizon': inc,
            'family_masses': value.tolist(),
            'worst_mass': float(value.min()),
        }
        strict_width = Decimal(0) if i is None else Decimal(str(edges[j])) - Decimal(str(edges[i]))
        key = (strict_width, -round(float(value.min()), 12), int(inc), -1 if i is None else i)
        if best is None or key < best[0]:
            best = (key, spec)
    for inc in (False, True):
        if inc:
            offer(None, None, True)
        j = 0
        for i in range(K):
            j = max(j, i + 1)
            while j <= K and np.min(p[:, j] - p[:, i] + (n if inc else 0)) < target - 1e-12:
                j += 1
                operations += len(a)
            if j > K:
                break
            offer(i, j, inc)
            operations += len(a)
    return {
        'set': None if best is None else best[1],
        'operations_bound_counter': operations,
        'complexity': 'O(RK) nonnegative masses; two category modes',
    }

def central_set(edges, row, alpha=0.05):
    cumulative = np.r_[0, np.cumsum(row['finite_mass'])]
    lo = int(np.searchsorted(cumulative, alpha / 2, side='left'))
    hi = int(np.searchsorted(cumulative, 1 - alpha / 2, side='left'))
    K = len(edges) - 1
    i = max(0, min(K - 1, lo - 1))
    inc = hi > K
    j = K if inc else min(K, hi)
    mass = float(sum(row['finite_mass'][i:j]) + (row['no_hit_mass'] if inc else 0))
    return {
        'finite_interval_s': [float(edges[i]), float(edges[j])],
        'finite_width_s': float(edges[j] - edges[i]),
        'includes_no_hit_by_horizon': bool(inc),
        'nominal_mass': mass,
    }

def build(regime, chi):
    original = load_default_forcing()
    progress = OUT / f'progress_{regime}.json'
    records = {r['d'].hex(): r for r in json.loads(progress.read_text())} if progress.exists() else {}
    start = time.perf_counter()

    def query(d, role):
        d = float(d)
        key = d.hex()
        if key not in records:
            p = VelocityForcing(original.t, np.maximum(0, original.u + d), identity='Q3_v2_coherent_offset', metadata={'offset_m_s': d})
            tick = time.perf_counter()
            r = run_member(0.5, chi, p, start_s=27, horizon_s=73, trajectory_points=0)
            records[key] = {
                'd': d,
                'status': r['status'],
                'tau_s': r['tau_s'],
                'role': role,
                'wall_s': time.perf_counter() - tick,
                'path_sha256': p.path_sha256,
                'solver_status': r['solver_status'],
                'dense_event_audit': r['dense_event_audit'],
            }
            save(OUT / f'progress_{regime}.json', list(records.values()))
        return records[key]
    for x in np.linspace(-1.5, 1.5, 33):
        query(x, 'initial_monotonicity_grid')
    for x in np.linspace(-1.5, 1.5, 65)[1::2]:
        query(x, 'midpoint_monotonicity_check')
    rows = sorted(records.values(), key=lambda r: r['d'])
    hit = [r for r in rows if r['status'] == 'hit']
    non = [r for r in rows if r['status'] == 'no_hit']
    if non and hit:
        lo = max((r['d'] for r in non))
        hi = min((r['d'] for r in hit))
        while hi - lo > 1e-07:
            mid = (lo + hi) / 2
            r = query(mid, 'class_boundary_bisection')
            if r['status'] == 'hit':
                hi = mid
            elif r['status'] == 'no_hit':
                lo = mid
            else:
                break
    for _ in range(40):
        rows = sorted(records.values(), key=lambda r: r['d'])
        h = [r for r in rows if r['status'] == 'hit']
        x = np.array([r['d'] for r in h])
        y = np.array([r['tau_s'] for r in h])
        p = PchipInterpolator(x, y)
        mids = (x[:-1] + x[1:]) / 2
        score = np.abs(p(mids) - (y[:-1] + y[1:]) / 2)
        score[np.diff(x) < 0.0001] = -1
        j = int(np.argmax(score))
        if score[j] < 1e-07:
            break
        query(mids[j], 'adaptive_curvature_split')
    rows = sorted(records.values(), key=lambda r: r['d'])
    model = BaselineModel(rows)
    save(OUT / f'baseline_model_{regime}.json', {
        'regime': regime,
        'chi': chi,
        'q0': 0.5,
        'rows': rows,
        'boundary': model.boundary,
        'maximum_sample_reversal_s': model.maximum_sample_reversal_s,
        'global_monotonicity_proved': False,
        'wall_s': time.perf_counter() - start,
    })
    return model

def edges_from(cfg, horizon):
    early = np.arange(0, cfg['early_end_s'] + cfg['early_step_s'] / 2, cfg['early_step_s'])
    late = np.arange(cfg['late_start_s'], horizon + cfg['late_step_s'] / 2, cfg['late_step_s'])
    edges = np.unique(np.round(np.r_[early, cfg['middle_edges_s'], late, horizon], 9))
    assert edges[0] == 0 and edges[-1] == horizon and np.all(np.diff(edges) > 0)
    return edges
CONFIG = dict(version='Q3-GP-v2.0', seed=20261030, budgets=[8, 16, 32], q0=0.5, conditions={'strong': 0.02, 'near': 0.0002265378784906}, support=[-1.5, 1.5], horizon_s=73.0, pool_count=601, min_query_distance=0.012, test_count=600, test_protocol='all 600 midpoint cells; training on 601 edges; min distance .0025 m/s', boundary_radius=0.08, alpha=0.05, delta=0.05, classifier_length_scale=0.45, regression_length_scale=0.45, classifier_fallback_probability=[0.35, 0.65], regression_fallback_std_s={'strong': 0.0005, 'near': 0.05}, final_test_use='once after freeze; no retuning; numerical scenario test, not iid coverage certificate', probability_identity='declared input families conditional on fixed Q2, not natural disaster probability')

class Oracle:

    def __init__(self, out):
        self.original = load_default_forcing()
        self.cache = {}
        self.rows = []
        self.out = out
        if (out / 'q2_queries.jsonl').exists():
            self.rows = [json.loads(s) for s in (out / 'q2_queries.jsonl').read_text().splitlines()]
            self.cache = {(r['regime'], round(r['d'], 12)): r for r in self.rows}

    def query(self, regime, d, role):
        key = (regime, round(float(d), 12))
        if key in self.cache:
            return self.cache[key]
        f = VelocityForcing(self.original.t, np.maximum(0, self.original.u + d), identity='fixed_Q2_coherent_offset_conditional', metadata={'offset_m_s': float(d)})
        t = time.perf_counter()
        r = run_member(0.5, CONFIG['conditions'][regime], f, start_s=27.0, horizon_s=73.0, trajectory_points=0)
        row = dict(regime=regime, d=float(d), role=role, status=r['status'], tau_s=r['tau_s'], wall_s=time.perf_counter() - t, input_path_sha256=f.path_sha256, solver_status=r['solver_status'], first_hit_valid=r['first_hit_valid'], whole_window_valid=r['whole_window_valid'], qdot_cross_s_inv=r['qdot_cross_s_inv'], failure=r['failure'], dense_event_audit=r['dense_event_audit'])
        self.rows.append(row)
        self.cache[key] = row
        with (self.out / 'q2_queries.jsonl').open('a') as h:
            h.write(json.dumps(row, allow_nan=False) + '\n')
        return row

class TwoPartGP:

    def fit(self, rows):
        self.rows = rows
        self.x = np.array([r['d'] for r in rows])
        self.lo = min(self.x)
        self.hi = max(self.x)
        known = [r for r in rows if r['status'] in ('hit', 'no_hit')]
        if not known:
            raise ValueError('no known future outcomes')
        X = np.array([r['d'] for r in known])[:, None]
        y = np.array([r['status'] == 'hit' for r in known], int)
        self.single = int(y[0]) if len(np.unique(y)) == 1 else None
        self.clf = None
        if self.single is None:
            self.clf = GaussianProcessClassifier(1.0 * Matern(0.45, nu=1.5), optimizer=None, random_state=CONFIG['seed'], max_iter_predict=100).fit(X, y)
        hits = [r for r in known if r['status'] == 'hit']
        self.reg = None
        if len(hits) >= 2:
            self.reg = GaussianProcessRegressor(ConstantKernel(1.0, (0.01, 100)) * Matern(0.45, (0.03, 3.0), nu=2.5), alpha=1e-08, normalize_y=True, n_restarts_optimizer=0, random_state=CONFIG['seed'])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always', ConvergenceWarning)
                self.reg.fit(np.array([r['d'] for r in hits])[:, None], [r['tau_s'] for r in hits])
            self.warnings = [str(w.message) for w in caught]
        else:
            self.warnings = []
        return self

    def predict(self, d):
        d = np.asarray(d, float)
        X = d[:, None]
        p = np.full(len(d), float(self.single)) if self.clf is None else self.clf.predict_proba(X)[:, 1]
        mu, std = (np.full(len(d), np.nan), np.full(len(d), np.inf)) if self.reg is None else self.reg.predict(X, return_std=True)
        status = np.where(p >= 0.5, 'hit', 'no_hit').astype(object)
        bad = (d < self.lo) | (d > self.hi) | ~np.isfinite(d)
        bad |= (status == 'hit') & (~np.isfinite(mu) | (mu <= 0) | (mu > 73))
        status[bad] = 'unknown'
        mu = np.where(status == 'hit', mu, np.nan)
        return (status, mu, p, std)

def train(oracle, regime, method, out):
    construction_started = time.perf_counter()
    pool = np.linspace(-1.5, 1.5, CONFIG['pool_count'])
    indices = np.rint(np.linspace(0, 600, 8)).astype(int).tolist()
    rows = [oracle.query(regime, pool[i], f'{method}_train') for i in indices]
    snapshots = []
    started = time.perf_counter()
    qtime = sum((r['wall_s'] for r in rows))
    fit_time = 0.0
    acq_time = 0.0
    for n in CONFIG['budgets']:
        while len(rows) < n:
            t = time.perf_counter()
            gp = TwoPartGP().fit(rows)
            fit_time += time.perf_counter() - t
            t = time.perf_counter()
            dist = np.min(abs(pool[:, None] - np.array([r['d'] for r in rows])), axis=1)
            if method == 'active':
                _, _, p, s = gp.predict(pool)
                finite = np.isfinite(s)
                scale = np.max(s[finite]) if finite.any() else 1.0
                uncertainty = np.where(finite, s / max(scale, 1e-15), 0)
                score = (4 * np.exp(-((p - 0.5) / 0.15) ** 2) + uncertainty + 0.05) * np.minimum(dist / 0.06, 1)
                score[dist < CONFIG['min_query_distance']] = -np.inf
            else:
                score = dist
            i = int(np.argmax(score))
            acq_time += time.perf_counter() - t
            row = oracle.query(regime, pool[i], f'{method}_train')
            qtime += row['wall_s']
            rows.append(row)
        t = time.perf_counter()
        gp = TwoPartGP().fit(rows)
        fit_time += time.perf_counter() - t
        snapshots.append((n, gp, dict(regime=regime, method=method, budget=n, q2_nominal_calls=n, q2_nominal_wall_s=qtime, fit_wall_s=fit_time, acquisition_wall_s=acq_time, total_nominal_training_wall_s=time.perf_counter() - construction_started + sum((r['wall_s'] for r in rows if r['role'] != f'{method}_train')), actual_elapsed_training_wall_s=time.perf_counter() - construction_started, single_class=gp.single, regression_kernel=str(gp.reg.kernel_) if gp.reg else None, warnings=gp.warnings)))
    save(out / f'{regime}_{method}_training.json', rows)
    return snapshots

def metrics(truth, status, tau, mask):
    mask = np.asarray(mask, bool)
    hit = truth.status.to_numpy() == 'hit'
    known = truth.status.isin(['hit', 'no_hit']).to_numpy()
    both = mask & hit & (status == 'hit') & np.isfinite(tau)
    err = tau[both] - truth.tau_s.to_numpy()[both]
    k = mask & known
    return dict(n=int(mask.sum()), n_true_hit=int((mask & hit).sum()), n_both_hit=int(both.sum()), mae_s=float(np.mean(abs(err))) if len(err) else None, max_abs_error_s=float(max(abs(err))) if len(err) else None, rmse_s=float(np.sqrt(np.mean(err ** 2))) if len(err) else None, false_negative=int((mask & hit & (status != 'hit')).sum()), false_positive=int((mask & (truth.status.to_numpy() == 'no_hit') & (status == 'hit')).sum()), unknown_predictions=int((mask & (status == 'unknown')).sum()), class_error_rate=float(np.mean(status[k] != truth.status.to_numpy()[k])) if k.any() else None)

def distribution_metrics(truth, status, tau):
    d = truth.d.to_numpy()
    z = (d + 1.5) / 3
    weights = {
        'nominal_triangular': np.where(z <= 0.5, 4 * z, 4 * (1 - z)),
        'uniform_support': np.ones(len(d)),
        'lower_mode_triangular': 2 * (1 - z),
    }
    edges = np.r_[np.linspace(0, 1.5, 301), np.linspace(60, 73, 261)]
    rows = []
    for family, w in weights.items():
        w = w / w.sum()
        truehit = truth.status.to_numpy() == 'hit'
        predhit = status == 'hit'
        actual = np.array([w[truehit & (truth.tau_s.to_numpy() <= t)].sum() for t in edges])
        pred = np.array([w[predhit & (tau <= t)].sum() for t in edges])

        def end(s, t):
            sel = s == 'hit'
            order = np.argsort(t[sel])
            times = t[sel][order]
            m = w[sel][order]
            cumulative = np.cumsum(m)
            return [float(times[min(len(times) - 1, np.searchsorted(cumulative, q))]) if len(times) and cumulative[-1] >= q - 1e-14 else None for q in [0.025, 0.5, 0.975]]
        a = end(truth.status.to_numpy(), truth.tau_s.to_numpy())
        b = end(status, tau)
        rows.append(dict(family=family, no_hit_mass_true=float(w[truth.status.to_numpy() == 'no_hit'].sum()), no_hit_mass_pred=float(w[status == 'no_hit'].sum()), unknown_mass_pred=float(w[status == 'unknown'].sum()), no_hit_mass_abs_error=float(abs(w[status == 'no_hit'].sum() - w[truth.status.to_numpy() == 'no_hit'].sum())), cdf_sup_error=float(max(abs(actual - pred))), finite_subcdf_quantiles_true=a, finite_subcdf_quantiles_pred=b, quantile_endpoint_abs_error_s=[abs(x - y) if x is not None and y is not None else None for x, y in zip(a, b)]))
    return rows

def finite(v, name, nonnegative=False):
    if isinstance(v, bool) or not isinstance(v, (float, int)) or (not math.isfinite(v)):
        raise ValueError(f'{name}: finite number required')
    if nonnegative and v < 0:
        raise ValueError(f'{name}: must be >=0')
    return float(v)

def text(v, name):
    if not isinstance(v, str) or not v.strip():
        raise ValueError(f'{name}: nonempty text required')
    return v

def connect(q3, q4):
    for document in (q3, q4):
        if document.get('time_unit') != 's' or document.get('distance_unit') != 'm':
            raise ValueError('explicit SI time/distance units required')
        text(document.get('time_origin_id'), 'time_origin_id')
    if q3['time_origin_id'] != q4['time_origin_id']:
        raise ValueError('time origins differ')
    entries = {}
    for x in q3['entries']:
        sid = text(x.get('scenario_id'), 'scenario_id')
        if sid in entries:
            raise ValueError('duplicate Q3 scenario')
        for k in ('scenario_hash', 'entry_id', 'entry_state_hash', 'bridge_evidence_id'):
            text(x.get(k), k)
        if x.get('stage') != 'propagation_entry':
            raise ValueError('Q3 must provide explicit propagation_entry, not local threshold')
        if x.get('status') not in ('entry_reached', 'not_entered', 'unknown'):
            raise ValueError('invalid entry status')
        if x['status'] == 'entry_reached':
            finite(x.get('entry_time_s'), 'entry_time_s', True)
        elif x.get('entry_time_s') is not None:
            raise ValueError('unknown/not-entered entry must not carry invented time')
        entries[sid] = x
    out, seen = ([], set())
    for y in q4['arrivals']:
        sid = text(y.get('scenario_id'), 'scenario_id')
        target = text(y.get('target_id'), 'target_id')
        key = (sid, target)
        if key in seen:
            raise ValueError('duplicate scenario-target')
        seen.add(key)
        if sid not in entries:
            raise ValueError('unmatched Q4 scenario')
        x = entries[sid]
        for k in ('scenario_hash', 'entry_id', 'entry_state_hash'):
            if x[k] != y.get(k):
                raise ValueError(f'{k} differs')
        if y.get('stage') != 'target_arrival':
            raise ValueError('Q4 must label target_arrival stage')
        for k in ('terrain_version', 'propagation_version', 'distance_definition', 'arrival_criterion'):
            text(y.get(k), k)
        distance = finite(y.get('distance_m'), 'distance_m', True)
        pub = finite(y.get('publication_time_s'), 'publication_time_s', True)
        act = finite(y.get('action_delay_s'), 'action_delay_s', True)
        if y.get('status') not in ('arrived', 'not_arrived', 'out_of_domain', 'failed', 'unknown'):
            raise ValueError('invalid arrival status')
        if y['status'] == 'arrived' and x['status'] == 'not_entered':
            raise ValueError('arrived contradicts not_entered')
        arrival, lead = (None, None)
        if y['status'] == 'arrived':
            travel = finite(y.get('travel_time_s'), 'travel_time_s', True)
            if x['status'] == 'entry_reached':
                arrival = x['entry_time_s'] + travel
                lead = arrival - pub - act
        elif y.get('travel_time_s') is not None:
            raise ValueError('non-arrival must not carry a finite travel time')
        out.append({
            **y,
            'entry_status': x['status'],
            'entry_time_s': x.get('entry_time_s'),
            'bridge_evidence_id': x['bridge_evidence_id'],
            'distance_m': distance,
            'arrival_time_s': arrival,
            'effective_lead_s': lead,
            'positive_effective_lead': None if lead is None else lead > 0,
            'connection_status': 'computed_conditional' if lead is not None else 'pending_or_nonarrival',
        })
    missing = sorted(set(entries) - {r['scenario_id'] for r in out})
    return {
        'time_origin_id': q3['time_origin_id'],
        'time_unit': 's',
        'distance_unit': 'm',
        'rows': out,
        'scenarios_without_q4': missing,
        'scope': 'conditional scenario calculation; no probability or warning radius inferred',
    }

def record_query(regime, d, role='numerical_test'):
    f0 = load_default_forcing()
    f = VelocityForcing(f0.t, np.maximum(0, f0.u + float(d)), identity='delivered_fixed_Q2_offset', metadata={'offset_m_s': float(d)})
    r = run_member(0.5, CONFIG['conditions'][regime], f, start_s=27.0, horizon_s=73.0, trajectory_points=0)
    return {
        'regime': regime,
        'd': float(d),
        'role': role,
        'status': r['status'],
        'tau_s': r['tau_s'],
        'qdot_cross_s_inv': r['qdot_cross_s_inv'],
        'wall_s': r['elapsed_wall_s'],
        'input_path_sha256': f.path_sha256,
        'solver_status': r['solver_status'],
        'failure': r['failure'],
        'whole_window_valid': r['whole_window_valid'],
    }

def grid_edges():
    return edges_from({
        'early_end_s': 1.5,
        'early_step_s': 0.005,
        'middle_edges_s': [2, 5, 10, 20, 30, 40, 50, 60],
        'late_start_s': 60.0,
        'late_step_s': 0.05,
    }, 73.0)

def categorical_masses(status, tau, d):
    edges = grid_edges()
    bins = []
    none = []
    unknown = []
    h = np.asarray(status) == 'hit'
    nn = np.asarray(status) == 'no_hit'
    uu = ~(h | nn)
    xedges = np.linspace(-1.5, 1.5, len(d) + 1)
    for family in FAMILIES:
        w = np.diff(cdf(xedges, family))
        b = np.histogram(np.asarray(tau)[h], bins=edges, weights=w[h])[0]
        bins.append(b)
        none.append(float(w[nn].sum()))
        unknown.append(float(w[uu].sum()))
    return (np.array(bins), np.array(none), unknown)

def common_probability_task(status, tau, d):
    start = time.perf_counter()
    a, n, u = categorical_masses(status, tau, d)
    edges = grid_edges()
    result = {
        'nominal_central': central_set(edges, {'finite_mass': a[0].tolist(), 'no_hit_mass': n[0]}),
        'nominal_shortest': shortest_set(edges, a[:1], n[:1]),
        'robust_shortest': shortest_set(edges, a, n),
        'cdf': np.cumsum(a, axis=1).tolist(),
        'no_hit_mass': n.tolist(),
        'unknown_mass': u,
    }
    return (result, time.perf_counter() - start)

def numerical_pipeline(results):
    global OUT
    OUT = results / 'baseline'
    OUT.mkdir()
    edges = grid_edges()
    baseline = {}
    basemodel = {}
    baselinecost = {}
    for name, chi in CONDITIONS.items():
        reg = {'Q2_strong_design': 'strong', 'Q2_near_boundary_comparison': 'near'}[name]
        t = time.perf_counter()
        model = build(name, chi)
        construction = time.perf_counter() - t
        t = time.perf_counter()
        prop = model.propagated(edges)
        a = np.array([r['finite_mass'] for r in prop])
        n = np.array([r['no_hit_mass'] for r in prop])
        baseline[name] = {
            'propagation': prop,
            'nominal_central': central_set(edges, prop[0]),
            'nominal_shortest': shortest_set(edges, a[:1], n[:1]),
            'robust_shortest': shortest_set(edges, a, n),
            'sample_count': len(model.rows),
            'boundary': model.boundary,
            'maximum_sample_reversal_s': model.maximum_sample_reversal_s,
        }
        integralwall = time.perf_counter() - t
        baselinecost[reg] = {
            'construction_wall_s': construction,
            'continuous_integral_set_wall_s': integralwall,
            'q2_calls': len(model.rows),
        }
        basemodel[reg] = model
        chosen = baseline[name]['robust_shortest']['set']
        candidate = next((x for x in FINAL['candidates'] if x['condition'] == reg))
        if not np.allclose(chosen['finite_interval_s'], candidate['finite_interval_s'], atol=1e-12, rtol=0) or chosen['includes_no_hit_by_horizon'] != candidate['includes_no_hit']:
            raise RuntimeError('frozen candidate differs from recomputed baseline; final validation blocked')
    save(results / 'baseline_results.json', {'edges_s': edges.tolist(), 'conditions': baseline})
    gpout = results / 'gp'
    gpout.mkdir()
    oracle = Oracle(gpout)
    models = {}
    training = []
    nested = {}
    for reg in CONFIG['conditions']:
        nested[reg] = {}
        for method in ['passive', 'active']:
            models[reg, method] = train(oracle, reg, method, gpout)
            nested[reg][method] = models[reg, method][-1][1].rows
            training.extend((x[2] for x in models[reg, method]))
            print('training', reg, method, 'complete', flush=True)
    save(results / 'gp_training.json', nested)
    save(results / 'gp_training_summary.json', training)
    save(results / 'FROZEN_GP_BEFORE_DIAGNOSTIC.json', {
        'utc': stamp(),
        'analysis_sha256': sha(__file__),
        'config': CONFIG,
        'training_sha256': sha(results / 'gp_training.json'),
        'test_identity': '600 fixed cell midpoints; not binomial certification',
    })
    alltruth = []
    allpred = []
    allmetrics = []
    basechecks = {}
    distributions = []
    costs = []
    setcomparisons = {}
    for reg in CONFIG['conditions']:
        d = -1.5 + (np.arange(600) + 0.5) * 3 / 600
        tt = time.perf_counter()
        truth = pd.DataFrame([oracle.query(reg, x, 'fixed_midpoint_diagnostic') for x in d])
        directwall = time.perf_counter() - tt
        truth['regime'] = reg
        alltruth.append(truth.drop(columns=['dense_event_audit'], errors='ignore'))
        trueprob, probwall = common_probability_task(truth.status.to_numpy(), truth.tau_s.to_numpy(), d)
        setcomparisons[reg] = {'direct': trueprob}
        costs.append({
            'regime': reg,
            'method': 'Direct Q2',
            'budget': 600,
            'mode': 'direct',
            'construction_wall_s': 0.0,
            'query_wall_s': directwall,
            'probability_set_wall_s': probwall,
            'first_use_total_wall_s': directwall + probwall,
            'repeat_query_total_wall_s': directwall + probwall,
            'construction_q2_calls': 0,
            'fallback_q2_calls': 600,
        })
        t = time.perf_counter()
        bs, bt = basemodel[reg].predict(d)
        querywall = time.perf_counter() - t
        basechecks[reg] = {'all': metrics(truth, bs, bt, np.ones(600, bool))}
        dist = np.min(abs(d[:, None] - basemodel[reg].x), axis=1)
        basechecks[reg]['separated_002'] = metrics(truth, bs, bt, dist >= 0.002)
        basechecks[reg]['exact_overlap_count'] = int((dist < 1e-12).sum())
        basechecks[reg]['within_002_count'] = int((dist < 0.002).sum())
        bp, bpwall = common_probability_task(bs, bt, d)
        setcomparisons[reg]['baseline'] = bp
        bcost = baselinecost[reg]
        costs.append({
            'regime': reg,
            'method': '1D baseline',
            'budget': bcost['q2_calls'],
            'mode': 'interpolation',
            'construction_wall_s': bcost['construction_wall_s'],
            'query_wall_s': querywall,
            'probability_set_wall_s': bpwall,
            'first_use_total_wall_s': bcost['construction_wall_s'] + querywall + bpwall,
            'repeat_query_total_wall_s': querywall + bpwall,
            'construction_q2_calls': bcost['q2_calls'],
            'fallback_q2_calls': 0,
            'continuous_integral_extra_wall_s': bcost['continuous_integral_set_wall_s'],
            'primary_continuous_first_use_total_wall_s': bcost['construction_wall_s'] + bcost['continuous_integral_set_wall_s'],
            'primary_continuous_repeat_use_wall_s': bcost['continuous_integral_set_wall_s'],
        })
        for method in ['passive', 'active']:
            for n, gp, timing in models[reg, method]:
                t = time.perf_counter()
                status, tau, p, std = gp.predict(d)
                querywall = time.perf_counter() - t
                trainx = np.array([r['d'] for r in gp.rows])
                minimum = float(np.min(abs(d[:, None] - trainx)))
                if minimum < 0.002 - 1e-12:
                    raise RuntimeError('GP near-duplicate split guard failed')
                center = None if gp.single is not None else float(np.linspace(-1.5, 1.5, 601)[np.argmin(abs(gp.predict(np.linspace(-1.5, 1.5, 601))[2] - 0.5))])
                boundary = np.zeros(600, bool) if center is None else abs(d - center) <= 0.08
                for mode in ['raw', 'fallback']:
                    fallbackwall = 0.0
                    calls = 0
                    if mode == 'fallback':
                        fallback_loop_started = time.perf_counter()
                        mask = (status == 'unknown') | (p > 0.35) & (p < 0.65) | (status == 'hit') & (std > CONFIG['regression_fallback_std_s'][reg])
                        for j in np.flatnonzero(mask):
                            r = record_query(reg, d[j], 'GP_fallback')
                            fallbackwall += r['wall_s']
                            calls += 1
                            status[j] = r['status']
                            tau[j] = np.nan if r['tau_s'] is None else r['tau_s']
                            with (results / 'fallback_queries.jsonl').open('a') as h:
                                h.write(json.dumps({**r, 'method': method, 'budget': n}) + '\n')
                        fallback_solver_wall_s = fallbackwall
                        fallbackwall = time.perf_counter() - fallback_loop_started
                    t = time.perf_counter()
                    distribution = distribution_metrics(truth, status, tau)
                    sets, pwall = common_probability_task(status, tau, d)
                    pwalltotal = time.perf_counter() - t
                    commonwall = pwall
                    total = timing['total_nominal_training_wall_s'] + querywall + fallbackwall + commonwall
                    fields = {
                        'regime': reg,
                        'method': method,
                        'budget': n,
                        'mode': mode,
                        'fallback_calls': calls,
                        'fallback_wall_s': fallbackwall,
                        'total_nominal_q2_calls': n + calls,
                        'total_nominal_wall_s': total,
                    }
                    for subset, mask in [('all', np.ones(600, bool)), ('boundary', boundary)]:
                        allmetrics.append({**fields, 'subset': subset, **metrics(truth, status, tau, mask)})
                    distributions.extend(({**fields, **r} for r in distribution))
                    costs.append({
                        'regime': reg,
                        'method': method + ' GP',
                        'budget': n,
                        'mode': mode,
                        'construction_wall_s': timing['total_nominal_training_wall_s'],
                        'query_wall_s': querywall + fallbackwall,
                        'probability_set_wall_s': commonwall,
                        'first_use_total_wall_s': total,
                        'repeat_query_total_wall_s': querywall + fallbackwall + commonwall,
                        'construction_q2_calls': n,
                        'fallback_q2_calls': calls,
                        'metric_diagnostic_extra_wall_s': pwalltotal - commonwall,
                    })
                    if n == 32:
                        setcomparisons[reg][method + '_' + mode] = sets
                    for j, x in enumerate(d):
                        allpred.append({
                            'regime': reg,
                            'method': method,
                            'budget': n,
                            'mode': mode,
                            'd': x,
                            'true_status': truth.status.iloc[j],
                            'true_tau_s': truth.tau_s.iloc[j],
                            'status': status[j],
                            'tau_s': tau[j],
                            'p_hit': p[j],
                            'std_s': std[j],
                            'boundary': bool(boundary[j]),
                        })
        print('numerical diagnostic', reg, 'complete', flush=True)
    pd.concat(alltruth, ignore_index=True).to_csv(results / 'test_truth.csv', index=False)
    pd.DataFrame(allpred).to_csv(results / 'gp_predictions.csv.gz', index=False)
    pd.DataFrame(allmetrics).to_csv(results / 'gp_metrics.csv', index=False)
    pd.DataFrame(costs).to_csv(results / 'cost_comparison.csv', index=False)
    save(results / 'distribution_metrics.json', distributions)
    save(results / 'baseline_test_metrics.json', basechecks)
    save(results / 'set_comparison.json', setcomparisons)
    return {
        'baseline_cost': baselinecost,
        'gp_unique_train_and_diagnostic_calls': len(oracle.rows),
        'baseline_calls': sum((x['q2_calls'] for x in baselinecost.values())),
        'gp_fallback_calls': sum((r['fallback_calls'] for r in allmetrics if r['subset'] == 'all')),
    }

def contains_candidate(row, candidate):
    if row['status'] == 'no_hit':
        return bool(candidate['includes_no_hit'])
    if row['status'] != 'hit' or row['tau_s'] is None:
        return False
    lo, hi = candidate['finite_interval_s']
    return lo <= float(row['tau_s']) <= hi

def validate_final(results, workers):
    frozen = json.loads((DATA / 'frozen_validation.json').read_text())
    children = np.random.SeedSequence(frozen['seed']).spawn(6)
    inputs = []
    k = 0
    for c in frozen['candidates']:
        for family in frozen['families']:
            rng = np.random.default_rng(children[k])
            k += 1
            n = frozen['n_per_candidate_family']
            if family == 'nominal_triangular':
                d = rng.triangular(-1.5, 0.0, 1.5, n)
            elif family == 'uniform_support':
                d = rng.uniform(-1.5, 1.5, n)
            else:
                d = rng.triangular(-1.5, -1.5, 1.5, n)
            inputs.extend(({'candidate': c['name'], 'regime': c['condition'], 'family': family, 'member': j, 'd': float(x)} for j, x in enumerate(d)))
    inputfile = results / 'validation_inputs.csv'
    pd.DataFrame(inputs).to_csv(inputfile, index=False)
    bindings = {
        'utc_before_labels': stamp(),
        'analysis_sha256': sha(__file__),
        'frozen_config_sha256': sha(DATA / 'frozen_validation.json'),
        'validation_inputs_sha256': sha(inputfile),
        'source_binding_sha256': sha(DATA / 'SOURCE_BINDING.json'),
        'forcing_sha256': sha(DATA / 'forcing.csv'),
        'q2_config_sha256': sha(DATA / 'q2_config.json'),
        'label_generator': 'delivered sol.py; equivalence to original Q2 required before certification release',
        'criterion': 'each one-sided CP lower>=.95, all fixed n denominator; any unknown or initial blocks its claim',
        'optional_stopping': False,
        'distribution_parameters': {
            'nominal_triangular': [-1.5, 0.0, 1.5],
            'uniform_support': [-1.5, 1.5],
            'lower_mode_triangular': [-1.5, -1.5, 1.5],
        },
        'input_construction': 'max(0,U_ref_node+d) applied to every provided node, then original piecewise linear interpolation; d constant over complete path',
    }
    save(results / 'FROZEN_VALIDATION_BINDINGS.json', bindings)
    jobs = [(r['regime'], r['d']) for r in inputs]
    start = time.perf_counter()
    rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=configure, initargs=(str(DATA),)) as pool:
        for i, row in enumerate(pool.map(validation_worker, jobs, chunksize=32)):
            r = {**inputs[i], **row}
            rows.append(r)
            with (results / 'validation_members.jsonl').open('a') as h:
                h.write(json.dumps(r, allow_nan=False) + '\n')
            if (i + 1) % 1024 == 0:
                print('random validation', i + 1, '/', len(jobs), flush=True)
    pd.DataFrame(rows).to_csv(results / 'validation_members.csv.gz', index=False)
    claims = []
    for c in frozen['candidates']:
        for family in frozen['families']:
            selected = [r for r in rows if r['candidate'] == c['name'] and r['family'] == family]
            n = len(selected)
            if n != frozen['n_per_candidate_family']:
                raise RuntimeError('incomplete fixed validation')
            accepted = sum((contains_candidate(r, c) for r in selected))
            unknown = sum((r['status'] == 'unknown' for r in selected))
            initial = sum((r['status'] == 'initial' for r in selected))
            lower = 0.0 if accepted == 0 else float(beta.ppf(frozen['per_claim_delta'], accepted, n - accepted + 1))
            claims.append({
                'regime': c['condition'],
                'candidate': c['name'],
                'family': family,
                'n': n,
                'k': accepted,
                'empirical_mass': accepted / n,
                'lower_bound': lower,
                'unknown': unknown,
                'initial': initial,
                'passed': bool(lower >= 1 - frozen['alpha'] and unknown == 0 and (initial == 0)),
            })
    answer = {
        'version': FINAL['version'],
        'status': 'completed_fixed_n',
        'alpha': frozen['alpha'],
        'familywise_delta': frozen['familywise_delta'],
        'per_claim_delta': frozen['per_claim_delta'],
        'claims': claims,
        'all_passed': all((x['passed'] for x in claims)),
        'wall_s': time.perf_counter() - start,
        'bindings_sha256': sha(results / 'FROZEN_VALIDATION_BINDINGS.json'),
        'identity': 'new independent declared-input iid Q2 simulation validation; not natural-event coverage',
    }
    save(results / 'validation.json', answer)
    return answer

def validation_worker(job):
    return record_query(job[0], job[1], 'fixed_n_random_validation')

def lead_outputs(results, validation):
    baseline = json.loads((results / 'baseline_results.json').read_text())['conditions']
    rows = []
    for c in FINAL['candidates']:
        reg = c['condition']
        key = 'Q2_strong_design' if reg == 'strong' else 'Q2_near_boundary_comparison'
        unknown = max((x['boundary_unresolved_mass'] for x in baseline[key]['propagation']))
        claims = [x for x in validation.get('claims', []) if x['regime'] == reg]
        certified = None if not claims else all((x['passed'] for x in claims))
        for pub in FINAL['delays']['publication_s']:
            for act in FINAL['delays']['action_s']:
                lo, hi = c['finite_interval_s']
                a = lo - pub - act
                b = hi - pub - act
                rows.append({
                    'regime': reg,
                    'candidate': c['name'],
                    'publication_delay_s': pub,
                    'action_delay_s': act,
                    'finite_lower_s': a,
                    'finite_upper_s': b,
                    'includes_no_hit': c['includes_no_hit'],
                    'unknown_mass': unknown,
                    'finite_width_s': hi - lo,
                    'finite_all_positive': a > 0,
                    'finite_all_nonpositive': b <= 0,
                    'certification_passed': certified,
                    'identity': 'translated frozen local set; delays hypothetical, no infinite lead assigned to no_hit',
                })
    pd.DataFrame(rows).to_csv(results / 'lead_margins.csv', index=False)

def execution_controls(results):
    checks = {}
    m = BaselineModel([{'d': -1.0, 'status': 'hit', 'tau_s': 2.0}, {'d': 1.0, 'status': 'hit', 'tau_s': 1.0}])
    s, t = m.predict(np.array([-1.5, 1.5]))
    checks['ALG01_incomplete_support_unknown'] = bool(np.all(s == 'unknown') and np.isnan(t).all())
    tiny = shortest_set([0, 1e-11, 2e-11], [[0.5, 0.5]], [0], target=0.5)['set']
    checks['ALG02_strict_tiny_width_primary'] = tiny['finite_width_s'] == 1e-11 and tiny['worst_mass'] == 0.5
    checks['no_benefit_no_hit_category'] = not shortest_set([0, 1, 2], [[0.5, 0.5]], [0], target=0.5)['set']['includes_no_hit_by_horizon']
    candidate = {'finite_interval_s': [0.775, 0.815], 'includes_no_hit': False}
    checks['closed_lower_and_upper'] = all((contains_candidate({'status': 'hit', 'tau_s': t}, candidate) for t in [0.775, 0.815]))
    checks['closed_interval_outside_excluded'] = all((not contains_candidate({'status': 'hit', 'tau_s': t}, candidate) for t in [np.nextafter(0.775, -np.inf), np.nextafter(0.815, np.inf)]))
    checks['nonhit_categories_not_future_hits'] = all((not contains_candidate({'status': s, 'tau_s': None}, candidate) for s in ['no_hit', 'initial', 'unknown']))
    for name, q0, chi, d, expect in [('strong_hit', 0.5, 0.02, 0, 'hit'), ('near_hit', 0.5, 0.0002265378784906, 0.5, 'hit'), ('no_hit', 0.5, 0.0002265378784906, -1, 'no_hit'), ('initial', 0.1, 0.02, 0, 'initial'), ('invalid_initial', -0.1, 0.02, 0, 'unknown')]:
        f = load_default_forcing()
        vf = VelocityForcing(f.t, np.maximum(0, f.u + d))
        r = run_member(q0, chi, vf, trajectory_points=0)
        checks[name] = r['status'] == expect
        checks[name + '_nonhit_no_time'] = r['status'] == 'hit' or r['tau_s'] is None
    try:
        build_case(float('nan'), 0.02)
        checks['nan_initial_rejected'] = False
    except ValueError:
        checks['nan_initial_rejected'] = True

    def tangent(t, y):
        return np.array([2 * (t - 1) / 3])
    sol = solve_ivp(tangent, (0, 2), [1 / 6 + 1 / 3], method='Radau', dense_output=True, first_step=2, max_step=2, rtol=1e-08, atol=1e-12)
    audit = audit_dense_first_hit(sol, stress_scale_E=1.0, qcrit=1 / 6, rtol=1e-08, atol_stress=1e-12, state_index=0)
    checks['tangent_not_no_hit'] = audit['audit_status'] == 'contact_candidate'
    save(results / 'execution_controls.json', {'checks': checks, 'all_passed': all(checks.values()), 'independent_review': 'pending'})
    if not all(checks.values()):
        raise RuntimeError('execution controls failed')

def join_q4_same_scenarios(q3_entries, q4_responses, expected_targets):
    import math
    if not isinstance(expected_targets, dict) or not expected_targets or any((not v for v in expected_targets.values())):
        raise ValueError('expected_targets must map target id to frozen geometry hash')
    targets = list(expected_targets)
    entries = {}
    for row in q3_entries:
        sid = row['scenario_id']
        if sid in entries:
            raise ValueError('duplicate Q3 scenario')
        if row['entry_stage'] != 'propagation_entry':
            raise ValueError('local first hit is not automatically propagation entry')
        if not row.get('bridge_evidence_id'):
            raise ValueError('explicit stage bridge evidence required')
        if row['status'] not in ('entry_reached', 'not_entered', 'unknown'):
            raise ValueError('unknown Q3 entry status')
        entries[sid] = row
    response = {}
    identities = ['scenario_hash', 'joint_input_hash', 'entry_id', 'entry_geometry_hash', 'entry_state_hash', 'time_origin_id']
    for row in q4_responses:
        key = (row['scenario_id'], row['target_id'])
        if key in response:
            raise ValueError('duplicate Q4 scenario/target')
        if key[0] not in entries or key[1] not in targets:
            raise ValueError('unexpected scenario/target')
        q = entries[key[0]]
        if any((q.get(k) != row.get(k) or q.get(k) is None for k in identities)):
            raise ValueError('joint scenario, entry state or clock mismatch')
        if row.get('target_geometry_hash') != expected_targets[key[1]]:
            raise ValueError('target geometry mismatch')
        if row['time_unit'] != 's' or row['distance_unit'] != 'm':
            raise ValueError('unit mismatch')
        if row['target_arrival_status'] not in ('arrived', 'valid_no_arrival_by_horizon', 'unknown', 'not_run'):
            raise ValueError('invalid target arrival status')
        travel = row.get('travel_time_s')
        if row['target_arrival_status'] == 'arrived':
            if travel is None or not math.isfinite(travel) or travel < 0:
                raise ValueError('invalid travel time')
            if q['status'] == 'not_entered':
                raise ValueError('not-entered member cannot have arrived')
            if row.get('prefix_review_status') != 'pass':
                raise ValueError('arrived prefix needs explicit validity review')
            valid_until = row.get('valid_until_s')
            if valid_until is None or not math.isfinite(valid_until) or valid_until < travel:
                raise ValueError('valid prefix does not include arrival')
        elif travel is not None:
            raise ValueError('non-arrived travel time must remain null')
        if row['target_arrival_status'] == 'valid_no_arrival_by_horizon':
            horizon = row.get('horizon_s', row.get('propagation_horizon_s'))
            if horizon is None or not math.isfinite(horizon) or horizon <= 0:
                raise ValueError('valid-noarrival needs a positive finite propagation horizon')
            if 'horizon_s' in row and 'propagation_horizon_s' in row and (row['horizon_s'] != row['propagation_horizon_s']):
                raise ValueError('conflicting propagation horizon aliases')
            valid_until = row.get('valid_until_s')
            if row.get('whole_window_valid') is not True or valid_until is None or (not math.isfinite(valid_until)) or (valid_until < horizon):
                raise ValueError('valid-noarrival requires a valid complete propagation window')
        response[key] = row
    out = []
    for sid, q in entries.items():
        for target in targets:
            r = response.get((sid, target))
            row = {
                'scenario_id': sid,
                'target_id': target,
                'travel_time_s': None,
                'absolute_arrival_s': None,
                'effective_lead_s': None,
                'target_arrival_status': 'unknown',
                'reason': 'missing scenario-target response',
            }
            if r is None:
                out.append(row)
                continue
            row.update(travel_time_s=r.get('travel_time_s'), target_arrival_status=r['target_arrival_status'], run_termination_status=r.get('run_termination_status'), whole_window_valid=bool(r.get('whole_window_valid', False)), prefix_review_status=r.get('prefix_review_status'), horizon_s=r.get('horizon_s', r.get('propagation_horizon_s')), reason='relative Q4 response only')
            t_in = q.get('T_in_s')
            pub = q.get('t_pub_s')
            action = q.get('d_act_s')
            if r['target_arrival_status'] == 'arrived' and q['status'] == 'entry_reached' and (t_in is not None):
                if not math.isfinite(t_in):
                    raise ValueError('nonfinite entry time')
                row['absolute_arrival_s'] = t_in + r['travel_time_s']
                row['reason'] = 'same-scenario arrival; bridge source must still be scientifically reviewed'
                if pub is not None and action is not None:
                    if not math.isfinite(pub) or not math.isfinite(action) or action < 0:
                        raise ValueError('invalid publication/action time')
                    row['effective_lead_s'] = row['absolute_arrival_s'] - pub - action
            out.append(row)
    return out

def q4_five_point_conditional_demo(context, candidates, publication_delays_s, action_delays_s):
    import math
    demo = context['existing_five_point_ODE_demo']
    points = demo['points']
    bridge = float(demo['synthetic_bridge']['T_in_equals_T_L_plus_s'])
    if bridge != 0:
        raise ValueError('this fixed demonstration declares T_in=T_L+0 only')
    if len(points) != 5:
        raise ValueError('expected five fixed pre-existing Q4 ODE checks')
    out = []
    for spec in candidates:
        interval = spec.get('finite_interval_s')
        if interval is not None and (len(interval) != 2 or not all((math.isfinite(x) for x in interval)) or interval[0] > interval[1]):
            raise ValueError('invalid local candidate interval')
        for p in points:
            status = p['target_arrival_status']
            travel = p['travel_time_s']
            if status == 'arrived' and (travel is None or not math.isfinite(travel) or travel < 0 or (travel > p['horizon_s'])):
                raise ValueError('invalid saved Q4 finite arrival')
            if status != 'arrived' and travel is not None:
                raise ValueError('Q4 nonarrival cannot have finite travel time')
            for pub in publication_delays_s:
                for action in action_delays_s:
                    if not math.isfinite(pub) or not math.isfinite(action) or min(pub, action) < 0:
                        raise ValueError('invalid hypothetical delay')
                    arrival = None if interval is None or travel is None else [v + bridge + travel for v in interval]
                    lead = None if arrival is None else [v - pub - action for v in arrival]
                    out.append({
                        'candidate': spec['name'],
                        'condition': spec['condition'],
                        'q4_scenario_id': p['scenario_id'],
                        'target_id': p['target_id'],
                        'target_distance_m': p['target_distance_m'],
                        'distance_definition': p['distance_definition'],
                        'source_local_finite_interval_s': interval,
                        'hypothetical_bridge_T_in_minus_T_L_s': bridge,
                        'publication_delay_s': pub,
                        'action_delay_s': action,
                        'q4_target_arrival_status': status,
                        'q4_horizon_from_T_in_s': p['horizon_s'],
                        'q4_run_termination_status': p['run_termination_status'],
                        'target_finite_arrival_interval_from_shear_s': arrival,
                        'target_finite_lead_interval_s': lead,
                        'includes_local_no_hit': bool(spec['includes_no_hit']),
                        'includes_local_unknown': bool(spec.get('includes_unknown', False)),
                        'local_no_hit_lead_s': None,
                        'unknown_lead_s': None,
                        'conditional_finite_censored_arrival_lower_bound_from_shear_s': interval[0] + bridge + p['horizon_s'] if interval is not None and status == 'valid_no_arrival_by_horizon' else None,
                        'negative_finite_lead_interpretation': 'publication/action cannot finish before target arrival under this hypothesis',
                        'category_interpretation': 'Q3 no-hit and unknown remain separate; Q4 noarrival only by180s from entry, not never',
                        'joint_probability_weight': None,
                        'statistical_identity': 'translation of fixed local candidate within each unweighted synthetic Q4 scenario; not new joint95% or natural-event coverage',
                        'geographic_warning_radius_m': None,
                        'scope': 'explicit conditional scalar-slope demonstration with assumed stage bridge; not actual Gyirong propagation',
                    })
    return out

def connect(q3, q4):
    for doc in (q3, q4):
        if doc.get('time_unit') != 's' or doc.get('distance_unit') != 'm':
            raise ValueError('Q4 explicit SI units required')
    if q3.get('expected_targets') != q4.get('expected_targets'):
        raise ValueError('expected target geometry contracts differ')
    for entry in q3['entries']:
        if entry['status'] == 'entry_reached':
            finite(entry.get('T_in_s'), 'T_in_s')
        elif entry.get('T_in_s') is not None:
            raise ValueError('unknown/not-entered must not carry invented entry time')
        for field in ('t_pub_s', 'd_act_s'):
            if entry.get(field) is not None:
                finite(entry[field], field, field == 'd_act_s')
    for response in q4['arrivals']:
        horizon = finite(response.get('propagation_horizon_s'), 'propagation_horizon_s', True)
        if horizon <= 0:
            raise ValueError('positive propagation horizon required')
        if type(response.get('whole_window_valid')) is not bool:
            raise ValueError('whole_window_valid must be Boolean')
        if response['target_arrival_status'] == 'valid_no_arrival_by_horizon' and (not response['whole_window_valid']):
            raise ValueError('failed/incomplete window cannot certify no arrival')
        if response['target_arrival_status'] == 'arrived':
            travel = finite(response.get('travel_time_s'), 'travel_time_s', True)
            if travel > horizon:
                raise ValueError('arrival exceeds supported horizon')
    rows = join_q4_same_scenarios(q3['entries'], q4['arrivals'], q3['expected_targets'])
    return {'rows': rows, 'scope': 'strict scenario-target connection; stage bridge still needs scientific evidence'}

def q4_outputs(results, context):
    rows = q4_five_point_conditional_demo(context, FINAL['candidates'], FINAL['delays']['publication_s'], FINAL['delays']['action_s'])
    save(results / 'q4_conditional_demo.json', rows)
    table = []
    for r in rows:
        a = r['target_finite_arrival_interval_from_shear_s']
        l = r['target_finite_lead_interval_s']
        table.append({
            'regime': r['condition'],
            'scenario_id': r['q4_scenario_id'],
            'target_id': r['target_id'],
            'position_m': r['target_distance_m'],
            'travel_time_s': None if a is None else a[0] - r['source_local_finite_interval_s'][0],
            'arrival_lower_s': None if a is None else a[0],
            'arrival_upper_s': None if a is None else a[1],
            'publication_delay_s': r['publication_delay_s'],
            'action_delay_s': r['action_delay_s'],
            'effective_lead_lower_s': None if l is None else l[0],
            'effective_lead_upper_s': None if l is None else l[1],
            'status': r['q4_target_arrival_status'],
            'horizon_s': r['q4_horizon_from_T_in_s'],
            'includes_no_hit': r['includes_local_no_hit'],
            'unknown_lead_s': None,
            'terrain_identity': 'synthetic constant8degree slope h=2m, target1000m',
            'stage_bridge': 'hypothetical T_in=T_L+0; not observed Gyirong bridge',
            'tested_target_actionable_for_entire_finite_set': None if l is None else l[0] > 0,
            'natural_warning_radius_m': None,
        })
    pd.DataFrame(table).to_csv(results / 'q4_spatial_demo.csv', index=False)
SENSITIVITY_DESIGN = {
    'version': 'Q3-SENSITIVITY-v4.0',
    'q0_design': [0.17, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.48, 0.49, 0.495, 0.5, 0.505, 0.51, 0.52, 0.55, 0.6, 0.65, 0.7, 0.75],
    'offset_grid_count': 17,
    'offset_support_m_s': [-1.5, 1.5],
    'conditions': {'strong': 0.02, 'near': 0.0002265378784906},
    'chi_multipliers': [0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1],
    'chi_offset_design_m_s': [-1.5, 0.0, 1.5],
    'derivative_steps': {'offset_m_s': [0.0005, 0.001, 0.002], 'q0': [1e-06, 1e-05, 0.0001], 'chi_relative': [0.0001, 0.001, 0.01]},
    'support_halfwidths_m_s': [0.5, 1.0, 1.5],
    'probability_targets': [0.9, 0.95, 0.99],
    'families': ['nominal_triangular', 'uniform_support', 'lower_mode_triangular'],
    'q0_chi_identity': 'unweighted mathematical design conditions, not empirical credible ranges or random parameters',
    'scope': 'fixed Q2 material/equations/threshold/time window; no new natural-event coverage or changed final candidates',
    'state_rule': 'finite tau only for hit; initial/unknown separate; no_hit never filled with horizon',
    'certification_rule': 'new support/target/initial/chi sensitivities do not inherit frozen final CP certificate',
}

def sens_cdf(x, family, halfwidth):
    a = float(halfwidth)
    if not 0 < a <= 1.5:
        raise ValueError('sensitivity support must remain inside original[-1.5,1.5]')
    u = np.clip((np.asarray(x, dtype=float) + a) / (2 * a), 0.0, 1.0)
    if family == 'uniform_support':
        return u
    if family == 'lower_mode_triangular':
        return 2 * u - u * u
    if family == 'nominal_triangular':
        return np.where(u <= 0.5, 2 * u * u, 1 - 2 * (1 - u) ** 2)
    raise ValueError('unknown declared sensitivity family')

def sens_propagate(model, edges, halfwidth):
    from scipy.optimize import brentq
    edges = np.asarray(edges, dtype=float)
    a = float(halfwidth)
    if not 0 < a <= 1.5:
        raise ValueError('sensitivity support must remain inside original[-1.5,1.5]')
    if edges.ndim != 1 or len(edges) < 2 or (not np.isfinite(edges).all()) or (edges[0] != 0) or (not np.all(np.diff(edges) > 0)):
        raise ValueError('sensitivity probability grid requires finite1D increasing edges starting at0')
    if model.x[0] > -a or model.x[-1] < a:
        raise ValueError('incomplete support for sensitivity propagation')
    finite_left = max(float(model.hx[0]), -a)
    finite_right = min(float(model.hx[-1]), a)
    supported_max = None
    if finite_right >= finite_left:
        times = [float(model.interp(finite_left)), float(model.interp(finite_right))]
        times.extend((float(v) for v in model.ht[(model.hx > finite_left) & (model.hx < finite_right)]))
        supported_max = max(times)
        if edges[-1] < supported_max - 1e-12:
            raise ValueError('sensitivity grid omits finite hits inside its current input support')
    out = []
    for family in SENSITIVITY_DESIGN['families']:
        values = []
        for target in edges:
            total = 0.0
            for rawl, rawr in zip(model.hx[:-1], model.hx[1:]):
                left = max(float(rawl), -a)
                right = min(float(rawr), a)
                if right <= left:
                    continue
                yl = float(model.interp(left))
                yr = float(model.interp(right))
                mass = float(sens_cdf(right, family, a) - sens_cdf(left, family, a))
                if max(yl, yr) <= target:
                    total += mass
                elif min(yl, yr) <= target:
                    if yl == target:
                        root = left
                    elif yr == target:
                        root = right
                    else:
                        root = float(brentq(lambda z: float(model.interp(z)) - target, left, right, xtol=1e-12))
                    total += float(sens_cdf(right, family, a) - sens_cdf(root, family, a)) if yr < yl else float(sens_cdf(root, family, a) - sens_cdf(left, family, a))
            values.append(total)
        values = np.asarray(values)
        masses = np.diff(values)
        if np.min(masses) < -1e-10:
            raise ValueError('nonmonotone finite CDF')
        roundoff = float(max(0.0, -masses.min()))
        masses = np.maximum(0.0, masses)
        if model.boundary is None:
            none = unknown = 0.0
        else:
            lo, hi = model.boundary
            none = float(sens_cdf(lo, family, a))
            unknown = float(sens_cdf(hi, family, a) - sens_cdf(lo, family, a))
        total = float(masses.sum() + none + unknown)
        if abs(total - 1) > 1e-08:
            raise ValueError('sensitivity probability mass not closed')
        out.append({
            'family': family,
            'support_halfwidth_m_s': a,
            'finite_mass': masses.tolist(),
            'cdf': values.tolist(),
            'no_hit_mass': none,
            'unknown_mass': unknown,
            'total_mass': total,
            'negative_roundoff_removed': roundoff,
            'finite_support_max_tau_s': supported_max,
            'probability_identity': 'declared narrowed-support distribution propagated through fixed numerical map; no new certificate',
        })
    return out

def sens_certificate_identity(regime, halfwidth, target, method, spec, q0, chi_m_inv, families):
    import math
    expected = {'strong': (0.02, [0.775, 0.815], False), 'near': (0.0002265378784906, [71.4, 72.95], True)}
    if regime not in expected or spec is None:
        return 'diagnostic_no_new_certificate'
    chi, interval, none = expected[regime]
    same = method == 'robust_shortest' and halfwidth == 1.5 and (target == 0.95) and (q0 == 0.5) and (chi_m_inv == chi) and (set(families) == {'nominal_triangular', 'uniform_support', 'lower_mode_triangular'}) and (spec.get('finite_interval_s') is not None) and (len(spec['finite_interval_s']) == 2) and all((x == y for x, y in zip(spec['finite_interval_s'], interval))) and (spec.get('includes_no_hit_by_horizon') is none)
    return 'prior_frozen_robust_only' if same else 'diagnostic_no_new_certificate'

def run_sensitivity(data_dir, core_results_dir, out_dir):
    import json, time, hashlib, platform, sys, inspect
    from pathlib import Path
    from datetime import datetime, timezone
    import numpy as np
    import pandas as pd
    data_dir = Path(data_dir)
    core_results_dir = Path(core_results_dir)
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise FileExistsError('sensitivity output exists; use a new directory')
    out_dir.mkdir(parents=True)

    def digest(p):
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()

    def dump(p, o):
        Path(p).write_text(json.dumps(o, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    start = time.perf_counter()
    started = datetime.now(timezone.utc).isoformat()
    design = json.loads(json.dumps(SENSITIVITY_DESIGN))
    source_path = Path(__file__)
    engine_source_path = Path(inspect.getsourcefile(run_member))
    source_start_sha = digest(source_path)
    engine_start_sha = digest(engine_source_path)
    bound_inputs = {n: digest(data_dir / n) for n in ['q2_config.json', 'forcing.csv', 'frozen_validation.json']}
    core_paths = {reg: core_results_dir / 'baseline' / ('baseline_model_' + ('Q2_strong_design' if reg == 'strong' else 'Q2_near_boundary_comparison') + '.json') for reg in design['conditions']}
    if any((not p.is_file() for p in core_paths.values())):
        raise FileNotFoundError('complete delivered baseline models required for density-only sensitivity')
    freeze = {
        'utc_before_queries': started,
        'design': design,
        'analysis_source_sha256': source_start_sha,
        'engine_source_sha256': engine_start_sha,
        'engine_source_file': engine_source_path.name,
        'input_sha256': bound_inputs,
        'core_baseline_sha256': {r: digest(p) for r, p in core_paths.items()},
        'raw_scenario_draws': 'none; fixed design grid, no Monte Carlo probability assigned',
        'scientific_boundary': 'certified main candidate endpoints/seed/n/material/threshold unchanged',
    }
    dump(out_dir / 'FROZEN_SENSITIVITY_DESIGN.json', freeze)
    original = load_default_forcing()
    cache = {}
    queries = []
    cache_requests = 0

    def query(q0, chi, d, role):
        nonlocal cache_requests
        q0 = float(q0)
        chi = float(chi)
        d = float(d)
        key = (q0.hex(), chi.hex(), d.hex())
        if key in cache:
            cache_requests += 1
            return cache[key]
        forcing = VelocityForcing(original.t, np.maximum(0.0, original.u + d), identity='fixed_Q2_sensitivity_design', metadata={'offset_m_s': d})
        result = run_member(q0, chi, forcing, start_s=27.0, horizon_s=73.0, trajectory_points=0)
        row = {
            'query_id': len(queries),
            'q0': q0,
            'chi_m_inv': chi,
            'offset_m_s': d,
            'role': role,
            'status': result['status'],
            'tau_s': result['tau_s'],
            'qdot_cross_s_inv': result['qdot_cross_s_inv'],
            'first_hit_valid': result['first_hit_valid'],
            'whole_window_valid': result['whole_window_valid'],
            'solver_status': result['solver_status'],
            'wall_s': result['elapsed_wall_s'],
            'input_path_sha256': forcing.path_sha256,
            'failure': result['failure'],
        }
        hit = result.get('hit_diagnostics')
        if hit is not None and row['status'] == 'hit':
            physics = load_frozen_q2()
            case = physics.Case(**result['case'])
            gamma = case.chi_m_inv * forcing.speed(case, hit['elapsed_s'])
            alpha, N, meq, k = physics.constitutive(case, hit['m'], hit['s_Pa'], gamma)
            H, e, _, Jfactor = physics.closure(case, hit['m'], hit['s_Pa'])
            J = k * Jfactor / case.mu
            shear = gamma * (hit['m'] - meq) / (alpha * case.E)
            drainage = J / H / (alpha * case.E)
            row.update(hit_m=hit['m'], hit_s_Pa=hit['s_Pa'], hit_gamma_s_inv=gamma, shear_qdot_s_inv=shear, drainage_qdot_s_inv=drainage, mechanism_sum_qdot_s_inv=shear + drainage, mechanism_sum_abs_difference=abs(shear + drainage - row['qdot_cross_s_inv']))
        if row['status'] != 'hit' and row['tau_s'] is not None:
            raise RuntimeError('non-hit received finite time')
        cache[key] = row
        queries.append(row)
        with (out_dir / 'q2_queries.jsonl').open('a') as h:
            h.write(json.dumps(row, allow_nan=False) + '\n')
        return row
    phase = []
    for reg, chi in design['conditions'].items():
        for q0 in design['q0_design']:
            for d in np.linspace(-1.5, 1.5, design['offset_grid_count']):
                phase.append({'regime': reg, **query(q0, chi, d, 'q0_offset_phase_grid')})
        print('sensitivity phase map', reg, 'complete', flush=True)
    pd.DataFrame(phase).to_csv(out_dir / 'phase_map.csv', index=False)
    chi_rows = []
    for reg, chi in design['conditions'].items():
        for mult in design['chi_multipliers']:
            for d in design['chi_offset_design_m_s']:
                chi_rows.append({'regime': reg, 'chi_multiplier': mult, **query(0.5, chi * mult, d, 'chi_local_design')})
    pd.DataFrame(chi_rows).to_csv(out_dir / 'chi_perturbation.csv', index=False)
    crossing = []
    derivatives = []
    for reg, chi in design['conditions'].items():
        base = query(0.5, chi, 0.0, 'nominal_crossing')
        slope = base['qdot_cross_s_inv']
        crossing.append({
            'regime': reg,
            'status': base['status'],
            'tau_s': base['tau_s'],
            'qdot_cross_s_inv': slope,
            'seconds_per_0_001_q': None if slope is None or slope == 0 else 0.001 / abs(slope),
            'horizon_slack_s': None if base['status'] != 'hit' else 73.0 - base['tau_s'],
            'positive_0_001_q_linear_delay_exceeds_horizon': None if slope is None or base['status'] != 'hit' else bool(0.001 / abs(slope) > 73.0 - base['tau_s']),
            'shear_qdot_s_inv': base.get('shear_qdot_s_inv'),
            'drainage_qdot_s_inv': base.get('drainage_qdot_s_inv'),
            'mechanism_sum_qdot_s_inv': base.get('mechanism_sum_qdot_s_inv'),
            'mechanism_sum_abs_difference': base.get('mechanism_sum_abs_difference'),
            'drainage_to_abs_shear_ratio': None if base.get('shear_qdot_s_inv') in (None, 0) else base['drainage_qdot_s_inv'] / abs(base['shear_qdot_s_inv']),
            'hit_gamma_s_inv': base.get('hit_gamma_s_inv'),
            'identity': 'first-order transversal response to a hypothetical state error; no state-error distribution assigned',
        })
        for parameter, steps in design['derivative_steps'].items():
            for h in steps:
                if parameter == 'offset_m_s':
                    minus = query(0.5, chi, -h, 'local_derivative')
                    plus = query(0.5, chi, h, 'local_derivative')
                    delta = h
                    center = 0.0
                    unit = 's/(m/s)'
                elif parameter == 'q0':
                    minus = query(0.5 - h, chi, 0.0, 'local_derivative')
                    plus = query(0.5 + h, chi, 0.0, 'local_derivative')
                    delta = h
                    center = 0.5
                    unit = 's per unit q0'
                else:
                    minus = query(0.5, chi * (1 - h), 0.0, 'local_derivative')
                    plus = query(0.5, chi * (1 + h), 0.0, 'local_derivative')
                    delta = chi * h
                    center = chi
                    unit = 's/(m^-1)'
                finite = minus['status'] == 'hit' and plus['status'] == 'hit'
                derivative = (plus['tau_s'] - minus['tau_s']) / (2 * delta) if finite else None
                elasticity = derivative * center / base['tau_s'] if finite and base['status'] == 'hit' and (parameter != 'offset_m_s') else None
                derivatives.append({
                    'regime': reg,
                    'parameter': parameter,
                    'step': h,
                    'absolute_parameter_step': delta,
                    'minus_status': minus['status'],
                    'plus_status': plus['status'],
                    'minus_tau_s': minus['tau_s'],
                    'plus_tau_s': plus['tau_s'],
                    'derivative': derivative,
                    'derivative_unit': unit,
                    'elasticity': elasticity,
                    'parameter_reference_scale': 1.5 if parameter == 'offset_m_s' else center,
                    'scaled_sensitivity': None if not finite or base['status'] != 'hit' else derivative * (1.5 if parameter == 'offset_m_s' else center) / base['tau_s'],
                    'scale_identity': 'offset uses declared1.5m/s halfwidth, q0 uses.5, chi uses nominalchi; not a probability-based global ranking',
                    'method': 'central finite difference, two finite hits only' if finite else 'category transition/unknown; finite derivative not defined',
                    'minus_query_id': minus['query_id'],
                    'plus_query_id': plus['query_id'],
                })
    pd.DataFrame(crossing).to_csv(out_dir / 'crossing_amplification.csv', index=False)
    pd.DataFrame(derivatives).to_csv(out_dir / 'local_derivatives.csv', index=False)
    stability = []
    for reg in design['conditions']:
        for parameter in design['derivative_steps']:
            ds = [r['derivative'] for r in derivatives if r['regime'] == reg and r['parameter'] == parameter and (r['derivative'] is not None)]
            stability.append({
                'regime': reg,
                'parameter': parameter,
                'finite_step_count': len(ds),
                'total_step_count': 3,
                'derivative_min': None if not ds else min(ds),
                'derivative_max': None if not ds else max(ds),
                'relative_step_spread': None if len(ds) < 2 or np.mean(np.abs(ds)) == 0 else (max(ds) - min(ds)) / float(np.mean(np.abs(ds))),
                'identity': 'step variation diagnostic, not a rigorous derivative error bound',
            })
    pd.DataFrame(stability).to_csv(out_dir / 'derivative_step_stability.csv', index=False)
    edges = grid_edges()
    probabilities = {}
    dist_rows = []
    set_rows = []
    reweight_checks = []
    for reg, path in core_paths.items():
        model = BaselineModel.from_json(path)
        probabilities[reg] = {}
        for a in design['support_halfwidths_m_s']:
            propagated = sens_propagate(model, edges, a)
            probabilities[reg][str(a)] = propagated
            masses = np.array([r['finite_mass'] for r in propagated])
            none = np.array([r['no_hit_mass'] for r in propagated])
            targets = design['probability_targets'] if a == 1.5 else [0.95]
            for target in targets:
                nominal = shortest_set(edges, masses[:1], none[:1], target=target)['set']
                robust = shortest_set(edges, masses, none, target=target)['set']
                for method, spec in [('nominal_shortest', nominal), ('robust_shortest', robust)]:
                    set_rows.append({
                        'regime': reg,
                        'support_halfwidth_m_s': a,
                        'target_mass': target,
                        'method': method,
                        'finite_lower_s': None if spec is None or spec['finite_interval_s'] is None else spec['finite_interval_s'][0],
                        'finite_upper_s': None if spec is None or spec['finite_interval_s'] is None else spec['finite_interval_s'][1],
                        'finite_width_s': None if spec is None else spec['finite_width_s'],
                        'includes_no_hit': None if spec is None else spec['includes_no_hit_by_horizon'],
                        'worst_model_mass': None if spec is None else spec['worst_mass'],
                        'independent_certification': sens_certificate_identity(reg, a, target, method, spec, 0.5, design['conditions'][reg], design['families']),
                    })
                for row in propagated:
                    dist_rows.append({
                        'regime': reg,
                        'support_halfwidth_m_s': a,
                        'family': row['family'],
                        'target_mass': target,
                        'hit_mass': float(sum(row['finite_mass'])),
                        'no_hit_mass': row['no_hit_mass'],
                        'unknown_mass': row['unknown_mass'],
                        'total_mass': row['total_mass'],
                    })
            if a == 1.5:
                orig = model.propagated(edges)
                delta = max((float(np.max(np.abs(np.array(r['cdf']) - np.array(s['cdf'])))) for r, s in zip(propagated, orig)))
                reweight_checks.append({'regime': reg, 'original_fullsupport_CDF_max_abs_difference': delta, 'pass': delta < 1e-10})
    pd.DataFrame(dist_rows).to_csv(out_dir / 'distribution_sensitivity.csv', index=False)
    pd.DataFrame(set_rows).to_csv(out_dir / 'set_sensitivity.csv', index=False)
    dump(out_dir / 'probability_sensitivity.json', {'edges_s': edges.tolist(), 'conditions': probabilities})
    reuse_checks = []
    for reg, path in core_paths.items():
        prior = json.loads(path.read_text())['rows']
        for d in [-1.5, 0.0, 1.5]:
            exact = [r for r in prior if abs(float(r['d']) - d) < 1e-14]
            new = query(0.5, design['conditions'][reg], d, 'reuse_check')
            if exact:
                p = exact[0]
                delta = None if new['tau_s'] is None or p['tau_s'] is None else abs(new['tau_s'] - p['tau_s'])
                ok = new['status'] == p['status'] and (new['tau_s'] is None and p['tau_s'] is None or (delta is not None and delta < 1e-05))
                reuse_checks.append({
                    'regime': reg,
                    'offset_m_s': d,
                    'prior_status': p['status'],
                    'current_status': new['status'],
                    'tau_abs_difference_s': delta,
                    'pass': ok,
                })
    checks = {
        'no_nonhit_finite_times': all((r['tau_s'] is None for r in queries if r['status'] != 'hit')),
        'phase_grid_complete': len(phase) == 646,
        'chi_design_complete': len(chi_rows) == 42,
        'mass_closure': all((abs(r['total_mass'] - 1) < 1e-08 for r in dist_rows)),
        'original_fullsupport_reweight_equal': all((r['pass'] for r in reweight_checks)),
        'reused_nominal_reference_equal': all((r['pass'] for r in reuse_checks)),
        'mechanism_qdot_sum_matches': all((r.get('mechanism_sum_abs_difference', 0.0) < 1e-12 for r in queries)),
        'source_unchanged_during_sensitivity': digest(source_path) == source_start_sha and digest(engine_source_path) == engine_start_sha,
    }
    if not all(checks.values()):
        raise RuntimeError('sensitivity self-check failure')
    q = pd.DataFrame(queries)
    q.to_csv(out_dir / 'query_summary.csv', index=False)
    strong = crossing[0]['seconds_per_0_001_q']
    near = crossing[1]['seconds_per_0_001_q']
    phase_states = {reg: {str(k): int(v) for k, v in pd.DataFrame(phase).query('regime==@reg').status.value_counts().items()} for reg in design['conditions']}
    summary = {
        'version': design['version'],
        'conditions': design['conditions'],
        'unique_new_Q2_calls': len(queries),
        'cache_reuse_requests': cache_requests,
        'phase_rows': len(phase),
        'phase_states': phase_states,
        'chi_rows': len(chi_rows),
        'derivative_rows': len(derivatives),
        'crossing_amplification': crossing,
        'near_to_strong_amplification_ratio': near / strong,
        'nominal_horizon_slack_s': {r['regime']: r['horizon_slack_s'] for r in crossing},
        'linearization_boundary': 'state-error time amplification applies to transverse local roots only; if delay exceeds horizon slack retain category transition, not a finite in-window confidence band',
        'derivative_step_stability': stability,
        'checks': checks,
        'reweight_checks': reweight_checks,
        'nominal_reference_checks': reuse_checks,
        'independent_review': 'pending; execution self-check only',
        'probability_boundary': design['certification_rule'],
        'interpretation': 'q0/chi/offset grid is unweighted, not event probability; finite derivatives only where two sides hit',
    }
    dump(out_dir / 'summary.json', summary)
    receipt = {
        'version': design['version'],
        'started_utc': started,
        'completed_utc': datetime.now(timezone.utc).isoformat(),
        'wall_s': time.perf_counter() - start,
        'python': platform.python_version(),
        'executable': sys.executable,
        'analysis_source_sha256': source_start_sha,
        'engine_source_sha256': engine_start_sha,
        'input_sha256': bound_inputs,
        'core_baseline_sha256': freeze['core_baseline_sha256'],
        'unique_new_Q2_calls': len(queries),
        'cache_reuse_requests': cache_requests,
        'Q2_solver_wall_sum_s': float(q.wall_s.sum()),
        'density_reweight_extra_Q2_calls': 0,
        'outputs': {p.name: digest(p) for p in out_dir.iterdir() if p.is_file()},
        'checks': checks,
        'independent_review': 'pending',
    }
    dump(out_dir / 'SENSITIVITY_RECEIPT.json', receipt)
    print(json.dumps({
        'sensitivity_output': str(out_dir),
        'unique_Q2_calls': len(queries),
        'wall_s': receipt['wall_s'],
        'all_execution_checks_pass': True,
    }), flush=True)
    return summary
RESULT_VERSION = 'Q3-STRENGTHENED-v4.0'

def integrated_sensitivity(core_results, out_dir):
    answer = run_sensitivity(DATA, core_results, out_dir)
    required = ['legacy_unweighted_reference.csv', 'legacy_numerical_convergence.csv', 'LEGACY_SENSITIVITY_SOURCES.json']
    for name in required:
        if not (DATA / name).is_file():
            raise FileNotFoundError('required sensitivity reference missing: data/' + name)
        (out_dir / name).write_bytes((DATA / name).read_bytes())
    legacy = pd.read_csv(DATA / 'legacy_unweighted_reference.csv')
    phase = pd.read_csv(out_dir / 'phase_map.csv')
    aliases = {'Q2_strong_design': 'strong', 'Q2_near_boundary_comparison': 'near'}
    checks = []
    for _, row in legacy.iterrows():
        if row.regime not in aliases:
            continue
        sub = phase[(phase.regime == aliases[row.regime]) & np.isclose(phase.q0, float(row.q0), rtol=0, atol=1e-14) & np.isclose(phase.offset_m_s, float(row.offset_m_s), rtol=0, atol=1e-14)]
        if len(sub) != 1:
            raise RuntimeError('old reference cannot be matched to unique new phase input')
        new = sub.iloc[0]
        delta = None if row.status != 'hit' or new.status != 'hit' else abs(float(row.first_hit_s) - float(new.tau_s))
        ok = bool(row.status == new.status and (row.status != 'hit' or (delta is not None and delta < 1e-05)))
        checks.append({
            'regime': aliases[row.regime],
            'q0': float(row.q0),
            'offset_m_s': float(row.offset_m_s),
            'old_status': row.status,
            'new_status': new.status,
            'tau_abs_difference_s': delta,
            'pass': ok,
        })
    if len(checks) != 24 or not all((x['pass'] for x in checks)):
        raise RuntimeError('legacy24 scenario replay mismatch')
    oldnum = pd.read_csv(DATA / 'legacy_numerical_convergence.csv')
    if len(oldnum) != 6 or not oldnum.status_agrees.all():
        raise ValueError('legacy numerical reference malformed')
    oldchecks = {
        'role': 'historical numerical tolerance tests, not new physical sensitivity or new certification',
        'case_count': len(oldnum),
        'max_finite_tau_difference_s': float(oldnum.tau_abs_difference_s.max()),
        'all_statuses_equal': bool(oldnum.status_agrees.all()),
    }
    save(out_dir / 'LEGACY_REUSE_CHECKS.json', {
        'matched24_scenarios': checks,
        'old4_initial_admission_controls': 'preserved in source table, outside this phase grid',
        'old_numerical_controls': oldchecks,
    })
    result = {
        'version': RESULT_VERSION,
        'actual_sensitivity_receipt_sha256': sha(out_dir / 'SENSITIVITY_RECEIPT.json'),
        'analysis_sha256': sha(__file__),
        'legacy_input_sha256': {n: sha(DATA / n) for n in required},
        'new_unique_Q2_calls': answer['unique_new_Q2_calls'],
        'legacy_reference_only_new_Q2_calls': 0,
        'legacy_24_matching_checks_passed': True,
        'metadata_identity': 'new phase/chi/derivative sensitivity plus explicit historical numerical references; only original frozen robust objects retain old statistical evidence',
        'output_sha256': {str(p.relative_to(out_dir)): sha(p) for p in out_dir.rglob('*') if p.is_file()},
    }
    save(out_dir / 'INTEGRATED_SENSITIVITY_RECEIPT.json', result)
    return answer

def sensitivity_only_run(results, core_results, started, wall, source_sha):
    if core_results is None:
        raise ValueError('--stage sensitivity requires an explicit delivered --core-results directory')
    core_results = core_results.resolve()
    answer = integrated_sensitivity(core_results, results / 'sensitivity')
    summary = {
        'version': RESULT_VERSION,
        'execution_stage': 'sensitivity',
        'sensitivity': answer,
        'core_results_supplied_explicitly': str(core_results),
        'base_certificate': 'not rerun by this supplementary stage; no new scenarios inherit it',
        'independent_review': 'pending',
        'analysis_source_sha256': source_sha,
    }
    save(results / 'summary.json', summary)
    if sha(__file__) != source_sha:
        raise RuntimeError('source changed during supplementary run')
    save(results / 'RUN_RECEIPT.json', {
        'version': RESULT_VERSION,
        'stage': 'sensitivity',
        'started_utc': started,
        'completed_utc': stamp(),
        'wall_s': time.perf_counter() - wall,
        'command_arguments': sys.argv[1:],
        'python': platform.python_version(),
        'executable': sys.executable,
        'source_sha256': source_sha,
        'versions': {
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'pandas': pd.__version__,
            'sklearn': sklearn.__version__,
        },
        'input_sha256': {p.name: sha(p) for p in DATA.iterdir() if p.is_file()},
        'output_sha256': {str(p.relative_to(results)): sha(p) for p in results.rglob('*') if p.is_file()},
        'scientific_acceptance': 'pending independent review; supplementary stage is not a complete base workflow reproduction',
    })
    print(json.dumps({'results': str(results), 'stage': 'sensitivity', 'unique_new_Q2_calls': answer['unique_new_Q2_calls']}), flush=True)

def sol():
    global DATA, OUT
    ap = argparse.ArgumentParser(description='Q3 analysis')
    ap.add_argument('--data', type=Path, default=PACKAGE_ROOT / 'data')
    ap.add_argument('--results', type=Path, default=PACKAGE_ROOT / 'results')
    ap.add_argument('--stage', choices=['all', 'develop', 'sensitivity'], default='all')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--q3-entry', type=Path)
    ap.add_argument('--q4-arrivals', type=Path)
    ap.add_argument('--core-results', type=Path)
    args = ap.parse_args()
    if not 1 <= args.workers <= 4:
        raise ValueError('workers must be1..4')
    configure(args.data.resolve())
    results = args.results.resolve()
    if results.exists():
        raise FileExistsError('results directory exists: use a new output directory')
    results.mkdir(parents=True)
    started = stamp()
    wall = time.perf_counter()
    before = sha(__file__)
    save(results / 'RUN_CONFIG.json', {'final': FINAL, 'gp': CONFIG, 'stage': args.stage, 'workers': args.workers})
    if args.stage == 'sensitivity':
        sensitivity_only_run(results, args.core_results, started, wall, before)
        return
    if args.core_results is not None:
        raise ValueError('--core-results is only for the supplementary sensitivity stage')
    execution_controls(results)
    cost = numerical_pipeline(results)
    validation = validate_final(results, args.workers) if args.stage == 'all' else {'status': 'pending_development_only', 'claims': [], 'all_passed': False}
    if args.stage != 'all':
        save(results / 'validation.json', validation)
    lead_outputs(results, validation)
    q4 = json.loads((DATA / 'q4_context.json').read_text()) if (DATA / 'q4_context.json').exists() else {'status': 'pending', 'reason': 'q4_context.json not delivered'}
    if bool(args.q3_entry) != bool(args.q4_arrivals):
        raise ValueError('both q3-entry and q4-arrivals are required for connection')
    if 'existing_five_point_ODE_demo' in q4:
        q4_outputs(results, q4)
    if args.q3_entry:
        linked = connect(json.loads(args.q3_entry.read_text()), json.loads(args.q4_arrivals.read_text()))
        save(results / 'q4_connected.json', linked)
    sensitivity = integrated_sensitivity(results, results / 'sensitivity')
    summary = {
        'version': RESULT_VERSION,
        'base_certification_version': FINAL['version'],
        'sensitivity': sensitivity,
        'started_utc': started,
        'completed_utc': stamp(),
        'wall_s': time.perf_counter() - wall,
        'execution_stage': args.stage,
        'random_validation': validation,
        'conditions': CONFIG['conditions'],
        'candidate_sets': FINAL['candidates'],
        'cost': cost,
        'q4_summary': q4,
        'method_choice': 'adaptive_piecewise_integration_primary; GP actual comparison and optional accelerator',
        'independent_review': 'pending',
        'natural_event_coverage': 'not identified',
        'self_contained_source_sha256': before,
    }
    save(results / 'summary.json', summary)
    if sha(__file__) != before:
        raise RuntimeError('source changed during run')
    receipt = {
        'version': RESULT_VERSION,
        'base_certification_version': FINAL['version'],
        'started_utc': started,
        'completed_utc': stamp(),
        'wall_s': time.perf_counter() - wall,
        'command_arguments': sys.argv[1:],
        'python': platform.python_version(),
        'executable': sys.executable,
        'versions': {
            'numpy': np.__version__,
            'scipy': scipy.__version__,
            'pandas': pd.__version__,
            'sklearn': sklearn.__version__,
        },
        'source_sha256': before,
        'input_sha256': {p.name: sha(p) for p in DATA.iterdir() if p.is_file()},
        'output_sha256': {str(p.relative_to(results)): sha(p) for p in results.rglob('*') if p.is_file()},
        'scientific_acceptance': 'pending independent review; execution receipt is not acceptance',
    }
    save(results / 'RUN_RECEIPT.json', receipt)
    print(json.dumps({
        'results': str(results),
        'validation': validation['status'],
        'all_claims_pass': validation.get('all_passed'),
        'wall_s': receipt['wall_s'],
    }), flush=True)
if __name__ == '__main__':
    sol()

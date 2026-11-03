from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass, replace, asdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.integrate import solve_ivp, cumulative_simpson, simpson
from scipy.optimize import brentq

G = 9.81


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


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
        return self.initial_stress_ratio*self.E

    @property
    def qcrit(self):
        return (self.shape_n-1)/(2*(self.shape_n+1))

    @property
    def ki(self):
        return self.Kfs_transfer_m_s*self.water_conversion_mu_Pa_s/(self.rho_f*G)

    @property
    def Hs(self):
        return self.m0*self.H0

    @property
    def E(self):
        return (self.rho_s-self.rho_f)*G*math.cos(math.radians(self.theta_deg))*self.Hs


class Forcing:
    def __init__(self, path):
        d = pd.read_csv(path)
        self.t = d.time_since_refined_onset_s.to_numpy(float)
        self.u = d.velocity_m_s.to_numpy(float)
        if not (np.all(np.diff(self.t)>0) and np.all(self.u>=0)):
            raise ValueError('Forcing requires ordered times and nonnegative velocities')

    def speed(self, case, t):
        if case.start < self.t[0] or case.end > self.t[-1]:
            raise ValueError('Requested forcing window exceeds digitized data; no extrapolation')
        return max(0.0, float(np.interp(t+case.start, self.t, self.u))+case.curve_shift)


def constitutive(case, m, s, gamma):
    if not (0<m<1) or s <= -case.sigma0:
        raise ValueError('Constitutive state outside domain')
    alpha = case.a/(m*(s+case.sigma0))
    N = 0.0 if gamma == 0 else case.mu*gamma/(case.rho_s*gamma**2*case.delta**2+max(s,0.0))
    meq = case.mcrit/(1+math.sqrt(N))
    k = case.ki*math.exp((case.m0-m)/0.04)
    return alpha, N, meq, k


def closure(case, m, s):
    H = case.Hs/m
    zeta = s/case.E
    if zeta >= case.qcrit:
        e = (case.shape_n+1)/case.shape_n*(case.E/2-s)
        r = 1.0
        Jfactor = case.shape_n*e/H
    else:


        r = math.sqrt(max(0.0,zeta/case.qcrit))
        e = case.E*(1-(case.shape_n-1)*r/case.shape_n)
        Jfactor = case.E/H
    return H, e, r, Jfactor


def rhs(case, forcing, t, y):
    m,s,_=y
    U=forcing.speed(case,t)
    gamma=0.0 if case.mode=='zero_shear' else case.chi_m_inv*U
    alpha,N,meq,k=constitutive(case,m,s,gamma)
    H,e,r,Jfactor=closure(case,m,s)
    if case.mode=='drained':
        dm=-m*gamma*(m-meq)
        J=H*dm/m
        ds=0.0
    else:
        J=k*Jfactor/case.mu
        dm=m*J/H
        ds=(gamma*(m-meq)+J/H)/alpha
    return (dm,ds,J)


def profile(case,m,s,x):
    H,e,r,_=closure(case,m,s)
    if s/case.E >= case.qcrit:
        pe=e*(1-x**case.shape_n)
        effective=case.E*(1-x)-pe
    else:
        rn=max(r,1e-100)
        pe=np.where(x<=r,e-case.E*x**case.shape_n/(case.shape_n*rn**(case.shape_n-1)),case.E*(1-x))
        effective=np.where(x<=r,case.E*((case.shape_n-1)*r/case.shape_n-x+x**case.shape_n/(case.shape_n*rn**(case.shape_n-1))),0.0)
    hydro=case.rho_f*G*math.cos(math.radians(case.theta_deg))*H*(1-x)
    return pe, effective, hydro+pe


def initial_classification(case):
    q0=case.initial_stress_ratio
    B=case.rho_f*G*math.cos(math.radians(case.theta_deg))*case.H0
    qmax=.5+B/((case.shape_n+1)*case.E)
    if not math.isfinite(q0) or q0<0 or q0>qmax+1e-12:
        status='invalid_initial_state'
    elif math.isclose(q0,0,rel_tol=0,abs_tol=1e-12):
        status='initial_closure_limit'
    elif math.isclose(q0,qmax,rel_tol=0,abs_tol=1e-12):
        status='initial_saturation_boundary'
    elif q0<case.qcrit-1e-12:
        status='initial_positive_thickness_layer'
    elif math.isclose(q0,case.qcrit,rel_tol=0,abs_tol=1e-12):
        status='initial_first_layer_boundary'
    else:
        status='admissible_unreached_initial_state'
    f0=1-math.sqrt(q0/case.qcrit) if 0<=q0<case.qcrit else (0.0 if q0>=case.qcrit else None)
    if status=='invalid_initial_state':f0=None
    return {'q0':q0,'s_init_Pa':case.s_init,'qcrit':case.qcrit,'qmax':qmax,
            'initial_status':status,'initial_fluid_supported_fraction':f0,
            'initial_condition_satisfied':None if status in ('initial_closure_limit','invalid_initial_state') else status in ('initial_positive_thickness_layer','initial_first_layer_boundary'),
            'dynamic_first_hit_elapsed_s':None,
            'identity':'mathematical closure admission only; not local observed or credible range'}


def run_case(case, forcing, numerics, dense=True, keep_solution=False):
    if not (0<case.m0<1 and 0<case.mcrit<1 and case.H0>0 and case.ki>0 and case.mu>0 and case.shape_n>1):
        raise ValueError('Invalid physical input')
    admission=initial_classification(case)
    if admission['initial_status']!='admissible_unreached_initial_state':
        raise ValueError('Initial state is classified separately: '+admission['initial_status'])
    if case.mode=='drained' and case.initial_stress_ratio!=.5:
        raise ValueError('Historical ideal-drained control is restricted to hydrostatic initial state')
    def zero_stress(t,y): return y[1]
    zero_stress.terminal=True
    zero_stress.direction=-1
    def desaturation(t,y):
        H,e,r,_=closure(case,y[0],y[1])
        return e+case.rho_f*G*math.cos(math.radians(case.theta_deg))*H/case.shape_n
    desaturation.terminal=True
    desaturation.direction=-1
    def fluid_layer(t,y): return y[1]/case.E-case.qcrit
    fluid_layer.terminal=False
    fluid_layer.direction=-1
    duration=case.end-case.start
    sol=solve_ivp(lambda t,y:rhs(case,forcing,t,y),(0,duration),(case.m0,case.s_init,0.0),
                  method='Radau',rtol=numerics['rtol'],atol=(1e-11,1e-6,1e-10),
                  max_step=numerics['max_step_s'],dense_output=True,
                  events=(zero_stress,desaturation,fluid_layer))
    if not sol.success: raise RuntimeError(sol.message)
    state='completed_conditional_window'
    if len(sol.t_events[0]): state='mean_stress_zero_closure_limit'
    if len(sol.t_events[1]): state='saturation_hypothesis_boundary'
    if dense:
        ts=np.unique(np.r_[np.arange(0,sol.t[-1],numerics['output_step_s']),sol.t[-1],sol.t_events[2]])
    else:
        ts=np.unique(np.r_[sol.t,sol.t_events[2]])
    ys=sol.sol(ts)
    rows=[]
    residual=[]
    for t,(m,s,budget) in zip(ts,ys.T):
        H,e,r,_=closure(case,m,s)
        U=forcing.speed(case,t)
        gamma=0.0 if case.mode=='zero_shear' else case.chi_m_inv*U
        alpha,N,meq,k=constitutive(case,m,s,gamma)
        dm,ds,J=rhs(case,forcing,t,(m,s,budget))
        if case.mode=='drained': e=0.0; r=1.0
        resid=gamma*(m-meq)-alpha*ds+dm/m
        residual.append(resid)
        rows.append({'elapsed_s':t,'source_relative_s':t+case.start,'U_input_m_s':U,
                     'gamma_input_s_inv':gamma,'m':m,'H_m':H,'mean_effective_stress_Pa':s,
                     'mean_effective_ratio':s/case.E,'basal_excess_Pa':e,'basal_excess_ratio':e/case.E,
                     'fluid_supported_fraction':1-r,'drainage_flux_m_s':J,'cumulative_outflow_m':budget,
                     'skeleton_thickness_m':H*r,'skeleton_to_delta_ratio':H*r/case.delta,
                     'alpha_Pa_inv':alpha,'k_m2':k,'N':N,'meq':meq,
                     'solid_volume_residual_m':m*H-case.Hs,
                     'water_budget_residual_m':H-case.H0+budget,'constitutive_residual_s_inv':resid})
    frame=pd.DataFrame(rows)
    x=np.linspace(0,1,1001)
    profile_checks=[]
    for j in np.unique(np.linspace(0,len(frame)-1,min(21,len(frame))).astype(int)):
        row=frame.iloc[j]
        pe,eff,total=profile(case,row.m,row.mean_effective_stress_Pa,x)
        if case.mode=='drained':
            eff=case.E*(1-x);total=case.rho_f*G*math.cos(math.radians(case.theta_deg))*row.H_m*(1-x)
        profile_checks.append((np.min(eff),np.min(total),abs(np.trapezoid(eff,x)-row.mean_effective_stress_Pa)/case.E))
    onset=float(sol.t_events[2][0]) if len(sol.t_events[2]) else None
    summary={'q0':case.initial_stress_ratio,'s_init_Pa':case.s_init,'initial_status':admission['initial_status'],'status':state,'onset_elapsed_s':onset,'onset_source_relative_s':None if onset is None else onset+case.start,
             'max_fluid_fraction':float(frame.fluid_supported_fraction.max()),
             'min_mean_effective_ratio':float(frame.mean_effective_ratio.min()),
             'final_m':float(frame.m.iloc[-1]),'final_H_m':float(frame.H_m.iloc[-1]),
             'end_elapsed_s':float(sol.t[-1]),'ki_m2':case.ki,'E_Pa':case.E,
             'max_solid_residual_m':float(frame.solid_volume_residual_m.abs().max()),
             'max_water_residual_m':float(frame.water_budget_residual_m.abs().max()),
             'max_constitutive_residual_s_inv':float(max(abs(x) for x in residual)),
             'min_profile_effective_stress_Pa':float(min(x[0] for x in profile_checks)),
             'min_profile_total_pore_Pa':float(min(x[1] for x in profile_checks)),
             'max_profile_mean_relative_error':float(max(x[2] for x in profile_checks)),
             'min_skeleton_to_delta_ratio':float(frame.skeleton_to_delta_ratio.min()),
             'continuum_scale_warning':bool(frame.skeleton_to_delta_ratio.min()<10),
             'continuum_warning_identity':'10 delta is a disclosed resolution diagnostic, not a calibrated formation boundary'}
    if keep_solution: return frame,summary,sol
    return frame,summary


def threshold(case,forcing,cfg,numerics):
    admission=initial_classification(case)
    if admission['initial_status']!='admissible_unreached_initial_state':
        return {'threshold_status':admission['initial_status'],'chi_critical_m_inv':None,'critical_identity':'no new dynamic threshold for classified initial states'},[]
    search=cfg['threshold_search']

    maxu=float(max(forcing.u))+max(case.curve_shift,0)
    lo=search['chi_min_m_inv']; hi=search['max_gamma_s_inv']/maxu
    grid=np.geomspace(lo,hi,search['bracket_points'])
    values=[];records=[]
    def objective(chi):
        _,s=run_case(replace(case,chi_m_inv=float(chi)),forcing,numerics,False)
        if s['status']!='completed_conditional_window':
            return None,s
        return s['min_mean_effective_ratio']-case.qcrit,s
    for chi in grid:
        try:
            val,s=objective(chi)
        except (ValueError,RuntimeError) as error:
            val=None;s={'status':'numerical_or_physical_domain_failure','end_elapsed_s':None,'onset_elapsed_s':None,'failure':str(error)}
        values.append(val)
        records.append({'record_type':'search_grid','chi_m_inv':float(chi),'margin_mean_stress_ratio':val,
                        'status':s['status'],'end_elapsed_s':s['end_elapsed_s'],'onset_elapsed_s':s['onset_elapsed_s']})
    roots=[];root_records=[]
    for i in range(len(grid)-1):
        if values[i] is not None and values[i+1] is not None and values[i]*values[i+1]<0:
            def valid_objective(chi):
                val,s=objective(chi)
                if val is None: raise ValueError('Threshold bracket entered an incomplete/invalid forcing window')
                return val
            root=brentq(valid_objective,grid[i],grid[i+1],xtol=search['xtol_m_inv'])
            roots.append(float(root))
            eps=max(search['xtol_m_inv']*10,root*1e-5)
            left,ls=objective(max(grid[i],root-eps));right,rs=objective(min(grid[i+1],root+eps))
            root_records.append({'record_type':'valid_root_bracket','chi_m_inv':float(root),
                                 'bracket_lo_m_inv':float(grid[i]),'bracket_hi_m_inv':float(grid[i+1]),
                                 'bracket_margin_lo':values[i],'bracket_margin_hi':values[i+1],
                                 'neighbor_left_margin':left,'neighbor_right_margin':right,
                                 'neighbor_left_chi_m_inv':float(max(grid[i],root-eps)),
                                 'neighbor_right_chi_m_inv':float(min(grid[i+1],root+eps)),
                                 'neighbor_left_status':ls['status'],'neighbor_right_status':rs['status'],
                                 'root_event':'first-layer stress boundary within complete prescribed window'})
    valid_values=[v for v in values if v is not None]
    invalid_count=len(values)-len(valid_values)
    monotone=bool(np.all(np.diff(valid_values)<=1e-5)) if len(valid_values)>1 else None
    if roots: threshold_state='root_found_in_valid_bracket'
    elif invalid_count: threshold_state='domain_exit_or_failure_prevents_complete_window_search'
    elif valid_values and min(valid_values)>0: threshold_state='no_onset_at_valid_sampled_chi_values'
    else: threshold_state='onset_at_search_lower_bound_or_no_crossing_found'
    out={'search_chi_min_m_inv':float(lo),'search_chi_max_m_inv':float(hi),
         'sampled_margin_nonincreasing':monotone,'crossings_found':len(roots),
         'chi_critical_m_inv':roots[0] if roots else None,
         'equivalent_length_at_critical_m':1/roots[0] if roots else None,
         'valid_complete_window_grid_points':len(valid_values),'incomplete_or_failed_grid_points':invalid_count,
         'critical_identity':'first valid sampled bracket crossing; finite grid does not exclude missed tangencies or all other nonmonotone intervals',
         'threshold_status':threshold_state}
    return out,records+root_records


def speed_integral(forcing, start, end):
    if start<forcing.t[0] or end>forcing.t[-1] or end<start:
        raise ValueError('Integral outside actual curve support')
    t=np.unique(np.r_[start,forcing.t[(forcing.t>start)&(forcing.t<end)],end])
    return float(np.trapezoid(np.interp(t,forcing.t,forcing.u),t))


def diagnostic_grid(case,forcing,cut,max_step):
    nodes=np.unique(np.r_[0,forcing.t[(forcing.t>case.start)&(forcing.t<case.start+cut)]-case.start,cut])
    return np.unique(np.concatenate([np.linspace(a,b,max(2,int(math.ceil((b-a)/max_step))+1)) for a,b in zip(nodes[:-1],nodes[1:])]))


def shear_diagnostics(case,forcing,sol,state,case_id,regime,max_step):
    hit=state['onset_elapsed_s'];cut=state['end_elapsed_s'] if hit is None else hit
    ts=diagnostic_grid(case,forcing,cut,max_step);ys=sol.sol(ts)
    m=ys[0];s=ys[1];gamma=np.asarray([(0.0 if case.mode=='zero_shear' else case.chi_m_inv*forcing.speed(case,t)) for t in ts])
    N=[];psi=[];jh=[]
    for t,mi,si,ga in zip(ts,m,s,gamma):
        alpha,ni,meq,k=constitutive(case,mi,si,ga);H,e,r,jfac=closure(case,mi,si)
        N.append(ni);psi.append(mi-meq);jh.append(k*jfac/(case.mu*H))
    N=np.asarray(N);psi=np.asarray(psi);jh=np.asarray(jh)
    shear=m*gamma*psi/case.a;drain=m*jh/case.a
    Is=cumulative_simpson(shear,x=ts,initial=0);Id=cumulative_simpson(drain,x=ts,initial=0)
    lhs=np.log((s+case.sigma0)/(case.s_init+case.sigma0));volume=(m-case.m0)/case.a
    residual=lhs-Is-volume
    Delta=case.mcrit-case.m0
    IU=speed_integral(forcing,case.start,case.end)
    Gamma_asy=case.a/(case.m0*Delta)*math.log((case.s_init+case.sigma0)/(case.qcrit*case.E+case.sigma0)) if Delta>0 else None
    Lgoal=math.log((case.s_init+case.sigma0)/(case.qcrit*case.E+case.sigma0)) if case.initial_stress_ratio>case.qcrit else None
    Ghit=None if hit is None else case.chi_m_inv*speed_integral(forcing,case.start,case.start+hit)
    def finite_or_none(value): return float(value) if np.isfinite(value) else None
    eps_m=float(np.max(np.abs(m-case.m0))/Delta) if Delta>0 else None
    eps_N=float(np.max(case.mcrit*np.sqrt(N)/((1+np.sqrt(N))*Delta))) if Delta>0 else None
    record={'case_id':case_id,'regime':regime,'shape_n':case.shape_n,'qcrit':case.qcrit,
            'criterion_origin':'closure_derived','hypothesis_id':'H_SATURATED_MIXED_MATRIX_CHI_ASL',
            'forcing_origin':'post_event_inversion','forcing_window':f'[{case.start:g},{case.end:g}]',
            'q0':case.initial_stress_ratio,'s_init_Pa':case.s_init,'m0':case.m0,'mcrit':case.mcrit,'H0_m':case.H0,'chi_trajectory_m_inv':case.chi_m_inv,
            'status':state['status'],'diagnostic_cut_elapsed_s':float(cut),
            'diagnostic_cut_identity':'first_hit' if hit is not None else ('legal_window_end' if state['status']=='completed_conditional_window' else 'hypothesis_exit'),
            'I_U_window_m':IU,'Gamma_window_prescribed':case.chi_m_inv*IU,
            'Gamma_until_actual_stop':case.chi_m_inv*speed_integral(forcing,case.start,case.start+state['end_elapsed_s']),
            'Gamma_at_first_hit':Ghit,'Gamma_c_asy':Gamma_asy,'chi_asy_m_inv':None if Gamma_asy is None else Gamma_asy/IU,
            'Gamma_hit_relative_error_vs_asy':None if Ghit is None or Gamma_asy is None else (Ghit-Gamma_asy)/Gamma_asy,
            'max_N_before_cut':float(np.max(N)),'eps_m':eps_m,'eps_N':eps_N,
            'signed_shear_log_contribution':float(Is[-1]),'signed_volume_log_contribution':float(volume[-1]),
            'direct_drain_log_integral':float(Id[-1]),'absolute_shear_log_integral':float(simpson(np.abs(shear),x=ts)),
            'absolute_drain_log_integral':float(simpson(np.abs(drain),x=ts)),
            'absolute_gamma_psi_integral':float(simpson(np.abs(gamma*psi),x=ts)),
            'absolute_J_over_H_integral':float(simpson(np.abs(jh),x=ts)),
            'max_identity_residual':float(np.max(np.abs(residual))),
            'final_identity_residual':float(residual[-1]),'drain_quadrature_minus_delta_m_over_a':float(Id[-1]-volume[-1]),
            'L_goal':Lgoal,'eps_J':None if Lgoal is None else float(simpson(np.abs(drain),x=ts)/Lgoal),
            'absolute_drain_to_shear_ratio':None if simpson(np.abs(shear),x=ts)==0 else float(simpson(np.abs(drain),x=ts)/simpson(np.abs(shear),x=ts)),
            'diagnostic_samples':len(ts),'diagnostic_max_step_s':float(np.max(np.diff(ts))),
            'approximation_status':'not_applicable_Delta_nonpositive' if Delta<=0 else ('pre_layer_comparison_no_prespecified_success_cutoff' if hit is not None else 'diagnostics_only_no_first_hit'),
            'approximation_reason':'requires Delta>0, weak drainage, nearly fixed m and N correction small relative to Delta; report diagnostics without a pass percentage'}
    Gamma=np.asarray([case.chi_m_inv*speed_integral(forcing,case.start,case.start+t) for t in ts])
    sasy=np.full_like(s,np.nan) if Delta<=0 else (case.s_init+case.sigma0)*np.exp(-case.m0*Delta*Gamma/case.a)-case.sigma0
    frame=pd.DataFrame({'case_id':case_id,'regime':regime,'elapsed_s':ts,'m':m,'s_Pa':s,'s_ratio':s/case.E,
                        'gamma_s_inv':gamma,'N':N,'tan_psi':psi,'J_over_H_s_inv':jh,
                        'Gamma_sh':Gamma,'s_asy_ratio':sasy/case.E,'shear_log_integrand':shear,
                        'drain_log_integrand':drain,'shear_log_integral':Is,'direct_drain_log_integral':Id,
                        'volume_log_term':volume,'log_stress_ratio':lhs,'identity_residual':residual})
    return record,frame


def reached_state(state):
    hit=state['onset_elapsed_s'] is not None
    if state['status']=='completed_conditional_window':
        return 'later_first_hit' if hit else 'valid_no_hit_in_window'
    return 'first_hit_then_domain_exit' if hit else 'domain_exit_before_first_hit'


def sol():
    parser=argparse.ArgumentParser(description='Q2 v3 initial-effective-stress / cumulative-shear joint conditional boundary')
    parser.add_argument('--delivery-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args();started=datetime.now(timezone.utc).isoformat()
    delivery=args.delivery_dir.resolve();out=(args.output_dir or delivery/'results').resolve()
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((delivery/'config.json').read_text());en=cfg['initial_state_enhancement'];num=cfg['numerics']
    forcing=Forcing(delivery/cfg['forcing']['source_file']);ref=Case(**cfg['reference'])
    baseline_cfg=json.loads((delivery/'data/v2_reference_config.json').read_text())
    baseline_binding=json.loads((delivery/'data/BASELINE_BINDING.json').read_text())
    baseline_hashes={r['path'].split('/',1)[1]:r['sha256'] for r in baseline_binding['protected_checks']
                     if r['path'].startswith('Q2_active_data_20261029_v2/')}
    if fingerprint(delivery/'data/v2_reference_config.json')!=baseline_hashes['config.json']:
        raise RuntimeError('Inherited reference configuration does not match reviewed v2 binding')
    if ref!=Case(**baseline_cfg['reference']):
        raise ValueError('Frozen v3 reference differs from reviewed v2; inherited outputs require re-audit')
    if fingerprint(delivery/cfg['forcing']['source_file'])!=baseline_hashes[cfg['forcing']['source_file']]:
        raise ValueError('Forcing differs from reviewed v2; historical results cannot be silently reused')
    import shutil
    reuse=delivery/en['v2_reused_results']

    reused=[]
    for source in sorted(reuse.iterdir()):
        if source.is_file() and source.name not in ('ANALYSIS_RUN_RECEIPT.json','summary.json'):
            if fingerprint(source)!=baseline_hashes['results/'+source.name]:
                raise RuntimeError('Historical result reuse input has changed: '+source.name)
            shutil.copyfile(source,out/source.name);reused.append(source.name)
    oldcum=pd.read_csv(reuse/'cumulative_shear_check.csv')
    oldsum=pd.read_csv(reuse/'enhancement_case_summary.csv')
    oldbr=pd.read_csv(reuse/'enhancement_search_brackets.csv')
    reference_old=oldsum[(oldsum.shape_n==2)&(oldsum.m0==.54)&(oldsum.mcrit==.60)&(oldsum.T_s==73)].iloc[0]
    summaries=[];trajectories=[];diagnostics=[];series=[];bounds=[];brackets=[];regressions=[];relax=[]
    IU=speed_integral(forcing,ref.start,ref.end)
    for i,q0 in enumerate(en['q0_values'],1):
        case=replace(ref,initial_stress_ratio=float(q0));cid=f'V3_Q{i:02d}'
        admission=initial_classification(case)
        if admission['initial_status']!='admissible_unreached_initial_state':
            raise ValueError('Frozen main matrix contains non-interior initial state')
        for regime,chi in [('common_design',en['common_chi_design_m_inv']),('common_comparison',en['common_chi_comparison_m_inv'])]:
            c=replace(case,chi_m_inv=float(chi))
            f,state,sol=run_case(c,forcing,num,keep_solution=True)
            f.insert(0,'case_id',cid);f['regime']=regime;f['q0']=q0;trajectories.append(f)
            rec,ser=shear_diagnostics(c,forcing,sol,state,cid,regime,en['diagnostic_max_step_s'])
            diagnostics.append(rec);series.append(ser)
            H,e,r,jfactor=closure(c,c.m0,c.s_init)
            gamma0=chi*forcing.speed(c,0)
            alpha,N,meq,k=constitutive(c,c.m0,c.s_init,gamma0)
            J=k*jfactor/c.mu
            dm,ds,_=rhs(c,forcing,0,(c.m0,c.s_init,0))
            summaries.append({'case_id':cid,'regime':regime,'shape_n':c.shape_n,
                'closure_id':'no_tension_power_n2','criterion_origin':'closure_derived',
                'hypothesis_id':cfg['hypothesis_id'],'forcing_origin':'post_event_inversion',
                'initial_state_origin':en['initial_state_identity'],'forcing_window':'[27,100]',
                'q0':q0,'s_init_Pa':c.s_init,'qmax':admission['qmax'],'qcrit':c.qcrit,
                'm0':c.m0,'mcrit':c.mcrit,'H0_m':c.H0,'chi_m_inv':chi,
                'T_s':c.end-c.start,'I_U_m':IU,'initial_basal_excess_Pa':e,'J0_m_s':J,
                'm_dot0_s_inv':dm,'s_dot0_Pa_s':ds,'initial_s_dot_sign':int(np.sign(ds)),
                'gamma0_s_inv':gamma0,'N0':N,'tan_psi0':c.m0-meq,
                'initial_shear_driver_s_inv':gamma0*(c.m0-meq),'initial_J_over_H_s_inv':J/H,
                'initial_driver_sum_s_inv':gamma0*(c.m0-meq)+J/H,
                'top_liquid_supply_required_initially':bool(J<0),
                'Gamma_window_prescribed':chi*IU,'Gamma_at_first_hit':rec['Gamma_at_first_hit'],
                'Gamma_c_asy':rec['Gamma_c_asy'],'chi_asy_m_inv':rec['chi_asy_m_inv'],
                'reached_status':reached_state(state),'origin':'computed_v3; shared hydrostatic cases are regression reruns',**state})
            if q0==.5:
                oldregime='design' if regime=='common_design' else 'critical_upper_neighbor'
                original=oldcum[(oldcum.case_id=='V2_01')&(oldcum.regime==oldregime)].iloc[0]
                onset_error=abs(state['onset_elapsed_s']-float(original.diagnostic_cut_elapsed_s))
                gamma_error=abs(rec['Gamma_at_first_hit']-float(original.Gamma_at_first_hit))
                regressions.append({'case_id':cid,'regime':regime,'source_case':'V2_01',
                    'source_regime':oldregime,'onset_difference_s':onset_error,
                    'Gamma_hit_difference':gamma_error,'pass':onset_error<1e-6 and gamma_error<1e-7})
            print(cid,q0,regime,reached_state(state),'hit',state['onset_elapsed_s'],flush=True)
        if q0==.5:
            root=float(reference_old.chi_star_m_inv)
            oldroot=oldbr[(oldbr.case_id=='V2_01')&(oldbr.record_type=='valid_root_bracket')].iloc[0]
            grid=oldbr[(oldbr.case_id=='V2_01')&(oldbr.record_type=='search_grid')].drop(columns=['case_id','origin']).to_dict('records')
            eps=max(cfg['threshold_search']['xtol_m_inv']*10,root*1e-5)
            neighbor={}
            for side,chi in [('left',root-eps),('right',root+eps)]:
                _,state=run_case(replace(case,chi_m_inv=chi),forcing,num,False)
                neighbor[f'neighbor_{side}_chi_m_inv']=chi
                neighbor[f'neighbor_{side}_margin']=state['min_mean_effective_ratio']-case.qcrit if state['status']=='completed_conditional_window' else None
                neighbor[f'neighbor_{side}_status']=state['status']
            if not (neighbor['neighbor_left_margin']>0 and neighbor['neighbor_right_margin']<0):
                raise RuntimeError('Reused hydrostatic root lost its valid opposite-sign neighbors')
            grid.append({'record_type':'valid_root_bracket','chi_m_inv':root,
                'bracket_lo_m_inv':float(oldroot.bracket_lo_m_inv),'bracket_hi_m_inv':float(oldroot.bracket_hi_m_inv),
                'bracket_margin_lo':float(oldroot.bracket_margin_lo),'bracket_margin_hi':float(oldroot.bracket_margin_hi),
                'neighbor_check_origin':'v2 exact root reused; two neighbors actually recomputed v3',**neighbor})
            tr={'chi_critical_m_inv':root,'threshold_status':'root_found_in_valid_bracket',
                'critical_identity':'v2 hydrostatic root reused; both actual neighbors rechecked v3',
                'valid_complete_window_grid_points':int(reference_old.valid_complete_window_grid_points),
                'incomplete_or_failed_grid_points':int(reference_old.incomplete_or_failed_grid_points)}
            origin='v2_exact_root_reused_with_actual_v3_neighbor_regression'
        else:
            tr,grid=threshold(case,forcing,cfg,num);origin='computed_v3_new_initial_state_boundary'
        for row in grid:brackets.append({'case_id':cid,'q0':q0,'origin':origin,**row})
        root=tr['chi_critical_m_inv'];valid=[g for g in grid if g.get('record_type')=='valid_root_bracket']
        Delta=case.mcrit-case.m0
        Gasy=case.a/(case.m0*Delta)*math.log((case.s_init+case.sigma0)/(case.qcrit*case.E+case.sigma0))
        bound={'case_id':cid,'q0':q0,'s_init_Pa':case.s_init,'qcrit':case.qcrit,'qmax':admission['qmax'],
            'm0':case.m0,'mcrit':case.mcrit,'H0_m':case.H0,'shape_n':case.shape_n,
            'initial_state_origin':en['initial_state_identity'],'criterion_origin':'closure_derived',
            'forcing_origin':'post_event_inversion','forcing_window':'[27,100]','T_s':73,'I_U_m':IU,
            'Gamma_c_asy':Gasy,'chi_asy_m_inv':Gasy/IU,'chi_star_m_inv':root,
            'Gamma_star_window':None if root is None else root*IU,
            'chi_star_relative_error_vs_asy':None if root is None else (root-Gasy/IU)/(Gasy/IU),
            'chi_star_bracket':json.dumps([valid[0]['bracket_lo_m_inv'],valid[0]['bracket_hi_m_inv']]) if valid else None,
            'origin':origin,**tr}
        if root is not None:
            eps=max(cfg['threshold_search']['xtol_m_inv']*10,root*1e-5)
            c=replace(case,chi_m_inv=root+eps)
            f,state,sol=run_case(c,forcing,num,keep_solution=True)
            f.insert(0,'case_id',cid);f['regime']='own_critical_upper_neighbor';f['q0']=q0;trajectories.append(f)
            rec,ser=shear_diagnostics(c,forcing,sol,state,cid,'own_critical_upper_neighbor',en['diagnostic_max_step_s'])
            rec['chi_star_m_inv']=root;rec['chi_eval_minus_star_m_inv']=eps
            diagnostics.append(rec);series.append(ser)
            bound.update({'chi_neighbor_eval_m_inv':root+eps,'neighbor_status':reached_state(state),
                'neighbor_completion_status':state['status'],'neighbor_hit_elapsed_s':state['onset_elapsed_s'],
                'Gamma_neighbor_at_first_hit':rec['Gamma_at_first_hit'],
                'eps_m_at_own_neighbor':rec['eps_m'],'eps_N_at_own_neighbor':rec['eps_N'],
                'drain_to_shear_at_own_neighbor':rec['absolute_drain_to_shear_ratio'],
                'L_goal':rec['L_goal'],'eps_J_at_own_neighbor':rec['eps_J']})
        bounds.append(bound)

        c=replace(case,mode='zero_shear',chi_m_inv=0)
        f,state,sol=run_case(c,forcing,num,keep_solution=True)
        f.insert(0,'case_id',cid);f['regime']='zero_shear_relaxation';f['q0']=q0;trajectories.append(f)
        rec,ser=shear_diagnostics(c,forcing,sol,state,cid,'zero_shear_relaxation',en['diagnostic_max_step_s'])
        rec['approximation_status']='not_applicable_no_shear_relaxation_control'
        diagnostics.append(rec);series.append(ser)
        dm,ds,J=rhs(c,forcing,0,(c.m0,c.s_init,0))
        consistent=(abs(state['final_m']-case.m0)<1e-12 and abs(f.mean_effective_stress_Pa.iloc[-1]-case.s_init)<1e-8) if q0==.5 else bool((state['final_m']-case.m0)*dm>0 and (f.mean_effective_stress_Pa.iloc[-1]-case.s_init)*ds>0)
        relax.append({'case_id':cid,'q0':q0,'initial_m_dot_s_inv':dm,'initial_s_dot_Pa_s':ds,'J0_m_s':J,
            'delta_m':state['final_m']-case.m0,'delta_s_Pa':float(f.mean_effective_stress_Pa.iloc[-1]-case.s_init),
            'direction_or_static_check_pass':consistent,
            'identity':'numerical zero-shear pressure-relaxation / budget check; not external validation',**state})
        print(cid,q0,'own boundary',tr['threshold_status'],root,flush=True)
    if not all(r['pass'] for r in regressions):raise RuntimeError('Hydrostatic v2 regression failed')
    if not all(r['direction_or_static_check_pass'] for r in relax):raise RuntimeError('Zero-shear pressure relaxation check failed')
    qmax=initial_classification(ref)['qmax']
    admissions=[initial_classification(replace(ref,initial_stress_ratio=q)) for q in [0,.1,ref.qcrit,.17,.30,.5,.75,qmax,qmax+1e-4,-.01]]
    expected=['initial_closure_limit','initial_positive_thickness_layer','initial_first_layer_boundary',
        *(['admissible_unreached_initial_state']*4),'initial_saturation_boundary',
        'invalid_initial_state','invalid_initial_state']
    if [row['initial_status'] for row in admissions]!=expected:raise RuntimeError('Initial admission classification check failed')
    for row in admissions:
        row['root_search_in_this_classification_record']=False
        row['separate_dynamic_boundary_table']='initial_state_boundary.csv' if row['q0'] in en['q0_values'] else None
        row['Gamma_initial_condition']=0.0 if row['initial_condition_satisfied'] else None
    tables={'initial_state_response.csv':pd.DataFrame(summaries),'initial_state_boundary.csv':pd.DataFrame(bounds),
        'initial_state_search_brackets.csv':pd.DataFrame(brackets),'initial_state_admission.csv':pd.DataFrame(admissions),
        'initial_state_trajectories.csv':pd.concat(trajectories,ignore_index=True),
        'initial_state_cumulative_check.csv':pd.DataFrame(diagnostics),
        'initial_state_cumulative_series.csv':pd.concat(series,ignore_index=True),
        'initial_state_relaxation_checks.csv':pd.DataFrame(relax),'v2_hydrostatic_regression.csv':pd.DataFrame(regressions)}
    for name,frame in tables.items():frame.to_csv(out/name,index=False)
    save_json(out/'summary.json',{'version':cfg['version'],'main_response_settings':8,'initial_state_boundaries':4,
        'new_root_searches':3,'v2_exact_roots_reused':1,'numerical_zero_shear_controls':4,
        'admission_classifications':len(admissions),'diagnostic_trajectories':len(diagnostics),
        'diagnostic_saved_rows':sum(len(x) for x in series),'qmax':qmax,
        'historical_result_files_reused':reused,'actual_event_unknowns':cfg['actual_event_unknowns'],
        'scientific_status':'conditional initial-state sensitivity; local initial state remains unobserved'})
    inputs={str(p.relative_to(delivery)):fingerprint(p) for p in sorted((delivery/'data').rglob('*')) if p.is_file()}
    for name in ('config.json',):inputs[name]=fingerprint(delivery/name)
    save_json(out/'ANALYSIS_RUN_RECEIPT.json',{'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat(),
        'command':sys.argv,'code_sha256':fingerprint(__file__),'input_sha256':inputs,
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
        'result_sha256':{p.name:fingerprint(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='ANALYSIS_RUN_RECEIPT.json'},
        'execution':'completed','checkpoint_B':'pending','independent_review':'pending','local_validation':'unresolved'})
    print('Completed frozen v3 initial-state matrix, joint boundaries and saved-state diagnostics',flush=True)


if __name__=='__main__':sol()

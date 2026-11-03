import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STYLE={'font':'Arial Unicode MS','blue':'#246A8A','orange':'#C57934','purple':'#77608D',
       'grey':'#777E86','green':'#467F69','background':'#FFFFFF','dpi':300,
       'width_in':7.2,'height_in':3.65,'font_size':9}


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pair():
    fig,axes=plt.subplots(1,2,figsize=(STYLE['width_in'],STYLE['height_in']))
    fig.subplots_adjust(left=.10,right=.98,bottom=.23,top=.84,wspace=.40)
    for i,ax in enumerate(axes):
        ax.spines[['top','right']].set_visible(False)
        ax.tick_params(labelsize=8)
        ax.set_gid(chr(97+i))
    return fig,axes


def export(fig,axes,out,name,title,note):
    fig.suptitle(title,y=.975,fontsize=11)
    fig.text(.10,.055,note,ha='left',va='bottom',fontsize=8,color='#47505A')
    fig.canvas.draw()


    fig.savefig(out/(name+'.png'),facecolor='white',dpi=300)
    fig.savefig(out/(name+'.svg'),facecolor='white')
    fig.savefig(out/(name+'.pdf'),facecolor='white')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--delivery-dir',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--results-dir',type=Path)
    parser.add_argument('--output-dir',type=Path)
    a=parser.parse_args();d=a.delivery_dir.resolve();r=(a.results_dir or d/'results').resolve();out=(a.output_dir or d/'figures').resolve()
    out.mkdir(parents=True,exist_ok=True)
    started=datetime.now(timezone.utc).isoformat()
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':[STYLE['font'],'DejaVu Sans'],
                         'font.size':STYLE['font_size'],'axes.labelsize':9,'axes.titlesize':10,
                         'legend.fontsize':8,'axes.unicode_minus':False,'svg.fonttype':'none',
                         'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':STYLE['dpi']})
    b=STYLE['blue'];o=STYLE['orange'];g=STYLE['grey'];p=STYLE['purple']
    curve=pd.read_csv(d/'data/cghc_f7b_digitized.csv');gsd=pd.read_csv(d/'data/blatten_gsd_long.csv')
    assert np.all(gsd.diameter_mm>0), 'Log grain-size axis requires strictly positive diameters'
    f,ax=pair()
    t=curve.time_since_refined_onset_s;u=curve.velocity_m_s
    ax[0].fill_between(t,np.maximum(u-1.5,0),u+1.5,color=b,alpha=.16,lw=0)
    ax[0].plot(t,u,'o-',color=b,ms=3,lw=1.3)
    ax[0].axvspan(27,100,color=o,alpha=.10,zorder=0)
    ax[0].set(xlabel='相对源释放时间 / s',ylabel='作者反演速度 / (m/s)',xlim=(0,100),ylim=(0,78))
    ax[0].set_title('a  同事件运动曲线',loc='left')
    for color,(name,group) in zip((b,o,p,STYLE['green']),gsd.groupby('sample',sort=False)):
        group=group.sort_values('diameter_mm')
        ax[1].semilogx(group.diameter_mm,group.coarse_adjusted_passing_percent,'o-',color=color,ms=2.5,lw=1.1,label=name.replace('V-',''))
    ax[1].axhline(10,color=g,ls=':',lw=.8)
    ax[1].set(xlabel='粒径 / mm',ylabel='粗料校正累计通过率 / %',ylim=(0,105),xlim=(.001,100))
    ax[1].set_title('b  可比材料的四份真实样本',loc='left');ax[1].legend(loc='upper left',frameon=False)
    export(f,ax,out,'fig01_new_evidence','新增证据：外部驱动与具名材料类比',
           '施剪窗27–100 s；速度带仅为读图 ±1.5 m/s；粒度为2025 Blatten沉积样本。')

    tr=pd.read_csv(r/'conditional_trajectories.csv')
    f,ax=pair()
    for name,label,color,ls in [('reference','有限排水参考',b,'-'),('ideal_drained','理想自由排水',g,'--')]:
        z=tr[(tr.case_id==name)&(tr.elapsed_s<=10)]
        ax[0].plot(z.elapsed_s,z.mean_effective_ratio,color=color,ls=ls,label=label)
    ax[0].axhline(1/6,color=o,ls=':',lw=1.2,label='首层应力边界 1/6')
    ax[0].set(xlabel='外部施剪开始后 / s',ylabel='平均有效应力比 s/E',xlim=(0,10),ylim=(-.01,.55))
    ax[0].set_title('a  排水对照分离孔压反馈',loc='left');ax[0].legend(frameon=False,loc='center right')
    for name,label,color,ls in [('reference','m0=0.54',b,'-'),('finite_H1_m0.57','m0=0.57',o,'-'),('finite_H1_m0.63','m0=0.63（出域前）',g,'--')]:
        z=tr[(tr.case_id==name)&(tr.elapsed_s<=10)]
        ax[1].plot(z.elapsed_s,z.fluid_supported_fraction,color=color,ls=ls,label=label)
        if name.endswith('0.63'):
            ax[1].plot(z.elapsed_s.iloc[-1],z.fluid_supported_fraction.iloc[-1],'x',color=g,ms=6)
    ax[1].set(xlabel='外部施剪开始后 / s',ylabel='条件孔压支承厚度比例 f',xlim=(0,10),ylim=(-.03,1.03))
    ax[1].set_title('b  初始结构改变条件响应',loc='left');ax[1].legend(frameon=False,loc='center right')
    export(f,ax,out,'fig02_response_control','条件状态响应：相同驱动下的排水与结构差异',
           'H0=1 m、χ=0.02 /m为设计参考；“出域”是饱和假说退出，不能判为安全。')

    crit=pd.read_csv(r/'conditional_critical_chi.csv');grid=pd.read_csv(r/'critical_search_brackets.csv')
    assert np.all(crit.H0_m>0) and np.all(grid.chi_m_inv>0), 'Log axes require strictly positive input values'
    f,ax=pair()
    for m,color in ((.54,b),(.57,o)):
        z=crit[(crit.case_id.str.startswith('finite_'))&np.isclose(crit.m0,m)].sort_values('H0_m')
        ax[0].plot(z.H0_m,z.chi_critical_m_inv*1e4,'o-',color=color,label=f'm0={m:.2f}')
    ax[0].set_xscale('log');ax[0].set(xlabel='设计初始厚度 H0 / m',ylabel='有效括区临界 χ / (10⁻⁴ m⁻¹)',xlim=(.07,14))
    ax[0].set_title('a  12个有限成员中的有效反求',loc='left');ax[0].legend(frameon=False,loc='upper left')
    z=grid[(grid.case_id=='finite_H1_m0.54')&(grid.record_type=='search_grid')].sort_values('chi_m_inv')
    ax[1].semilogx(z.chi_m_inv,z.margin_mean_stress_ratio,'o-',color=b,ms=3,lw=1)
    ax[1].axhline(0,color=o,lw=1,ls=':')
    root=float(crit.loc[crit.case_id=='finite_H1_m0.54','chi_critical_m_inv'].iloc[0])
    ax[1].plot(root,0,'D',color=o,ms=5,label='已核有效括区的根')
    ax[1].set(xlabel='施剪转移组合 χ / m⁻¹',ylabel='窗内最小 s/E − 1/6')
    ax[1].set_title('b  有限搜索保留非单调响应',loc='left');ax[1].legend(frameon=False,loc='upper right')
    export(f,ax,out,'fig03_critical_conditions','反求临界组合：有限括区成立，全域单调性未证',
           '较密成员含提前出域；χ及1/χ为组合条件，均非实测水深或剪切厚度。')

    raw=pd.read_csv(d/'data/iverson2000_sparse_process_keypoints.csv')
    obs=pd.concat([pd.read_csv(d/'data/iverson2000_matched_depth67_extra_keypoints.csv'),
                   raw[raw.experiment_id.str.contains('dense')]],ignore_index=True)
    f,ax=pair()
    for axis,(exp,z),color,label in zip(ax,obs.groupby('experiment_id',sort=False),(b,o),('a  松土：67 cm深传感器','b  密土：67 cm深传感器')):
        z=z.sort_values('time_since_water_application_s')
        time=z.time_since_water_application_s-z.time_since_water_application_s.iloc[0]
        axis.errorbar(time,z.pore_water_pressure_head_cm,xerr=z.time_reading_error_s,yerr=z.head_reading_error_cm,
                      fmt='o-',color=color,ms=4,capsize=2,lw=1.2)
        axis.set(xlabel='所读过程首点后 / s',ylabel='实测孔压水头 / cm',ylim=(0,72))
        axis.set_title(label,loc='left')
    export(f,ax,out,'fig04_external_process','外部过程约束：松土增压与密土降压',
           '同一研究的两次代表过程、共8个稀疏读点；误差为读图支持，非独立样本置信区间。')
    closure=pd.read_csv(r/'closure_sensitivity.csv');profiles=pd.read_csv(r/'closure_profile_samples.csv')
    f,ax=pair()
    for n,color in ((2,b),(3,o)):
        z=profiles[np.isclose(profiles.shape_n,n)]
        ax[0].plot(z.effective_stress_over_E,z.height_fraction,color=color,label=f'n={n}')
    ax[0].set(xlabel='有效法向应力 / E',ylabel='相对高度 z/H',xlim=(0,.72),ylim=(0,1))
    ax[0].set_title('a  边界剖面随闭合形状改变',loc='left');ax[0].legend(frameon=False,loc='upper right')
    for n,color in ((2,b),(3,o)):
        z=closure[np.isclose(closure.shape_n,n)].iloc[0]
        ax[1].plot(n,z.chi_star_m_inv*1e4,'o',color=color,ms=7)
        ax[1].text(n,z.chi_star_m_inv*1e4+.09,f'sc/E={z.qcrit:.3f}',ha='center',fontsize=8)
    ax[1].set(xlabel='数学压力形状 n',ylabel='参考成员临界 χ* / (10⁻⁴ m⁻¹)',xlim=(1.5,3.5),ylim=(1.1,2.6),xticks=[2,3])
    ax[1].set_title('b  同一参考成员的条件反求',loc='left')
    export(f,ax,out,'fig05_closure_sensitivity','压力形状改变模型内边界，n=2仍为主模型',
           'n=3仅作数学压力测试；阈值比、应力映射、排水及停止边界成对改变。')

    cum=pd.read_csv(r/'cumulative_shear_series.csv');mw=pd.read_csv(r/'material_window_sensitivity.csv')
    f,ax=pair()
    reference_id=str(closure.loc[np.isclose(closure.shape_n,2),'case_id'].iloc[0])
    for regime,label,color,ls in [('design','设计施剪 χ=0.02',b,'-'),('critical_upper_neighbor','临界根上邻侧',o,'-')]:
        z=cum[(cum.case_id==reference_id)&(cum.regime==regime)]
        ax[0].plot(z.Gamma_sh,z.s_ratio,color=color,ls=ls,lw=1.4,label=label)
    z=cum[(cum.case_id==reference_id)&(cum.regime=='design')]
    ax[0].plot(z.Gamma_sh,z.s_asy_ratio,color=g,ls='--',lw=1.2,label='弱排水、小N近似')
    ax[0].axhline(1/6,color=p,ls=':',lw=.8)
    ax[0].set(xlabel='模型输入累计剪切 Γsh',ylabel='首层前平均应力比 s/E',xlim=(0,1.04),ylim=(.145,.52))
    ax[0].set_title('a  相同累计剪切下的差别',loc='left');ax[0].legend(frameon=False,loc='upper right')
    z=mw[mw.block=='window'].sort_values('T_s')
    ax[1].plot(z.T_s,z.chi_star_m_inv*1e4,'o-',color=b,label='完整模型有效括区根')
    ax[1].plot(z.T_s,z.chi_asy_m_inv*1e4,'s--',color=o,ms=4,label='累计剪切解析近似')
    ax[1].set(xlabel='同起点前缀窗长 T / s',ylabel='临界 χ* / (10⁻⁴ m⁻¹)',xticks=[33,53,73])
    ax[1].set_title('b  前缀窗口的临界组合',loc='left');ax[1].legend(frameon=False,loc='upper right')
    export(f,ax,out,'fig06_cumulative_shear','累计剪切解释条件门槛，同时保留近似误差',
           'Γsh无量纲；设计响应与临界根分开。窗口同时改变累计驱动和排水时间。')

    boundary=pd.read_csv(r/'initial_state_boundary.csv').sort_values('q0')
    response=pd.read_csv(r/'initial_state_response.csv')
    diagnostics=pd.read_csv(r/'initial_state_cumulative_check.csv')
    initial_series=pd.read_csv(r/'initial_state_cumulative_series.csv')
    assert len(boundary)==4 and len(response)==8, 'Frozen four-state / two-rate matrix required'
    f,ax=pair()
    ax[0].plot(boundary.q0,boundary.Gamma_c_asy,'--',color=g,lw=1,label='一般初值解析近似')
    ax[0].plot(boundary.q0,boundary.Gamma_star_window,'o',color=b,ms=6,label='完整窗有效根（有限点）')
    comparison=float(response.loc[response.regime=='common_comparison','Gamma_window_prescribed'].iloc[0])
    ax[0].axhline(comparison,color=o,ls=':',lw=1,label='共同比较档的窗末 Γ')
    ax[0].set(xlabel='设计初始有效应力比 q0=s_init/E',ylabel='自身临界组合 χ*I_U（无量纲）',
              xlim=(.13,.79),ylim=(-.04,1.62),xticks=boundary.q0)
    ax[0].set_title('a  初态—累计剪切联合边界',loc='left')
    ax[0].legend(frameon=False,loc='upper left',fontsize=7)
    for rowindex,(regime,color) in enumerate([('common_comparison',o),('common_design',b)]):
        z=response[response.regime==regime].sort_values('q0')
        for j,row in enumerate(z.itertuples()):
            hit=pd.notna(row.onset_elapsed_s)
            ax[1].plot(j,rowindex,marker='o' if hit else 's',ms=6,
                       markerfacecolor=color if hit else 'white',markeredgecolor=color,lw=0)
            label=f'{row.onset_elapsed_s:.3g} s' if hit else '窗内未达'
            if hit and regime=='common_comparison' and row.q0==.5:
                label=f'{row.onset_elapsed_s:.3f} s'
            ax[1].text(j,rowindex+.20,label,ha='center',va='bottom',fontsize=7.5,color=color)
    ax[1].set(xlabel='设计初始应力比 q0',xlim=(-.5,3.5),ylim=(-.45,1.65),
              xticks=range(4),xticklabels=[f'{x:g}' for x in boundary.q0],
              yticks=[0,1],yticklabels=['共同比较档','设计档'])
    ax[1].set_title('b  两个共同 χ 下的实际状态',loc='left')
    export(f,ax,out,'fig07_initial_state_boundary','约0.98的参考累计剪切条件依赖初始状态',
           '四点均为设计初态；数值根只显示实际有限点。空方块为73 s内未达，不赋到达时间。')

    f,ax=pair()
    near=str(boundary.iloc[0].case_id)
    for regime,label,color in [('own_critical_upper_neighbor','自身临界上邻侧',b),
                                ('common_comparison','共同比较档',o)]:
        z=initial_series[(initial_series.case_id==near)&(initial_series.regime==regime)]
        rec=diagnostics[(diagnostics.case_id==near)&(diagnostics.regime==regime)].iloc[0]
        assert rec.Gamma_c_asy>0 and rec.L_goal>0, 'Positive pre-layer targets required'
        progress=-z.log_stress_ratio/rec.L_goal
        ax[0].plot(z.Gamma_sh/rec.Gamma_c_asy,progress,color=color,lw=1.4,label=label)
    ax[0].plot([0,1],[0,1],'--',color=g,lw=1,label='弱排水解析近似')
    ax[0].set(xlabel='累计剪切 / 一般初值解析 Γc',ylabel='对数应力降 / 初始目标 L_goal',
              xlim=(0,1.13),ylim=(0,1.17))
    ax[0].set_title('a  q0=0.17：小目标放大排水作用',loc='left')
    ax[0].legend(frameon=False,loc='upper left',fontsize=7)
    for regime,label,color,marker in [('own_critical_upper_neighbor','自身临界上邻侧',b,'o'),
                                      ('common_comparison','共同比较档',o,'s'),
                                      ('common_design','设计档',p,'^')]:
        z=diagnostics[diagnostics.regime==regime].sort_values('q0')
        assert np.all(z.eps_J>0), 'Log drainage diagnostic requires positive values'
        ax[1].semilogy(z.q0,z.eps_J,marker+'-',color=color,lw=1,ms=4,label=label)
    ax[1].set(xlabel='设计初始应力比 q0',ylabel='绝对排水对数贡献 / L_goal（εJ）',
              xlim=(.13,.79),xticks=boundary.q0,ylim=(5e-6,.3))
    ax[1].set_title('b  排水须相对于目标间隔衡量',loc='left')
    ax[1].legend(frameon=False,loc='upper right',fontsize=7)
    export(f,ax,out,'fig08_initial_state_diagnostics','近边界的弱施剪可使排水改变累计剪切近似',
           '诊断由保存轨迹直接积分；不设置通过百分比。数值残差检验约化关系，不证明自然精度。')

    sources=[d/'data/cghc_f7b_digitized.csv',d/'data/blatten_gsd_long.csv',d/'data/iverson2000_sparse_process_keypoints.csv',
             d/'data/iverson2000_matched_depth67_extra_keypoints.csv',
             r/'conditional_trajectories.csv',r/'conditional_critical_chi.csv',r/'critical_search_brackets.csv',
             r/'closure_sensitivity.csv',r/'closure_profile_samples.csv',r/'cumulative_shear_series.csv',r/'material_window_sensitivity.csv',
             r/'initial_state_boundary.csv',r/'initial_state_response.csv',r/'initial_state_cumulative_check.csv',r/'initial_state_cumulative_series.csv']
    receipt={'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
             'input_sha256':{str(x):sha(x) for x in sources},'style':STYLE,'solver_started':False,
             'output_sha256':{x.name:sha(x) for x in sorted(out.iterdir()) if x.is_file() and x.name!='FIGURES_RUN_RECEIPT.json'},
             'visual_inspection':'pending'}
    (out/'FIGURES_RUN_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    print('Rendered 8 figures as PNG/SVG/PDF')


if __name__=='__main__': main()

from __future__ import annotations
import argparse
import gzip
import json
import math
import os
from pathlib import Path
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'q1_final_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

STYLE={'width_mm':156.,'dpi':400,'font':9.,'small':8.,'panel':9.5,'line':1.5,
       'target':'#164D68','accent':'#B66B32','neutral':'#62676D','band':'#DEE1E3',
       'grid':'#E6E8EA','missing':'#E4E6E8','colors':['#164D68','#B66B32','#718B9A','#6E637C','#8C8C66'],
       'display_seed':20261026,
       'fonts':['Noto Sans CJK SC','Source Han Sans SC','Microsoft YaHei','PingFang SC','Heiti SC','SimHei','Arial Unicode MS']}
SHORT={'rh2_mean':'2 m相对湿度','t20_t2_diff':'20−2 m气温差','ws10_mean':'10 m标量风速',
       'warm_rate':'24 h块均温变化率','ws20_ws2_diff':'20−2 m标量风差','pdd':'正温指数PDD',
       'rsd_mean':'Rsd原通道均值','vwc01':'0.1 m体积含水量','soilt01':'0.1 m土温',
       't2_max':'小时均温最大值代理','rld_mean':'Rld原通道均值','g01':'0.1 m板热流','p_mean':'MET气压均值'}
WINDOWS=(24,72,168,360,720)


def configure_style():
    fonts={f.name for f in font_manager.fontManager.ttflist}
    chinese=next((s for s in STYLE['fonts'] if s in fonts),None)
    if chinese is None:raise RuntimeError('需要可用中文字体，如Noto Sans CJK SC。')
    matplotlib.rcParams.update({'font.family':['DejaVu Sans',chinese],'font.sans-serif':[chinese,'DejaVu Sans'],
        'font.size':STYLE['font'],'axes.labelsize':STYLE['font'],'axes.titlesize':STYLE['panel'],
        'xtick.labelsize':STYLE['font'],'ytick.labelsize':STYLE['font'],'legend.fontsize':8.5,
        'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,
        'axes.unicode_minus':True,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
        'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
        'text.color':'#252A2E','axes.labelcolor':'#252A2E','xtick.color':'#44494D','ytick.color':'#44494D',
        'legend.frameon':False,'lines.linewidth':STYLE['line']})
    return chinese


def load_json(path):
    p=Path(path)
    with gzip.open(p,'rt',encoding='utf-8') if p.suffix=='.gz' else p.open(encoding='utf-8') as stream:
        return json.load(stream)


def array(values):
    return np.asarray([float(v) if v is not None and math.isfinite(float(v)) else np.nan for v in values],dtype=float)


def grid(rows,cols,height,left=16,right=4,bottom=16,top=14,hgap=18,vgap=21):

    w=STYLE['width_mm'];fig=plt.figure(figsize=(w/25.4,height/25.4))
    pw=(w-left-right-(cols-1)*hgap)/cols;ph=(height-top-bottom-(rows-1)*vgap)/rows
    axes=[]
    for row in range(rows):
        axes.append([fig.add_axes([(left+col*(pw+hgap))/w,(height-top-(row+1)*ph-row*vgap)/height,pw/w,ph/height]) for col in range(cols)])
    return fig,np.asarray(axes)


def panel(ax,letter,title):
    ax.text(0,1.20,f'({letter}) {title}',transform=ax.transAxes,ha='left',va='bottom',fontsize=STYLE['panel'])
    ax.tick_params(direction='out',length=3,width=.7,pad=3)
    ax.grid(axis='y',color=STYLE['grid'],linewidth=.5,zorder=0)
    ax.set_axisbelow(True)


def note(fig,text,y=.012):
    fig.text(.025,y,text,fontsize=STYLE['small'],color=STYLE['neutral'],va='bottom')


def time_axes(ax,start,label=True):
    ax.set_xlim(start,0)
    ax.set_xticks([-30,-20,-10,0] if start==-30 else [-7,-5,-3,-1,0] if start==-7 else [-3,-2,-1,0])
    if label:ax.set_xlabel('距截止时刻（日）',labelpad=3)


def support(ch,period):
    s=ch['summaries'][period]
    return f"{s['available_hours']}/{s['expected_hours']} h"


def background(ax,x,ch,start):
    mask=(x>start)&(x<=0);xx=x[mask]
    lo=array(ch['background_q25'])[mask];mid=array(ch['background_median'])[mask];hi=array(ch['background_q75'])[mask]
    finite=np.isfinite(lo)&np.isfinite(hi)
    if np.any(lo[finite]>hi[finite]):raise ValueError('背景分位次序错误。')
    ax.fill_between(xx,lo,hi,color=STYLE['band'],linewidth=0,zorder=1)
    ax.plot(xx,mid,color=STYLE['neutral'],lw=.9,ls=(0,(3,2)),zorder=2)


def target(ax,x,ch,start,color=None,label=None,linestyle='-',bg=True,lw=None):
    if bg:background(ax,x,ch,start)
    mask=(x>start)&(x<=0)
    ax.plot(x[mask],array(ch['values'])[mask],color=color or STYLE['target'],lw=lw or STYLE['line'],
            ls=linestyle,label=label,zorder=3)


def context_legend(fig,y=.99,ncol=3):
    handles=[Line2D([],[],color=STYLE['target'],lw=1.5,label='目标过程'),
             Line2D([],[],color=STYLE['neutral'],ls=(0,(3,2)),lw=.9,label='历史中位数'),
             Patch(facecolor=STYLE['band'],edgecolor='none',label='历史25%—75%')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,y),ncol=ncol,
               handlelength=1.8,columnspacing=1.0,handletextpad=.4,fontsize=8.5)


def channel_map(process):
    return {c['id']:c for c in process['weather']['channels']}


def validate(data,process,pca):
    branch=data['primary_branch'];meta=process['metadata']
    if branch['endpoint']!='H-END' or branch['cutoff_hour']!=10 or meta['primary_cutoff']!='2026-08-26T10:00:00+08:00':
        raise ValueError('输入不是冻结的UTC+8主截止/端点分支。')
    if pca['metadata']['endpoint']!='H-END' or pca['metadata']['cutoff']!=meta['primary_cutoff']:
        raise ValueError('PCA与过程时间条件不一致。')
    ids=[r['id'] for r in data['primary']['candidate_ranking']]
    if ids!=pca['feature_ids'] or len(ids)!=13:raise ValueError('候选池/PCA输入不一致。')
    return ids


def figure_structure(data,process,pca):
    fig,axes=grid(2,2,130,top=22,bottom=20,vgap=24)
    near=process['near_surface'];heights=array(near['heights_m'])
    profiles={int(p['support_end_relative_hours']):p for p in near['profiles']}
    for col,(key,unit,title) in enumerate((('temperature','°C','气温剖面'),('wind','m/s','标量风速剖面'))):
        ax=axes[0,col]
        ref=profiles[0][key]
        ax.fill_betweenx(heights,array(ref['background_q25']),array(ref['background_q75']),color=STYLE['band'],lw=0)
        ax.plot(array(ref['background_median']),heights,ls=(0,(3,2)),color=STYLE['neutral'],lw=.9)
        for hour,color,marker in ((-24,STYLE['accent'],'s'),(0,STYLE['target'],'o')):
            ax.plot(array(profiles[hour][key]['values']),heights,color=color,marker=marker,lw=1.5,ms=3.5)
        ax.set_ylim(0,21);ax.set_yticks([2,10,20]);ax.set_ylabel('高度 / m');ax.set_xlabel(unit)
        panel(ax,chr(97+col),title)
    th=near['time_height'];x=array(th['support_end_relative_hours'])/24
    for col,(value,prefix,skey,unit,title) in enumerate((
        ('temperature_difference_c','pair_temperature','temperature_difference_summaries','°C','T20−T2：末7日'),
        ('wind_difference_ms','pair_wind','wind_difference_summaries','m/s','U20−U2：末7日'))):
        ch={'values':th[value],**{f'background_{q}':th[f'{prefix}_background_{q}'] for q in ('q25','median','q75')}}
        ax=axes[1,col];target(ax,x,ch,-7);ax.axhline(0,color=STYLE['neutral'],lw=.65)
        ax.set_ylabel(unit);time_axes(ax,-7);panel(ax,chr(99+col),title)
        s=th[skey]['7d'];ax.text(0,1.015,f"n={s['available_hours']}/{s['expected_hours']} h",transform=ax.transAxes,va='bottom',fontsize=8)
    fig.legend(handles=[Line2D([],[],color=STYLE['accent'],marker='s',ms=3.5,label='−24 h剖面'),
        Line2D([],[],color=STYLE['target'],marker='o',ms=3.5,label='截止剖面'),
        Line2D([],[],color=STYLE['neutral'],ls=(0,(3,2)),lw=.9,label='历史中位数'),
        Patch(facecolor=STYLE['band'],edgecolor='none',label='历史四分位')],loc='upper center',bbox_to_anchor=(.53,.995),ncol=4,
        fontsize=8,handlelength=1.3,columnspacing=.8,handletextpad=.35)
    note(fig,'北京时间8月26日10:00截止；剖面背景取截止同钟点，层差背景逐时匹配。')
    return fig


def add_extremes(ax,process,x,ch_id,start):
    for e in process['extremes']['series']:
        if e['channel_id']!=ch_id:continue
        hit=np.asarray([v is True for v in e['exceedance']])&(x>start)
        color=STYLE['accent'] if e['direction']=='high' else STYLE['neutral']
        ax.plot(x[hit],array(e['values'])[hit],ls='none',marker='^' if e['direction']=='high' else 'v',ms=3.0,color=color,zorder=5)


def figure_process(data,process,pca):
    fig,axes=grid(3,2,171,top=22,bottom=20,vgap=23)
    ch=channel_map(process);x=array(process['weather']['time_relative_hours'])/24
    plan=[('t2',-30,'气温：30日','°C'),('t2',-3,'气温：末72 h','°C'),
          ('rh2',-30,'相对湿度：30日','%'),('vwc01',-30,'0.1 m含水量：30日','m³/m³'),
          ('radiation',-30,'辐射原通道：30日','W/m²'),('ws10',-3,'标量风速：末72 h','m/s')]
    for i,(key,start,title,unit) in enumerate(plan):
        ax=axes.flat[i]
        if key=='radiation':
            for g,color in (('rsd',STYLE['target']),('rld',STYLE['accent'])):
                target(ax,x,ch[g],start,color=color,label=g.upper(),linestyle='-' if g=='rsd' else '--',bg=False)
            ns=[support(ch[g],'30d') for g in ('rsd','rld')]
            ntext='Rsd实 / Rld虚；各'+ns[0] if ns[0]==ns[1] else 'Rsd实 '+ns[0]+'；Rld虚 '+ns[1]
        else:
            target(ax,x,ch[key],start)
            if key in ('t2','ws10'):add_extremes(ax,process,x,key,start)
            ntext='n='+support(ch[key],'30d' if start==-30 else '72h')
        ax.set_ylabel(unit);time_axes(ax,start);panel(ax,chr(97+i),title)
        ax.text(0,1.015,ntext,transform=ax.transAxes,va='bottom',fontsize=8,color=STYLE['neutral'])

    limits=[axes[0,j].get_ylim() for j in range(2)];lo=min(a for a,b in limits);hi=max(b for a,b in limits)
    for ax in axes[0]:ax.set_ylim(lo,hi)
    context_legend(fig)
    note(fig,'△ 高于历史第95百分位；▽ 低于历史第5百分位；辐射仅列原值。')
    return fig


def figure_multiscale(data,process,pca):
    ids=[r['id'] for r in data['primary']['candidate_ranking']];rows=[];texts=[]
    for g in ids:
        source=data['primary']['factors'][g]['windows'];vals=[];txt=[]
        for hours in WINDOWS:
            r=source.get(str(hours))
            if r is None or r.get('a_mean') is None:vals.append(np.nan);txt.append('不适');continue
            a=float(r['a_mean']);lo=float(r['q_min']);hi=float(r['q_max'])
            direction='高' if lo>.5 else '低' if hi<.5 else '中' if lo==hi==.5 else '异'
            vals.append(a);txt.append(f'{a:.2f}{direction}')
        rows.append(vals);texts.append(txt)
    fig=plt.figure(figsize=(156/25.4,127/25.4));ax=fig.add_axes([.285,.14,.55,.82]);cm=plt.get_cmap('cividis').copy();cm.set_bad(STYLE['missing'])
    vals=np.asarray(rows);mesh=ax.pcolormesh(np.arange(6)-.5,np.arange(len(ids)+1)-.5,np.ma.masked_invalid(vals),cmap=cm,vmin=0,vmax=1,shading='flat',linewidth=0,rasterized=False)
    ax.set_ylim(len(ids)-.5,-.5);ax.set_xlim(-.5,4.5);ax.set_yticks(range(len(ids)),[SHORT[g] for g in ids]);ax.set_xticks(range(5),['1','3','7','15','30']);ax.set_xlabel('窗口长度 / 日')
    ax.tick_params(length=0,pad=5);[s.set_visible(False) for s in ax.spines.values()]
    for i,row in enumerate(texts):
        for j,txt in enumerate(row):
            value=vals[i,j];ax.text(j,i,txt,ha='center',va='center',fontsize=8.5,color='white' if np.isfinite(value) and value<.49 else '#202428')
    cax=fig.add_axes([.855,.265,.016,.56]);cb=fig.colorbar(mesh,cax=cax,ticks=[0,.5,1]);cb.ax.set_title('a',fontsize=9,pad=6)
    if cb.solids is not None:cb.solids.set_rasterized(False)
    note(fig,'高/低：全部配置q在0.5同侧；异：跨侧；灰格为不适用，真实0值仍写数字。')
    return fig


def figure_ranking(data,process,pca):
    rows=data['primary']['candidate_ranking'];ids=[r['id'] for r in rows];values=np.array([r['score_S'] for r in rows])
    fig=plt.figure(figsize=(156/25.4,108/25.4));ax=fig.add_axes([.285,.16,.60,.79]);y=np.arange(len(rows));top=set(data['primary']['top_candidates'])
    for i,(g,s) in enumerate(zip(ids,values)):
        color=STYLE['accent'] if g in top else STYLE['target'];ax.hlines(i,0,s,color=STYLE['band'],lw=1.0)
        ax.plot(s,i,'o',color=color,ms=5.5 if g in top else 4.4);ax.text(1.025,i,f'{s:.3f}',ha='left',va='center',fontsize=9,clip_on=False)
    cutoff=data['primary'].get('top_cutoff_score_S')
    if cutoff is not None:ax.axvline(cutoff,color=STYLE['accent'],ls=(0,(3,3)),lw=.8)
    ax.set_yticks(y,[SHORT[g] for g in ids]);ax.set_ylim(len(rows)-.45,-.55);ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('校准历史位置 S')
    ax.grid(axis='x',color=STYLE['grid'],lw=.5);ax.tick_params(axis='y',length=0,pad=5)
    fig.legend(handles=[Line2D([],[],color=STYLE['accent'],marker='o',ls='none',label='主前五'),Line2D([],[],color=STYLE['target'],marker='o',ls='none',label='其余候选')],loc='upper right',bbox_to_anchor=(.98,1),ncol=2,fontsize=8.5,handletextpad=.25,columnspacing=.8)
    note(fig,f"去重后{len(rows)}项竞争；虚线为第五项分数，不是概率、显著性或灾害阈值。")
    return fig


def figure_pca(data,process,pca):
    fig,axes=grid(2,2,130,top=21,bottom=20,vgap=23)
    rng=np.random.default_rng(STYLE['display_seed']);years=list(pca['support']['valid_years'])
    for i,(hours,metric) in enumerate(((168,'T2'),(720,'T2'),(168,'Q'),(720,'Q'))):
        ax=axes.flat[i];scale=pca['scales'][str(hours)];summary=scale['summaries']['main'];agg=summary['metrics'][metric]
        for j,year in enumerate(years):
            conf=[c for c in scale['configurations'] if c['status']=='ok' and c['k_branch']=='main' and c['calibration_year']==year]
            vals=[c['positions'][metric] for c in conf];valid=[v for v in vals if v is not None]
            jitter=rng.uniform(-.15,.15,len(valid));ax.scatter(j+jitter,valid,s=12,color='#A4ABB0',alpha=.65,zorder=2,linewidths=0)
        for row in agg['calibration_years']:
            j=years.index(row['year']);ax.plot(j,row['position_mean'],'o',mfc='white',mec=STYLE['target'],mew=1,ms=5,zorder=4)
        if agg['position_mean'] is not None:ax.plot(len(years),agg['position_mean'],'D',ms=6,color=STYLE['accent'],zorder=5)
        ax.axhline(.5,color=STYLE['grid'],lw=.7);ax.set_ylim(0,1.03);ax.set_yticks([0,.5,1]);ax.set_xlim(-.45,len(years)+.5)
        ax.set_xticks(range(len(years)+1),[y[-2:] for y in years]+['汇总']);ax.set_xlabel('校准年（2021—2025）');ax.set_ylabel('留出历史中秩位置')
        kvals='/'.join(map(str,summary['k_values']));title=f"{hours//24}日 {'T²' if metric=='T2' else 'Q'}，k={kvals}"
        panel(ax,chr(97+i),title)
        text='NA' if agg['position_mean'] is None else f"汇总 {agg['position_mean']:.3f}"
        ax.text(0,1.015,text,transform=ax.transAxes,va='bottom',fontsize=8.5)
    fig.legend(handles=[Line2D([],[],marker='o',ms=3,color='#A4ABB0',ls='none',label='各参考配置'),Line2D([],[],marker='o',mfc='white',mec=STYLE['target'],ls='none',label='校准年均值'),Line2D([],[],marker='D',color=STYLE['accent'],ls='none',label='年度等权汇总')],loc='upper center',bbox_to_anchor=(.52,.995),ncol=3,fontsize=8.5,columnspacing=.8,handletextpad=.3)
    counts={(s['summaries']['main']['metrics'][m]['valid_configurations'],s['summaries']['main']['metrics'][m]['valid_calibration_years']) for s in pca['scales'].values() for m in ('T2','Q')}
    count_text=f'各项{next(iter(counts))[0]}配置/{next(iter(counts))[1]}校准年' if len(counts)==1 else '有效分母见结果表'
    note(fig,count_text+f"；横向抖动种子{STYLE['display_seed']}，纵值不变；0.5非显著性阈值。")
    return fig


def figure_soil(data,process,pca):
    fig,axes=grid(2,2,143,top=22,bottom=20,vgap=24);ch=channel_map(process);x=array(process['weather']['time_relative_hours'])/24
    depths=('01','02','04','08','16');depth_labels=('0.1 m','0.2 m','0.4 m','0.8 m','1.6 m');styles=('-','--','-.',':',(0,(5,1,1,1)))
    for row,(prefix,unit,title) in enumerate((('soilt','°C','五深度土温'),('vwc','m³/m³','五深度含水量'))):
        for col,start in enumerate((-30,-3)):
            ax=axes[row,col];background(ax,x,ch[prefix+'01'],start)
            for j,d in enumerate(depths):target(ax,x,ch[prefix+d],start,color=STYLE['colors'][j],linestyle=styles[j],bg=False,lw=1.35)
            ax.set_ylabel(unit);time_axes(ax,start);panel(ax,chr(97+row*2+col),title+('：30日' if col==0 else '：末72 h'))
        limits=[axes[row,j].get_ylim() for j in range(2)];lo=min(a for a,b in limits);hi=max(b for a,b in limits)
        for ax in axes[row]:ax.set_ylim(lo,hi)
    fig.legend(handles=[Line2D([],[],color=STYLE['colors'][j],ls=styles[j],label=lab,lw=1.35) for j,lab in enumerate(depth_labels)],loc='upper center',bbox_to_anchor=(.52,.995),ncol=5,fontsize=8.5,handlelength=1.3,columnspacing=.65,handletextpad=.25)
    counts30={support(ch[p+d],'30d') for p in ('soilt','vwc') for d in depths};counts72={support(ch[p+d],'72h') for p in ('soilt','vwc') for d in depths}
    text='各深度n30='+next(iter(counts30))+'，n72='+next(iter(counts72)) if len(counts30)==len(counts72)==1 else '逐深度支持见结果表'
    note(fig,'灰线/灰带仅为0.1 m历史背景；'+text+'。')
    return fig


def figure_channels(data,process,pca):
    fig,axes=grid(4,2,232,left=20,top=13,bottom=20,hgap=19,vgap=23);ch=channel_map(process);x=array(process['weather']['time_relative_hours'])/24
    plan=[(('pressure',),'气压原列','hPa'),(('rsd','rsu'),'短波原通道','W/m²'),
          (('rld','rlu'),'长波原通道','W/m²'),(('h',),'H原通道','W/m²'),
          (('le',),'LE原通道','W/m²'),(('ustar',),'USTAR原通道','m/s'),
          (('g01','g02'),'板热流原通道','W/m²')]
    for i,(keys,title,unit) in enumerate(plan):
        ax=axes.flat[i]
        for j,key in enumerate(keys):target(ax,x,ch[key],-30,color=STYLE['colors'][j],label=key.upper() if key not in ('g01','g02') else '0.1 m' if key=='g01' else '0.2 m',linestyle='-' if j==0 else '--',bg=len(keys)==1,lw=1.35)
        ntext='；'.join((key.upper()+('实 ' if j==0 else '虚 ') if len(keys)>1 else 'n=')+support(ch[key],'30d') for j,key in enumerate(keys))
        ax.text(0,1.015,ntext,transform=ax.transAxes,va='bottom',fontsize=8,color=STYLE['neutral'])
        ax.set_ylabel(unit);time_axes(ax,-30);panel(ax,chr(97+i),title)
    ax=axes.flat[-1];ax.set_axis_off()
    ax.text(0,1,'全部保持有符号原值\n缺测处断线，不补零\n多通道面板只列目标\n单通道灰带为历史四分位\nH/LE/USTAR无源QC列\n不作净收支或蒸散换算',transform=ax.transAxes,va='top',fontsize=9,linespacing=1.7)
    note(fig,'各通道支持以有效/应有小时列示；UTC+8截止与历史背景口径同正文。')
    return fig


def figure_contributions(data,process,pca):
    ids=pca['feature_ids'];columns=((168,'T2'),(168,'Q'),(720,'T2'),(720,'Q'));values=np.full((len(ids),4),np.nan);labels=[]
    for j,(hours,metric) in enumerate(columns):
        rows={r['id']:r for r in pca['scales'][str(hours)]['summaries']['main']['contributions'][metric]}
        column=[]
        for i,g in enumerate(ids):
            r=rows[g];q=r['position_mean'];share=r['normalized_share_mean'];values[i,j]=q if q is not None else np.nan
            if q is None or share is None:column.append('NA');continue
            rank=str(r['rank_min']) if r['rank_min']==r['rank_max'] else f"{r['rank_min']}–{r['rank_max']}"
            column.append(f"{share*100:.1f}%/{q:.2f}\n秩{rank}")
        labels.append(column)
    fig=plt.figure(figsize=(156/25.4,160/25.4));ax=fig.add_axes([.285,.18,.55,.76]);cm=plt.get_cmap('cividis').copy();cm.set_bad(STYLE['missing'])
    mesh=ax.pcolormesh(np.arange(5)-.5,np.arange(len(ids)+1)-.5,np.ma.masked_invalid(values),cmap=cm,vmin=0,vmax=1,shading='flat',rasterized=False,linewidth=0)
    ax.set_ylim(len(ids)-.5,-.5);ax.set_xlim(-.5,3.5);ax.set_yticks(range(len(ids)),[SHORT[g] for g in ids]);ax.set_xticks(range(4),['7日 cT/T²','7日 cQ/Q','30日 cT/T²','30日 cQ/Q']);ax.xaxis.tick_top();ax.tick_params(length=0,pad=5,labelsize=8.5);[s.set_visible(False) for s in ax.spines.values()]
    for j in range(4):
        for i,txt in enumerate(labels[j]):ax.text(j,i,txt,ha='center',va='center',fontsize=8,color='white' if np.isfinite(values[i,j]) and values[i,j]<.49 else '#202428',linespacing=1.25)
    cax=fig.add_axes([.855,.30,.016,.50]);cb=fig.colorbar(mesh,cax=cax,ticks=[0,.5,1]);cb.ax.set_title('q',fontsize=9,pad=6)
    if cb.solids is not None:cb.solids.set_rasterized(False)
    note(fig,'颜色为贡献历史位置q；格内：份额均值/q，第二行为配置名次范围。\n负cT不截断，份额不等于因果占比；整体联合位置仍以正文图5为准。')
    return fig


FIGURES={'q1_fig01_structure':figure_structure,'q1_fig02_process':figure_process,
         'q1_fig03_multiscale':figure_multiscale,'q1_fig04_ranking':figure_ranking,'q1_fig05_pca':figure_pca,
         'q1_app01_soil':figure_soil,'q1_app02_channels':figure_channels,'q1_app03_contributions':figure_contributions}


def build_figures(data,process,pca,only=None):
    validate(data,process,pca);wanted=set(only or FIGURES)
    if wanted-set(FIGURES):raise ValueError('未知图名。')
    return {name:fn(data,process,pca) for name,fn in FIGURES.items() if name in wanted}


def require_matplotlib_panel_alignment(fig, tolerance_pt=1.5):


    fig.canvas.draw()
    width_pt,height_pt=fig.get_size_inches()*72
    boxes=[np.array([a.get_position().x0*width_pt,a.get_position().y0*height_pt,
                     a.get_position().x1*width_pt,a.get_position().y1*height_pt])
           for a in fig.axes if a.axison and a.get_position().width>.10]
    for i,a in enumerate(boxes):
        for b in boxes[i+1:]:
            if abs(a[1]-b[1])<tolerance_pt:
                if abs(a[3]-b[3])>tolerance_pt or abs((a[2]-a[0])-(b[2]-b[0]))>tolerance_pt:
                    raise ValueError('同排证据面板的高度或宽度不一致。')
            if abs(a[0]-b[0])<tolerance_pt and abs(a[2]-b[2])>tolerance_pt:
                raise ValueError('同列证据面板没有对齐。')
    for item in fig.findobj(matplotlib.text.Text):
        if item.get_visible() and item.get_text() and item.get_fontsize()<STYLE['small']:
            raise ValueError('存在低于约定8 pt的绘图文字。')


def save_figure(fig,path,dpi=400,alignment_audit=None):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig)
    if alignment_audit is not None:alignment_audit(fig)
    path=Path(path)
    fig.savefig(path.with_suffix('.png'),dpi=dpi)
    fig.savefig(path.with_suffix('.pdf'),metadata={'Creator':'Matplotlib','Author':'','CreationDate':None,'ModDate':None})
    fig.savefig(path.with_suffix('.svg'),metadata={'Creator':'Matplotlib','Date':None})


def main(argv=None):
    here=Path(__file__).resolve().parent;parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir',type=Path,default=here/'results');parser.add_argument('--output-dir',type=Path,default=here/'figures')
    parser.add_argument('--only',nargs='+',choices=list(FIGURES));parser.add_argument('--dpi',type=int,default=STYLE['dpi'])
    args=parser.parse_args(argv)
    if args.dpi<300:parser.error('PNG分辨率至少300 dpi。')
    data=load_json(args.results_dir/'analysis_results.json.gz');process=load_json(args.results_dir/'process_results.json.gz');pca=load_json(args.results_dir/'pca_results.json.gz')
    font=configure_style();args.output_dir.mkdir(parents=True,exist_ok=True)
    for name,fig in build_figures(data,process,pca,args.only).items():save_figure(fig,args.output_dir/name,args.dpi);plt.close(fig)
    print(f'已生成{len(args.only or FIGURES)}张图，宽156 mm；PNG {args.dpi} dpi，PDF/SVG；中文字体：{font}。')
    return 0


SHORT.update({'d1_era5_precip_sum':'格点累计总降水','d1_era5_precip_max':'格点最大小时降水',
              'd1_era5_rain_hours':'格点有降水小时','d1_era5_api':'格点前期降水指数'})


def figure_precipitation(record):

    fig,axes=grid(3,1,169,left=17,right=5,top=21,bottom=25,vgap=23)
    x=array(record['time_relative_hours'])/24
    values=array(record['values'])
    ymax=max(float(np.nanmax(values)),float(np.nanmax(array(record['background_q95']))))*1.12
    for i,(key,start,title) in enumerate((('30d',-30,'30日过程'),('7d',-7,'末7日'),('72h',-3,'末72小时'))):
        ax=axes[i,0];mask=(x>start)&(x<=0);xx=x[mask]
        lo=array(record['background_q25'])[mask];hi=array(record['background_q75'])[mask]
        ax.fill_between(xx,lo,hi,color=STYLE['band'],lw=0)
        ax.plot(xx,array(record['background_median'])[mask],color=STYLE['neutral'],ls='--',lw=.8)
        ax.plot(xx,values[mask],color=STYLE['target'],lw=1.05)
        exceed=np.asarray(record['high_p95_exceedance'],dtype=bool)[mask]
        ax.scatter(xx[exceed],values[mask][exceed],marker='^',s=14,color=STYLE['accent'],zorder=5)
        s=record['summaries'][key]
        panel(ax,chr(97+i),title)
        ax.text(0,1.015,f"完整{s['available_hours']}/{s['expected_hours']} h；累计{s['sum_precipitation_mm']:.1f} mm；高于历史第95百分位{s['high_p95_hours']} h",transform=ax.transAxes,fontsize=8,va='bottom')
        ax.set_ylim(0,ymax);ax.set_ylabel('小时总降水 / mm');time_axes(ax,start)
    handles=[Line2D([],[],color=STYLE['target'],lw=1.1,label='目标格点'),
             Line2D([],[],color=STYLE['neutral'],ls='--',lw=.8,label='历史中位数'),
             Patch(facecolor=STYLE['band'],label='历史25%—75%'),
             Line2D([],[],color=STYLE['accent'],marker='^',ls='none',ms=4,label='高于历史第95百分位')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.995),ncol=2,fontsize=8,handlelength=1.5,columnspacing=.9)
    note(fig,'固定ERA5格点总降水（含雨雪），各图采用同一量程。\n历史按年等权；重叠轨迹不是独立事件。格点结果不代表源坡实测。',y=.008)
    return fig


def enhanced_label(feature_id,features):
    if feature_id in SHORT:return SHORT[feature_id]
    return features.get(feature_id,{}).get('label',feature_id).replace('ERA5','格点')


def figure_enhanced_ranking(enhanced):
    rows=enhanced['primary']['candidate_ranking'];features={f['id']:f for f in enhanced['features']}
    fig=plt.figure(figsize=(156/25.4,139/25.4));ax=fig.add_axes([.34,.15,.53,.78])
    top=set(enhanced['primary']['top_candidates'])
    for i,row in enumerate(rows):
        g=row['id'];s=row['score_S'];color=STYLE['accent'] if g in top else STYLE['target']
        ax.hlines(i,0,s,color=STYLE['band'],lw=1)
        ax.plot(s,i,'s' if g.startswith('d1_') else 'o',ms=5,color=color)
        ax.text(1.035,i,f'{s:.3f}',va='center',ha='left',fontsize=8.5,clip_on=False)
    ax.set_yticks(range(len(rows)),[enhanced_label(r['id'],features) for r in rows]);ax.tick_params(axis='y',length=0,pad=4)
    ax.set_xlim(0,1);ax.set_ylim(len(rows)-.4,-.6);ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('校准历史位置 S')
    ax.grid(axis='x',color=STYLE['grid'],lw=.5)
    fig.legend(handles=[Line2D([],[],color=STYLE['target'],marker='o',ls='none',label='原站通道或代理'),
                        Line2D([],[],color=STYLE['target'],marker='s',ls='none',label='ERA5格点代理'),
                        Line2D([],[],color=STYLE['accent'],marker='o',ls='none',label='主前五')],
               loc='upper center',bbox_to_anchor=(.53,.995),ncol=3,fontsize=8.5,handletextpad=.3,columnspacing=.8)
    note(fig,f'加入格点降水后的去重池：{len(rows)}项；S为历史相对位置，不是灾害概率。\n站点与格点空间支持不同；方形标识不能解释成新增站点观测。',y=.012)
    return fig


def figure_enhanced_multiscale(enhanced):
    rows=enhanced['primary']['candidate_ranking'];features={f['id']:f for f in enhanced['features']}
    ids=[r['id'] for r in rows];vals=[];texts=[]
    for g in ids:
        v=[];t=[]
        for h in WINDOWS:
            r=enhanced['primary']['factors'][g]['windows'].get(str(h),{})
            a=r.get('a_mean')
            if a is None:v.append(np.nan);t.append('不适');continue
            lo,hi=r['q_min'],r['q_max'];direction='高' if lo>.5 else '低' if hi<.5 else '中' if lo==hi==.5 else '异'
            v.append(a);t.append(f'{a:.2f}{direction}')
        vals.append(v);texts.append(t)
    fig=plt.figure(figsize=(156/25.4,145/25.4));ax=fig.add_axes([.34,.13,.48,.82]);cm=plt.get_cmap('cividis').copy();cm.set_bad(STYLE['missing'])
    z=np.asarray(vals);mesh=ax.pcolormesh(np.arange(6)-.5,np.arange(len(ids)+1)-.5,np.ma.masked_invalid(z),cmap=cm,vmin=0,vmax=1,shading='flat')
    ax.set_xlim(-.5,4.5);ax.set_ylim(len(ids)-.5,-.5);ax.set_yticks(range(len(ids)),[enhanced_label(g,features) for g in ids]);ax.set_xticks(range(5),['1','3','7','15','30']);ax.set_xlabel('窗口长度 / 日');ax.tick_params(length=0,pad=5)
    for spine in ax.spines.values():spine.set_visible(False)
    for i,row in enumerate(texts):
        for j,t in enumerate(row):ax.text(j,i,t,ha='center',va='center',fontsize=8,color='white' if np.isfinite(z[i,j]) and z[i,j]<.49 else '#202428')
    cax=fig.add_axes([.86,.28,.018,.56]);cb=fig.colorbar(mesh,cax=cax,ticks=[0,.5,1]);cb.ax.set_title('a',fontsize=9,pad=6)
    if cb.solids is not None:cb.solids.set_rasterized(False)
    note(fig,'高/低：所有参考配置的q同侧；异：跨侧。\n逐窗双侧异常a与多尺度校准S分开，灰格为该尺度不适用。',y=.01)
    return fig


def figure_method_comparison(enhanced):
    rows=enhanced['comparisons'];features={f['id']:f for f in enhanced['features']}
    fig=plt.figure(figsize=(156/25.4,133/25.4));ax=fig.add_axes([.34,.13,.50,.81])
    z=np.asarray([[r[k]['rank_min'] for k in ('E2_rank','B_rank','C_rank')] for r in rows],dtype=float)
    mesh=ax.pcolormesh(np.arange(4)-.5,np.arange(len(rows)+1)-.5,z,cmap='cividis_r',vmin=1,vmax=len(rows),shading='flat')
    ax.set_xlim(-.5,2.5);ax.set_ylim(len(rows)-.5,-.5);ax.set_yticks(range(len(rows)),[enhanced_label(r['id'],features) for r in rows]);ax.set_xticks(range(3),['多尺度','熵权TOPSIS','固定7日']);ax.xaxis.tick_top();ax.tick_params(length=0,pad=6)
    for s in ax.spines.values():s.set_visible(False)
    for i in range(len(rows)):
        for j in range(3):ax.text(j,i,str(int(z[i,j])),ha='center',va='center',fontsize=8.5,color='white' if z[i,j]>=8 else '#202428')
    cax=fig.add_axes([.88,.28,.018,.56]);cb=fig.colorbar(mesh,cax=cax,ticks=[1,5,10,15]);cb.ax.set_title('名次',fontsize=8,pad=6)
    if cb.solids is not None:cb.solids.set_rasterized(False)
    note(fig,'共同15项原因素及降水代理；两项高度差代理没有原机制等级，未列入。\nTOPSIS在共同集重新拟合；名次重合不表示预测准确率。',y=.01)
    return fig


def enhanced_main(argv=None):
    here=Path(__file__).resolve().parent;parser=argparse.ArgumentParser(description='Q1 D0 + fixed ERA5 result-only figures')
    parser.add_argument('--results-dir',type=Path,default=here/'results');parser.add_argument('--output-dir',type=Path,default=here/'figures')
    parser.add_argument('--only',nargs='+');parser.add_argument('--dpi',type=int,default=400)
    args=parser.parse_args(argv)
    if args.dpi<300:parser.error('PNG分辨率至少300 dpi。')
    baseline=args.results_dir/'baseline';data=load_json(baseline/'analysis_results.json.gz');process=load_json(baseline/'process_results.json.gz');pca=load_json(baseline/'pca_results.json.gz')
    validate(data,process,pca);e=load_json(args.results_dir/'enhanced_results.json.gz');rain=load_json(args.results_dir/'d1_process.json')
    if rain['metadata']['primary_cutoff']!=process['metadata']['primary_cutoff']:raise ValueError('D0/D1截止不一致')
    epca=e.get('pca');functions={
        'q1_fig01_structure':lambda:figure_structure(data,process,pca),
        'q1_fig02_process':lambda:figure_process(data,process,pca),
        'q1_fig03_precipitation':lambda:figure_precipitation(rain),
        'q1_fig04_ranking':lambda:figure_enhanced_ranking(e),
        'q1_app01_multiscale':lambda:figure_enhanced_multiscale(e),
        'q1_app02_soil':lambda:figure_soil(data,process,pca),
        'q1_app03_channels':lambda:figure_channels(data,process,pca),
        'q1_app04_methods':lambda:figure_method_comparison(e)}
    if epca and epca.get('scales'):functions['q1_fig05_pca']=lambda:figure_pca(e,process,epca)
    chosen=args.only or list(functions)
    if set(chosen)-set(functions):parser.error('未知或当前不可计算图名：'+','.join(set(chosen)-set(functions)))
    font=configure_style();args.output_dir.mkdir(parents=True,exist_ok=True)
    audit_hook=globals().get('_external_alignment_hook')
    for name in chosen:
        fig=functions[name]()
        hook=(lambda fig,n=name:audit_hook(fig,n)) if audit_hook else None
        save_figure(fig,args.output_dir/name,args.dpi,alignment_audit=hook);plt.close(fig)
    print(f'实际结果绘图完成：{len(chosen)}图；156 mm，PNG{args.dpi} dpi与PDF/SVG；字体{font}')
    return 0


if __name__=='__main__':raise SystemExit(enhanced_main())

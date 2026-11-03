from pathlib import Path
import argparse,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sol():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
 inputs=[];outputs=[]
 def read(path):
  p=args.root/path;inputs.append({'path':path,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return p
 def save(fig,name):
  for ext in ['pdf','svg','png']:
   p=args.out/(name+'.'+ext);fig.savefig(p,dpi=600,bbox_inches='tight');outputs.append({'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  plt.close(fig)
 searches=json.loads(read('results/critical/CRITICAL_SEARCH.json').read_text());models=json.loads(read('results/model_seal/METHODS.json').read_text())
 colors=['#777777','#4477AA','#66AABB','#AA3377'];labels={'A_single_fill':'Single GP','B_MF_fill':'MF passive','Bm_MF_count_matched':'MF matched','C_MF_active':'MF active'}
 for search in searches:
  fig,ax=plt.subplots(figsize=(7.2,4.2),layout='constrained');cut=search['slice'];mu=.03+.12*np.arange(129)/128
  for meta,color in zip(models,colors):
   f=np.load(read('results/model_seal/'+meta['folder']+'/TEST_PREDICTIONS_FROZEN.npz'));ax.plot(mu,f['mean_g'][16+cut*129:16+(cut+1)*129],color=color,lw=1.2,label=labels[meta['method']])
  valid=[r for r in search['sequence'] if r['valid']];xs=[.03+.12*r['j']/128 for r in valid];ys=[np.log1p((r['h_m']-.1)/.2) for r in valid]
  ax.scatter(xs,ys,s=26,color='black',zorder=5,label='Queried 30 m PDE')
  interval=search['final_mu_interval']
  if interval:
   for value in interval:ax.axvline(value,color='#999999',ls=':',lw=.8)
  ax.axhline(0,color='#999999',ls='--',lw=.8);ax.set(xlabel='Coulomb coefficient μ',ylabel='Log depth margin g',xlim=(.025,.155));ax.legend(loc='upper right',fontsize=8)
  save(fig,'critical_slice'+str(cut))
 f=np.load(read('results/full_diagnostics/FULL_COMMON_ENERGY_CURVES.npz'));fig,ax=plt.subplots(figsize=(7.2,4.2),layout='constrained')
 reference=8825015023.410414
 for name,color in zip(['copernicus_30','copernicus_40','aster_40','double_40'],colors):
  a=f[name];ax.plot(a[:,0],100*a[:,3]/reference,lw=1.2,label=name.replace('_',' '),color=color)
 ax.axvline(960,color='#999999',ls=':',lw=.8);ax.set(xlabel='Time after prescribed entry (s)',ylabel='Energy ledger residual / fixed input (%)');ax.legend(loc='upper left',fontsize=8)
 save(fig,'full_energy_normalized')
 (args.out/'PLOT_RECEIPT.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':inputs,'outputs':outputs,'PDE_ran':False,'GP_fitted':False,'alignment':'single plot-area per figure; not applicable'},indent=2))
if __name__=='__main__':sol()

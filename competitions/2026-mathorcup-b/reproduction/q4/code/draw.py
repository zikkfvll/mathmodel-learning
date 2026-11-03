from pathlib import Path
import argparse,csv,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sol():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 p=a.root/'results/PUBLISH_COMPARISON.csv';rows=list(csv.DictReader(p.open()))
 plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
 names={'arrived':'Arrived','valid_no_arrival_by_horizon':'No arrival\nby 1800 s','unknown':'Unknown'}
 cells=[]
 for r in rows:
  cost=f"{float(r['new_PDE_child_wall_s']):.2f}" if r['new_PDE_child_wall_s'] else 'Cached history'
  cells.append([r['request'],f"{float(r['mu']):.5f}",names[r['GP_estimate']],names[r['physical_status']],names[r['published_status']],cost])
 fig,ax=plt.subplots(figsize=(9,3.2));ax.axis('off')
 table=ax.table(cellText=cells,colLabels=['Request','μ','GP estimate','Physical model','Published','Fresh PDE (s)'],cellLoc='center',colWidths=[.23,.1,.15,.17,.17,.18],bbox=[0,.16,1,.74]);table.auto_set_font_size(False);table.set_fontsize(8)
 for (r,c),cell in table.get_celld().items():
  cell.set_linewidth(.4);cell.set_edgecolor('#BBBBBB')
  if r==0:cell.set_facecolor('#EAF0F5');cell.set_text_props(weight='bold')
 ax.set_title('Physical verification controls publication; GP estimates remain unchanged',fontsize=10,pad=8)
 fig.text(.09,.065,'Four historical replays + two new PDE queries. This is interface consistency, not GP accuracy or natural-event validation.',fontsize=8)
 outputs=[]
 for ext in ['pdf','svg','png']:
  dest=a.out/('publication_guard.'+ext);fig.savefig(dest,dpi=600,bbox_inches='tight');outputs.append({'path':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
 plt.close(fig)
 (a.out/'PLOT_RECEIPT.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'outputs':outputs,'PDE_ran':False,'GP_fitted':False,'alignment':'single table plot-area'},indent=2))
if __name__=='__main__':sol()

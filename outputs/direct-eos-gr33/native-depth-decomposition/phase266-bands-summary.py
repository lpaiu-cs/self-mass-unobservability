import json
from pathlib import Path
R=Path('/home/lpaiu/work/native-retained-tail-runtime')
a=json.loads((R/'.phase266-depth-b128.json').read_text());b=json.loads((R/'.phase266-depth-b64.json').read_text())
print('errors',a.get('error'),b.get('error'))
ra={r['band']:r for r in a['rows']};rb={r['band']:r for r in b['rows']}
total=ra['all']['charge'];total64=rb['all']['charge']
parts=[k for k in ra if k!='all']
s=sum(ra[k]['charge'] for k in parts);s64=sum(rb[k]['charge'] for k in parts)
print('128 all %.10e sum %.10e rel %.2e'%(total,s,(s-total)/total))
print('64  all %.10e sum %.10e rel %.2e'%(total64,s64,(s64-total64)/total64))
print('%-16s %-22s %14s %9s %14s %10s'%('band','depth_km','charge128','share','charge64','64vs128'))
for k in parts:
    x=ra[k];y=rb[k];dep=x['depth_km']
    rel=(y['charge']-x['charge'])/x['charge'] if x['charge'] else float('nan')
    print('%-16s %-22s %14.6e %9.4f %14.6e %10.2e'%(k,('%.1f-%.1f'%(dep[0],dep[1])) if dep else '-',x['charge'],x['charge']/total,y['charge'],rel))
print('seconds per band ~',round(sum(r['setup_seconds']+r['propagate_seconds'] for r in a['rows'])/len(a['rows']),1))

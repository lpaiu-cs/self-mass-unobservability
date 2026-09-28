"""Phase269: PL-off vs PL/MHD compact charge on the same 2x grid and equation (registered rule), run summary, depth groups.

Rule (REQUEST269): the conditional negative charge is kept under this EOS/optics choice if the PL-off endpoint charge is
negative; the magnitude change against the PL/MHD 2x solution is reported separately (resolution context: 2x->4x -1.17%).
Usage: python3 compare269.py   (run in WSL; writes readout269-compare.json in the PL-off runtime)
"""
import glob, json, os
RE = '/home/lpaiu/work/native-eos269-runtime/'; R2 = '/home/lpaiu/work/native-refined267-runtime/'
ROOT = '/mnt/e/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/'
p267 = json.load(open(ROOT + 'native-refined-primary/result.json')); p268 = json.load(open(ROOT + 'native-quad-refined-primary/result.json'))
rows = {}
for n in ['61', '62', '63', '64']:
    f = RE + f'readout269-eos{n}-work/readout-field.json'
    if not os.path.exists(f): continue
    d = json.load(open(f)); qe = d['endpoint_compact_charge']; q2 = p267['charges_refined'][n]; q4 = p268['rows'][n]['q4']
    same_time = open(RE + f'.phase269-readout-eos{n}-time.txt').read().split()[-1] == 'True'
    rows[n] = dict(q_ploff_2x=qe, q_plmhd_2x=q2, q_plmhd_4x=q4, same_time_as_2x=same_time, evaluation_time=d['evaluation_time'],
                   magnitude_change_vs_plmhd_2x=(abs(qe) - abs(q2))/abs(q2), relative_change_vs_plmhd_2x=(qe - q2)/abs(q2),
                   inexact_inverse_max=d.get('inexact_inverse_max'))
decided = '64' in rows; negative = rows['64']['q_ploff_2x'] < 0 if decided else None
W = RE + 'primary269-eos64-work/'; seg = []; seconds = 0.; fallbacks = 0; lin = []; nonlin = []
for f in sorted(glob.glob(W + 'seg-*-driver.json'), key=lambda p: int(p.split('seg-')[1].split('-')[0])):
    d = json.load(open(f)); r = d.get('row', {}); seconds += d['seconds']; fallbacks += d.get('fallbacks') or 0
    lin += d.get('vector_gate_exceptions') or []; nonlin += d.get('nonlinear_exceptions') or []
    seg.append(dict(segment=os.path.basename(f)[:6], passed=d['passed'], seconds=round(d['seconds']), actual_steps=r.get('actual_completed_steps'),
                    linear=r.get('linear_relative'), energy=r.get('energy_balance_relative'), krylov=r.get('max_Krylov_iterations'),
                    fallbacks=d.get('fallbacks'), linear_exceptions=len(d.get('vector_gate_exceptions') or []), nonlinear_exceptions=len(d.get('nonlinear_exceptions') or [])))
exc = lambda es, k: max((e[k] for e in es), default=None)
run = dict(segments=seg, total_seconds=seconds, total_fallbacks=fallbacks,
           linear_exceptions=dict(count=len(lin), max_vector=exc(lin, 'vector_relative'), max_physical=exc(lin, 'physical_moment_max'), max_material=exc(lin, 'material_component_max')),
           nonlinear_exceptions=dict(count=len(nonlin), max_defect=exc(nonlin, 'vector_defect'), max_physical=exc(nonlin, 'physical_moment_max'),
                                     max_material=exc(nonlin, 'material_component_max'), rows=nonlin))
depth = None; DE = RE + 'readout269-eos64-work/depth-bands.json'
if os.path.exists(DE):
    be = {x['band']: x['charge'] for x in json.load(open(DE))['rows']}
    b2 = {x['band']: x['charge'] for x in json.load(open(R2 + 'readout267-refined64-work/depth-bands.json'))['rows']}
    two = lambda b, c: b[f'cells:{8+2*(c-8)}-{8+2*(c-8)}'] + b[f'cells:{9+2*(c-8)}-{9+2*(c-8)}']
    depth = dict(all=[b2['all'], be['all']], closure_ploff=(sum(v for k, v in be.items() if k != 'all') - be['all'])/abs(be['all']), **{'0-7': [b2['cells:0-7'], be['cells:0-7']]})
    for c in range(8, 16): depth[str(c)] = [two(b2, c), two(be, c)]
    for c in range(16, 19): depth[str(c)] = [b2[f'cells:{c+8}-{c+8}'], be[f'cells:{c+8}-{c+8}']]
    atm = ['cells:27-154', 'cells:155-282', 'cells:283-410', 'cells:411-538']
    depth['atmosphere'] = [sum(b2[k] for k in atm), sum(be[k] for k in atm)]; depth['boundary'] = [b2['boundary'], be['boundary']]
out = dict(classification='Counterexample candidate', rule='kept if the PL-off endpoint charge is negative; magnitude change vs PL/MHD 2x reported',
           decided=decided, endpoint_negative_ploff=negative, rows=rows, run=run, depth_T=depth)
open(RE + 'readout269-compare.json', 'w').write(json.dumps(out, indent=1) + '\n')
for n, r in rows.items():
    print('t%s PL-off 2x %.6e  PL/MHD 2x %.6e  (4x %.6e)  magnitude %+.3f%%  same time %s' % (n, r['q_ploff_2x'], r['q_plmhd_2x'], r['q_plmhd_4x'],
          100*r['magnitude_change_vs_plmhd_2x'], r['same_time_as_2x']))
print('decided', decided, 'endpoint negative', negative)
print('run seconds %.0f fallbacks %d linear exceptions %s nonlinear exceptions %d' % (seconds, fallbacks, run['linear_exceptions'], len(nonlin)))
if depth:
    for k, v in depth.items():
        if k != 'closure_ploff': print('depth %-10s PL/MHD %.4e  PL-off %.4e' % (k, v[0], v[1]))
    print('depth closure PL-off %.1e' % depth['closure_ploff'])

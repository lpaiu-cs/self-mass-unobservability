"""Phase268: 1x/2x/4x compact charge comparison at t61..T under the registered rule, plus run summary and depth groups.

Rule (REQUEST268): the 2x result is resolution-converged at the 2% level if the 2x->4x magnitude change is <= 2% at both
t61 and T. Observed order p = log2(|q2-q1|/|q4-q2|) and the Richardson value q4 + (q4-q2)/(2^p-1) are reported; they
are meaningful only when q2-q1 and q4-q2 have the same sign (monotone refinement sequence), which is flagged.
Usage: python3 compare268.py   (run in WSL; writes readout268-compare.json in the 4x runtime)
"""
import glob, json, math, os
R4 = '/home/lpaiu/work/native-refined268-runtime/'; R2 = '/home/lpaiu/work/native-refined267-runtime/'
PUB = '/mnt/e/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672/outputs/direct-eos-gr33/native-refined-primary/result.json'
pub = json.load(open(PUB)); q1s, q2s = pub['charges_original'], pub['charges_refined']
rows = {}
for n in ['61', '62', '63', '64']:
    f = R4 + f'readout268-quad{n}-work/readout-field.json'
    if not os.path.exists(f): continue
    d = json.load(open(f)); q4 = d['endpoint_compact_charge']; q1, q2 = q1s[n], q2s[n]
    same_time = open(R4 + f'.phase268-readout-quad{n}-time.txt').read().split()[-1] == 'True'
    d21, d42 = q2 - q1, q4 - q2; monotone = d21*d42 > 0
    p = math.log2(abs(d21)/abs(d42)) if d42 else math.inf
    rich = q4 + d42/(2**p - 1) if d42 and p > 0 else None
    rows[n] = dict(q1=q1, q2=q2, q4=q4, same_time_as_2x=same_time, evaluation_time=d['evaluation_time'],
                   change_2x_to_4x=(abs(q4) - abs(q2))/abs(q2), change_1x_to_4x=(abs(q4) - abs(q1))/abs(q1),
                   monotone=monotone, observed_order=p, richardson=rich, inexact_inverse_max=d.get('inexact_inverse_max'))
decided = all(n in rows for n in ['61', '64'])
converged = decided and all(abs(rows[n]['change_2x_to_4x']) <= 0.02 for n in ['61', '64'])
negative = rows.get('64', {}).get('q4', 0) < 0
# run summary (same fields as the phase-267 summary)
W = R4 + 'primary268-quad64-work/'; seg = []; seconds = 0.; fallbacks = 0; lin = []; nonlin = []
for f in sorted(glob.glob(W + 'seg-*-driver.json'), key=lambda p: int(p.split('seg-')[1].split('-')[0])):
    d = json.load(open(f)); r = d.get('row', {}); seconds += d['seconds']; fallbacks += d.get('fallbacks') or 0
    lin += d.get('vector_gate_exceptions') or []; nonlin += d.get('nonlinear_exceptions') or []
    seg.append(dict(segment=os.path.basename(f)[:6], passed=d['passed'], seconds=round(d['seconds']), actual_steps=r.get('actual_completed_steps'),
                    linear=r.get('linear_relative'), energy=r.get('energy_balance_relative'), krylov=r.get('max_Krylov_iterations'),
                    fallbacks=d.get('fallbacks'), linear_exceptions=len(d.get('vector_gate_exceptions') or []), nonlinear_exceptions=len(d.get('nonlinear_exceptions') or [])))
exc = lambda es, k: max((e[k] for e in es), default=None)
run = dict(segments=seg, total_seconds=seconds, total_fallbacks=fallbacks,
           linear_exceptions=dict(count=len(lin), max_vector=exc(lin, 'vector_relative'), max_physical=exc(lin, 'physical_moment_max'), max_material=exc(lin, 'material_component_max')),
           nonlinear_exceptions=dict(count=len(nonlin), max_defect=exc(nonlin, 'vector_defect') if nonlin and 'vector_defect' in nonlin[0] else None, rows=nonlin))
# depth at T grouped by original cell: 1x cell c <-> 2x cells 8+2(c-8), 9+2(c-8) <-> 4x cells 8+4(c-8) .. 11+4(c-8)
depth = None; D4 = R4 + 'readout268-quad64-work/depth-bands.json'
if os.path.exists(D4):
    b4 = {x['band']: x['charge'] for x in json.load(open(D4))['rows']}
    b2 = {x['band']: x['charge'] for x in json.load(open(R2 + 'readout267-refined64-work/depth-bands.json'))['rows']}
    b1 = {x['band']: x['charge'] for x in json.load(open('/home/lpaiu/work/native-retained-tail-runtime/.phase266-depth-b64.json'))['rows']}
    cells = lambda b, a, z: sum(b[f'cells:{i}-{i}'] for i in range(a, z + 1))
    depth = dict(all=[b1['all'], b2['all'], b4['all']],
                 closure_4x=(sum(v for k, v in b4.items() if k != 'all') - b4['all'])/abs(b4['all']))
    depth['0-7'] = [cells(b1, 0, 7), b2['cells:0-7'], b4['cells:0-7']]
    for c in range(8, 16):
        depth[str(c)] = [b1[f'cells:{c}-{c}'], cells(b2, 8 + 2*(c - 8), 9 + 2*(c - 8)), cells(b4, 8 + 4*(c - 8), 11 + 4*(c - 8))]
        depth[f'{c} 4x sub-cells'] = [b4[f'cells:{i}-{i}'] for i in range(8 + 4*(c - 8), 12 + 4*(c - 8))]
    for c in range(16, 19): depth[str(c)] = [b1[f'cells:{c}-{c}'], b2[f'cells:{c+8}-{c+8}'], b4[f'cells:{c+24}-{c+24}']]
    depth['atmosphere'] = [sum(b1[k] for k in ['cells:19-146', 'cells:147-274', 'cells:275-402', 'cells:403-530']),
                           sum(b2[k] for k in ['cells:27-154', 'cells:155-282', 'cells:283-410', 'cells:411-538']),
                           sum(b4[k] for k in ['cells:43-170', 'cells:171-298', 'cells:299-426', 'cells:427-554'])]
    depth['boundary'] = [b1['boundary'], b2['boundary'], b4['boundary']]
out = dict(classification='Counterexample candidate', rule='2x converged at 2% if |2x->4x magnitude change| <= 0.02 at t61 and T',
           decided=decided, converged_2x_at_2pct=converged if decided else None, endpoint_negative_4x=negative if '64' in rows else None,
           rows=rows, run=run, depth_T=depth)
open(R4 + 'readout268-compare.json', 'w').write(json.dumps(out, indent=1) + '\n')
for n, r in rows.items():
    print('t%s q1 %.6e q2 %.6e q4 %.6e  2x->4x %+.4f%%  1x->4x %+.4f%%  monotone %s  p %.3f  Richardson %s  same time %s' % (
        n, r['q1'], r['q2'], r['q4'], 100*r['change_2x_to_4x'], 100*r['change_1x_to_4x'], r['monotone'], r['observed_order'],
        '%.6e' % r['richardson'] if r['richardson'] is not None else '-', r['same_time_as_2x']))
print('decided', decided, 'converged_2x_at_2pct', out['converged_2x_at_2pct'], 'endpoint negative', out['endpoint_negative_4x'])
print('run seconds %.0f fallbacks %d linear exceptions %s nonlinear exceptions %d' % (seconds, fallbacks, run['linear_exceptions'], len(nonlin)))
if depth:
    for k, v in depth.items():
        if k != 'closure_4x': print('depth %-16s %s' % (k, ' '.join('%.4e' % x for x in v)))
    print('depth closure 4x %.1e' % depth['closure_4x'])

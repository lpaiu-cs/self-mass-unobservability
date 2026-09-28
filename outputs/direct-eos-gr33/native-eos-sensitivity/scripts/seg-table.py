"""Per-segment rows of a segmented primary run next to a reference run. Usage: python3 seg-table.py <work> <reference work>"""
import glob, json, os, sys
W, W2 = sys.argv[1].rstrip('/') + '/', sys.argv[2].rstrip('/') + '/'
for f in sorted(glob.glob(W + 'seg-*-driver.json'), key=lambda p: int(p.split('seg-')[1].split('-')[0])):
    d = json.load(open(f)); r = d.get('row', {}); name = os.path.basename(f)
    e = json.load(open(W2 + name)) if os.path.exists(W2 + name) else {}
    print(name[:6], d['passed'], round(d['seconds']), '(ref %s)' % round(e.get('seconds', 0)), 'steps', r.get('actual_completed_steps'), 'fb', d.get('fallbacks'),
          '(ref %s)' % e.get('fallbacks'), 'lin', '%.1e' % r.get('linear_relative', 0), 'energy', '%.1e' % r.get('energy_balance_relative', 0),
          'krylov', r.get('max_Krylov_iterations'), 'exc', len(d.get('vector_gate_exceptions') or []), len(d.get('nonlinear_exceptions') or []))

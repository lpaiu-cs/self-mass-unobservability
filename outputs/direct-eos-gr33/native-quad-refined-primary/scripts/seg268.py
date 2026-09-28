"""Phase268: per-segment rows of the 4x run next to the 2x run (seconds, steps, fallbacks, exceptions, gates)."""
import glob, json, os
W4 = '/home/lpaiu/work/native-refined268-runtime/primary268-quad64-work/'; W2 = '/home/lpaiu/work/native-refined267-runtime/primary267-refined64-work/'
for f in sorted(glob.glob(W4 + 'seg-*-driver.json'), key=lambda p: int(p.split('seg-')[1].split('-')[0])):
    d = json.load(open(f)); r = d.get('row', {}); name = os.path.basename(f)
    e = json.load(open(W2 + name)) if os.path.exists(W2 + name) else {}
    print(name[:6], d['passed'], round(d['seconds']), '(2x %s)' % round(e.get('seconds', 0)), 'steps', r.get('actual_completed_steps'), 'fb', d.get('fallbacks'),
          'lin', '%.1e' % r.get('linear_relative', 0), 'energy', '%.1e' % r.get('energy_balance_relative', 0), 'krylov', r.get('max_Krylov_iterations'),
          'exc', len(d.get('vector_gate_exceptions') or []), len(d.get('nonlinear_exceptions') or []))

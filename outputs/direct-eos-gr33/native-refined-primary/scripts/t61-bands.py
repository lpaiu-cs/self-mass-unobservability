"""Phase267: depth bands at t61 on both grids, refined sub-cells grouped by the original cell they split."""
import json
R = '/home/lpaiu/work/native-refined267-runtime/readout267-refined61-work/depth-bands.json'
G = '/home/lpaiu/work/native-retained-tail-runtime/readout267-identity2-work/depth-t61.json'
r = {x['band']: x for x in json.load(open(R))['rows']}; g = {x['band']: x for x in json.load(open(G))['rows']}
qr, qg = r['all']['charge'], g['all']['charge']
print('t61 total: original %.6e refined %.6e relative change %+.4f' % (qg, qr, (qr - qg)/abs(qg)))
print('band sums: refined %.3e (closure %.1e)' % (sum(x['charge'] for k, x in r.items() if k != 'all'), (sum(x['charge'] for k, x in r.items() if k != 'all') - qr)/abs(qr)))
groups = [('0-7', ['cells:0-7'], [f'cells:{i}-{i}' for i in range(8)])] + \
         [(f'{c}', [f'cells:{8+2*(c-8)}-{8+2*(c-8)}', f'cells:{9+2*(c-8)}-{9+2*(c-8)}'], [f'cells:{c}-{c}']) for c in range(8, 16)] + \
         [(f'{c}', [f'cells:{c+8}-{c+8}'], [f'cells:{c}-{c}']) for c in range(16, 19)] + \
         [('atm', ['cells:27-154', 'cells:155-282', 'cells:283-410', 'cells:411-538'], ['cells:19-146', 'cells:147-274', 'cells:275-402', 'cells:403-530']), ('boundary', ['boundary'], ['boundary'])]
print('%-9s %-22s %12s %12s %9s %9s  refined sub-cells' % ('orig cell', 'depth km', 'original', 'refined', 'share o', 'share r'))
for name, rk, gk in groups:
    a = sum(g[k]['charge'] for k in gk); b = sum(r[k]['charge'] for k in rk)
    dk = g[gk[0]].get('depth_km') or []
    print('%-9s %-22s %12.4e %12.4e %8.2f%% %8.2f%%  %s' % (name, ('%.1f-%.1f' % (dk[0], dk[1])) if dk else '-', a, b, 100*a/qg, 100*b/qr,
          ' '.join('%.3e' % r[k]['charge'] for k in rk) if len(rk) == 2 else ''))

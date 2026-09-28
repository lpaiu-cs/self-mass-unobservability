"""Independent native-input replay and exact sign extension of frozen bounds."""
from fractions import Fraction as F
import json, sys
import numpy as np
import exchange_curvature as e

g=e.g;OUT=e.OUT


def read(path): return json.loads(path.read_text())


def run():
    assert not (OUT/'audit.json').exists();plan=read(OUT/'plan.json')
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for name,digest in plan['trace_revision']['preserved'].items(): assert g.c.sha(OUT/name)==digest,name
    assert g.c.sha(e.LIB)==plan['trace_revision']['library_sha256']
    assert g.c.sha(e.c.s.LIB)==plan['library_sha256']
    assert not read(OUT/'before-cached-grouping-control.json')['passed']
    assert not read(OUT/'after-cached-grouping-control.json')['passed']
    assert read(OUT/'control.json')['passed'] and read(OUT/'result.json')['passed']
    eos=e.EOS();count=0;max_replay=0.;max_spectrum=0.;max_argument_gap=0.;changed_inputs=0;returned_mismatch=0.
    min_numerator=None;min_margin=1.;bitwise=True;max_Legendre=0.;values=[]
    for start in range(0,plan['cells'],plan['block_cells']):
        row=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert row['start']==count and row['stop']==min(count+128,plan['cells'])
        assert row['output_sha256']==g.c.sha(path) and row['plan_sha256']==g.c.sha(OUT/'plan.json')
        assert [r['cell'] for r in row['rows']]==list(range(count,row['stop']))
        a=dict(np.load(path));priorpath=e.c.OUT/f'block-{start}.npz';prior=dict(np.load(priorpath))
        assert g.c.sha(priorpath)==read(e.c.OUT/f'block-{start}.json')['output_sha256']
        for cell in plan['control_cells']:
            if start<=cell<row['stop']:
                control=dict(np.load(OUT/f'control-{cell}.npz'))
                for k in ['values','original','margins']: assert np.array_equal(a[k][cell-start],control[k])
        for j,b in enumerate(a['values']):
            old=a['original'][j];native=eos.probe_exchange(old[5],old[6])
            replay=abs(native[[0,1,2,4]]-old[:4])/np.maximum(abs(old[:4]),1e-100)
            max_replay=max(max_replay,float(replay.max()));bitwise &= np.array_equal(native[[0,1,2,4]],old[:4])
            returned=eos.probe_exchange(*a['arguments'][j]);assert np.array_equal(returned,b)
            assert np.array_equal(b[[4,5]],prior['values'][j,[10,12]])
            gap=float(old[5]-a['arguments'][j,0]);max_argument_gap=max(max_argument_gap,abs(gap));changed_inputs+=int(gap!=0)
            returned_mismatch=max(returned_mismatch,float((abs(b[[0,1,2,4]]-old[:4])/np.maximum(abs(old[:4]),1e-100)).max()))
            v=prior['values'][j]
            numerator=F(float(v[13]))+F(float(b[2]))
            assert numerator>0 and b[4]>0 and b[5]>0 and b[9]>0
            min_numerator=numerator if min_numerator is None else min(min_numerator,numerator)
            H=np.array([[v[4],v[5],v[6]],[v[5],v[7],v[8]],[v[6],v[8],v[9]]])/v[15]
            H[0,0]+=b[9]/b[4]
            spectrum=np.linalg.eigvals(prior['Gram'][j]@H);assert abs(spectrum.imag).max()<1e-10
            margin=1+min(0.,float(spectrum.real.min()));min_margin=min(min_margin,margin)
            max_spectrum=max(max_spectrum,abs(margin-a['margins'][j,1]))
            max_Legendre=max(max_Legendre,max(abs(b[13]),abs(b[14]))/max(1.,abs(b[1]),abs(b[12])))
            values.append(b)
        count=row['stop'];print('EXCHANGE INDEPENDENT AUDIT',count,'/',plan['cells'],flush=True)
    assert count==5735 and max_replay<plan['original_exchange_relative_tolerance'] and max_spectrum<1e-10
    assert max_Legendre<plan['normalized_Legendre_tolerance']
    result=read(OUT/'result.json');values=np.array(values)
    assert abs(min_margin-result['minimum_ideal_ions_Coulomb_electron_exchange_margin'])<1e-10
    assert [float(values[:,9].min()),float(values[:,9].max())]==result['electron_ne_H_over_kT_range']
    previous=read(e.c.OUT/'exact-matrix-bounds.json');bound=previous['minimum_lower_bound']
    assert F(bound)>0
    e.save('frozen-positive-extension.json',dict(classification='Proven',passed=True,cells=count,
        minimum_exact_electron_plus_exchange_numerator=str(min_numerator),
        inherited_component_lower_bound=bound,
        previous_certificate_sha256=g.c.sha(e.c.OUT/'exact-matrix-bounds.json'),
        proof='For every saved state, the exact rational numerator dpsi/dlnf + d(mu_exchange/kT)/dlnf is positive, and ne*dlnne/dlnf is positive. The combined ideal-electron/exchange Hessian is therefore a positive semidefinite rank-one correction on species space. Adding it cannot reduce the previously certified lower bound for ideal ions plus Coulomb on the same reconstructed positive support.',
        scope='Frozen coefficients and reconstructed populations only. Native continuous accuracy, returned versus last internal input differences, pressure ionization/excitation, missing populations and physical EOS errors are not certified.'))
    e.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,
        matched_cached_input_replay_bitwise=bool(bitwise),maximum_matched_input_relative_error=max_replay,
        returned_input_replay_bitwise=True,cells_with_returned_vs_cached_fl_difference=changed_inputs,
        maximum_abs_internal_fl_gap=max_argument_gap,maximum_original_returned_input_relative_mismatch=returned_mismatch,
        independent_spectrum_margin_error=max_spectrum,maximum_normalized_Legendre_error=max_Legendre,
        original_failed_controls_preserved=True,full_EOS_Hessian_certified=False,physical_EOS_certified=False))
    paths=[p for p in OUT.rglob('*') if p.is_file()]+[g.ROOT/'verification/exchange_curvature.py',g.ROOT/'verification/verify_exchange_curvature.py']
    e.save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        runtime={str(e.LIB):g.c.sha(e.LIB),str(e.CACHE/'exchange.so'):g.c.sha(e.CACHE/'exchange.so')},original_library_sha256=plan['library_sha256']))
    print('PASS EXCHANGE AUDIT',count,'matched input bitwise',bitwise,'exact component lower bound',bound,flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items(): assert g.c.sha(e.c.s.Path(path))==digest,path
    assert g.c.sha(e.c.s.LIB)==manifest['original_library_sha256']
    print('PASS EXCHANGE',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

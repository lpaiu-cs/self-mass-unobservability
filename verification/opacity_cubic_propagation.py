"""Carry the selected X-interpolator tangent through the actual Z/transport map."""
import json, sys
import numpy as np
import opacity_cubic as c

o=c.o;g=c.g;OUT=c.OUT/'propagation'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    prior=json.loads((c.OUT/'plan.json').read_text())
    assert g.c.sha(c.OUT/'before-context-native_opacity.py')==prior['inputs_sha256']['verification/native_opacity.py']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='55b5430',
        inputs_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            c.OUT/'result.json',o.OUT/'cubic-baseline-cubic.npz',o.OUT/'internal-baseline-captured.npz',
            o.OUT/'derivatives/finite-slopes.npz',g.ROOT/'verification/native_opacity.py',
            g.ROOT/'verification/opacity_cubic_propagation.py']},
        procedure='Read the actual Z grid and choices, verify original outputs bitwise, propagate the two captured X-table tangents through the original linear Z weights and harmonic radiative/conductive opacity. Preserve the original finite-derivative gate.',
        opacity_model_value_unchanged=True,continuous_native_or_physical_certificate=False))


def context():
    o.setup('context-baseline');o.trace('context-baseline',inspect_internal=True)
    base=dict(np.load(o.OUT/'internal-baseline-captured.npz'));actual=dict(np.load(o.OUT/'context-baseline-captured.npz'))
    equal={key:bool(np.array_equal(value,actual[key])) for key,value in base.items()};assert all(equal.values()),equal
    save('context-control.json',dict(classification='Counterexample candidate',passed=True,bitwise_equal=equal))


def run():
    assert json.loads((OUT/'context-control.json').read_text())['passed']
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    native=dict(np.load(o.OUT/'context-baseline-captured.npz'))
    controls=json.loads((o.OUT/'context-baseline-trace.json').read_text())['global_controls']
    assert controls['choices']==dict(cubic_interpolation_in_X=1,cubic_interpolation_in_Z=0,include_electron_conduction=1)
    assert np.all(native['inner'][:,3]>controls['kap_blend_logt_upper_bdy'])
    assert np.all(native['inner'][:,3]<controls['compton_blend_hi']-.5)
    grid=np.float32(controls['kap_z_tables']['Z']);assert np.all(np.diff(grid)>0)
    records=json.loads((c.OUT/'result.json').read_text())['rows'];bycell={}
    for row in records:
        assert row['axis']=='X';bycell.setdefault(row['cell'],[]).append(row)
    corrected=native['outputs'][:,1:].copy();rad_correction=np.zeros_like(corrected);reconstruction=[]
    factor=native['outputs'][:,0]/native['inner'][:,8];assert np.all((factor>0)&(factor<=1))
    for i,pair in bycell.items():
        assert len(pair)==2
        z=np.float32(native['parameters'][i,2]);j=int(np.searchsorted(grid,z,side='right'))-1
        assert 0<=j<len(grid)-1 and grid[j]<z<grid[j+1]
        beta=np.float32((z-grid[j+1])/(grid[j]-grid[j+1]));alpha=np.float32(1-beta)
        old=np.float32(beta*np.float32(pair[0]['original'])+alpha*np.float32(pair[1]['original'])).astype(float)
        reconstruction.append(float(abs(old-native['inner'][i,9:11]).max()))
        new=float(beta)*np.array(pair[0]['corrected'])+float(alpha)*np.array(pair[1]['corrected'])
        rad_correction[i]=new-native['inner'][i,9:11]
        corrected[i]+=factor[i]*rad_correction[i]
    assert max(reconstruction)<1e-12,max(reconstruction)
    finite=dict(np.load(o.OUT/'derivatives/finite-slopes.npz'))['slopes'];rows=[]
    for j,axis in enumerate(['lnd','lnT']):
        for k,h in enumerate([1e-3,5e-4]):
            score=abs(finite[j,k]-corrected[:,j])/np.maximum(1,abs(corrected[:,j]));worst=int(np.argmax(score))
            rows.append(dict(axis=axis,step=h,maximum_score=float(score[worst]),worst_cell=worst,
                failed_cells=int(np.count_nonzero(score>1e-3)),passed=bool(np.all(score<=1e-3))))
    np.savez_compressed(OUT/'corrected-tangents.npz',corrected=corrected,rad_correction=rad_correction)
    save('result.json',dict(classification='Counterexample candidate',completed=True,corrected_cubic_cells=len(bycell),
        native_radiative_derivative_reconstruction_maximum=max(reconstruction),
        maximum_total_log_derivative_changes=abs(corrected-native['outputs'][:,1:]).max(0).tolist(),
        finite_rows=rows,all_original_finite_gates_passed=all(row['passed'] for row in rows),
        physical_opacity_certified=False,continuous_full_derivatives_certified=False))
    print('PROPAGATED CUBIC TANGENTS',len(bycell),max(reconstruction),rows,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

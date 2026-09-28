"""Same-material GR reconstruction with the two-spectrum EOS and strict roots."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from types import FunctionType, MethodType, SimpleNamespace
import json, sys
import numpy as np
import gr_molecular_precision as precision
import gr_molecular_reference_runner as reference
import gr_subcell_precision_fallback as fallback

g=precision.g;model=precision.model;OUT=g.OUT/'gr-molecular-structure-defined'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();reference.verify()
    paths=[g.ROOT/'verification/gr_molecular_structure.py',g.ROOT/'verification/common_eos.py',
        g.ROOT/'verification/baryon_entropy.py',g.ROOT/'verification/audit_structured_enthalpy.py',
        g.ROOT/'verification/gr_subcell_precision_fallback.py',precision.OUT/'plan.json',
        precision.OUT/'build.json',reference.OUT/'manifest.json',g.OUT/'initial-adiabats-17.npz',
        g.OUT/'initial-structure-17-4.json',g.OUT/'initial-state-17-4.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='2c76707',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cells=5735,block_size=128,workers=4,points=17,subdivisions=[4,8],
        controls=[0,1175,1176,2972,5734],table_control_cells=[0,2972,5734],
        interface_tolerance=1e-8,finite_face_refinement_tolerance=1e-8,
        known_root_log_temperature_tolerance=1e-10,
        table='New EOS pressure and entropy at the frozen rho,T,X reference. Seventeen offsets in [-0.04,0.04] about each new pressure. Old table temperature differences about its midpoint are only Newton seeds added to the new reference temperature; every new entry is solved with the new EOS.',
        inverse='Reuse the frozen strict inverse and precision fallback without changing either energy-unit gate. The fallback uses the matching molecular-precision EOS, never the old physical model. Isolate and retain all failures and switches per block.',
        startup_gate='The molecular extended evaluator controls must complete and all pass. Freeze their manifest at table start.',
        GR='Reuse existing common_eos baryon TOV shooting, 17-point Pchip and direct strict-root fallback. Hold dm,X and the new reference entropy; the same finite photospheric pressure is the boundary. Centre pressure, R and gravitational M are matched outputs, with no baryon tuning. Compute subdivision 4 then 8.',
        scope='Numerical reconstruction of a new model, not physical time evolution. No physical EOS, continuous interpolation/root/geometry, atmosphere, transport or nonlinear observational certificate. Original models and failed paths remain unchanged.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def inverse(stats,folder):
    folder.mkdir(exist_ok=True)
    local_save=FunctionType(save.__code__,dict(save.__globals__,OUT=folder))
    fn=fallback.make_inverse
    return FunctionType(fn.__code__,dict(fn.__globals__,OUT=folder,save=local_save,
        quad=SimpleNamespace(EOS=precision.EOS,solve=precision.q.solve)))(stats)


def inputs():
    data=dict(np.load(reference.OUT/'reference-state.npz'))
    refs=np.concatenate([np.load(reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,len(data['X']),128)])
    old=dict(np.load(g.OUT/'initial-adiabats-17.npz'))
    assert len(refs)==len(data['X']) and np.array_equal(refs[:,3],data['s_B'])
    assert np.array_equal(old['X'],data['X'])
    return data,refs,old


def seed(data,old,i,j):return float(data['lnT'][i]+old['values'][i,j,1]-old['values'][i,8,1])


def block(start):
    plan=bindings();folder=OUT/f'block-{start:04}';folder.mkdir()
    data,refs,old=inputs();stop=min(start+plan['block_size'],plan['cells']);eos=model.EOS()
    stats=dict(calls=0,evaluations=0,maximum_score=0.,label=folder.name);inv=inverse(stats,folder)
    values=[];offset=np.linspace(-.04,.04,plan['points'])
    for i in range(start,stop):
        row=[]
        for j,dx in enumerate(offset):
            a,t,_=inv(eos,float(np.log(refs[i,1])+dx),refs[i,3],data['X'][i],seed(data,old,i,j))
            row.append([np.log(a[0]),t,a[2]*1e-4/g.c.gr.C**2])
        values.append(row)
    path=folder/'values.npz';np.savez_compressed(path,cells=np.arange(start,stop),values=np.array(values))
    record=dict(classification='Counterexample candidate',start=start,stop=stop,
        plan_sha256=g.c.sha(OUT/'plan.json'),output_sha256=g.c.sha(path),root_statistics=stats)
    (folder/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    print('MOLECULAR ADIABATS',start,stop,stats,flush=True);return record


def table():
    plan=bindings();precision.verify();assert json.loads((precision.OUT/'result.json').read_text())['all_passed']
    assert not any(OUT.glob('block-*'));save('startup-binding.json',dict(sha256=g.c.sha(precision.OUT/'manifest.json')))
    data,refs,old=inputs();eos=model.EOS();stats=dict(calls=0,evaluations=0,maximum_score=0.,label='known-roots')
    inv=inverse(stats,OUT/'known-roots');controls=[]
    for i in plan['controls']:
        lp=float(np.log(refs[i,1]));t=float(data['lnT'][i]);a=eos(1,lp,t,data['X'][i])
        _,actual,_=inv(eos,lp,a[3],data['X'][i],t+.01)
        controls.append(dict(cell=i,lnT_error=abs(actual-t),passed=bool(abs(actual-t)<plan['known_root_log_temperature_tolerance'])))
    save('known-roots.json',dict(classification='Counterexample candidate',rows=controls,statistics=stats))
    assert all(r['passed'] for r in controls)
    starts=list(range(0,plan['cells'],plan['block_size']));records=[];errors=[]
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        futures={pool.submit(block,i):i for i in starts}
        for done in as_completed(futures):
            try:records.append(done.result())
            except Exception as error:errors.append(dict(start=futures[done],error=repr(error)))
            save('table-progress.json',dict(classification='Counterexample candidate',completed_cells=sum(r['stop']-r['start'] for r in records),errors=errors))
    assert not errors,errors
    values=np.concatenate([np.load(OUT/f'block-{i:04}'/'values.npz')['values'] for i in starts])
    offset=np.linspace(-.04,.04,plan['points']);stats=dict(calls=0,evaluations=0,maximum_score=0.,label='table-controls')
    inv=inverse(stats,OUT/'table-controls');serial=[]
    for i in plan['table_control_cells']:
        for j in [0,8,16]:
            a,t,_=inv(eos,float(np.log(refs[i,1])+offset[j]),refs[i,3],data['X'][i],seed(data,old,i,j))
            answer=np.array([np.log(a[0]),t,a[2]*1e-4/g.c.gr.C**2])
            assert np.array_equal(answer,values[i,j]);serial.append([i,j])
    np.savez_compressed(OUT/'molecular-adiabats-17.npz',reference=refs,values=values,offset=offset,
        X=data['X'],lnT=data['lnT'],lnd=data['lnd'])
    save('table-result.json',dict(classification='Counterexample candidate',completed=True,cells=len(values),
        roots=int(np.prod(values.shape[:2])),serial_bitwise_controls=serial,
        output_sha256=g.c.sha(OUT/'molecular-adiabats-17.npz'),physical_EOS_certified=False))


class Structure(g.c.Structure):
    def __init__(self,label,data,points=17,subdivision=4):
        fn=g.c.Structure.__init__
        FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,EOS=model.EOS))(self,label,data,points,subdivision)
        self.stats=dict(calls=0,evaluations=0,maximum_score=0.,label=f'GR-{subdivision}')
        self.inverse=inverse(self.stats,OUT/f'GR-{subdivision}-roots')
        fn=g.c.Structure.state
        state=FunctionType(fn.__code__,dict(fn.__globals__,be=SimpleNamespace(invert=self.inverse)))
        self.state=MethodType(state,self)


def project():
    plan=bindings();t=json.loads((OUT/'table-result.json').read_text())
    assert t['completed'] and t['output_sha256']==g.c.sha(OUT/'molecular-adiabats-17.npz')
    data,_,_=inputs();old=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    parameters=json.loads((g.OUT/'initial-structure-17-4.json').read_text())['parameters'];reports=[];prior=None
    for sub in plan['subdivisions']:
        assert not (OUT/f'molecular-state-17-{sub}.npz').exists()
        holder=[]
        def create(*args):
            obj=Structure(*args);holder.append(obj);return obj
        def invert(eos,lp,entropy,X,guess):return holder[0].inverse(eos,lp,entropy,X,guess)
        fn=g.c.structure
        run=FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save,Structure=create,be=SimpleNamespace(invert=invert)),argdefs=fn.__defaults__)
        state,report=run('molecular',data,17,sub,parameters);parameters=report['parameters']
        assert np.array_equal(state['dm'],data['dm']) and np.array_equal(state['X'],data['X'])
        report.update(root_statistics=holder[0].stats,
            model_replacement_radius_change_m=float(state['radius_faces_m'][0]-old['radius_faces_m'][0]),
            model_replacement_mass_change_geom_m=float(state['mass_faces_geom'][0]-old['mass_faces_geom'][0]),
            maximum_reference_entropy_difference=float(abs(state['s_B']-data['s_B']).max()),
            physical_time_evolution=False)
        if prior is not None:
            gap=max(float(abs(state['radius_faces_m']-prior['radius_faces_m']).max()/state['radius_faces_m'][0]),
                float(abs(state['mass_faces_geom']-prior['mass_faces_geom']).max()/state['mass_faces_geom'][0]),
                float(abs(state['logP']-prior['logP']).max()))
            report.update(finite_face_refinement=gap,finite_face_refinement_passed=gap<plan['finite_face_refinement_tolerance'])
        prior=state;reports.append(report);save(f'GR-{sub}.json',report)
    save('result.json',dict(classification='Counterexample candidate',completed=True,reports=reports,
        all_interfaces_passed=all(r['interface_max']<plan['interface_tolerance'] for r in reports),
        finite_face_refinement_passed=reports[-1]['finite_face_refinement_passed'],
        physical_EOS_certified=False,continuous_errors_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.rglob('*') if p.is_file()}))
    verify()


def verify():
    plan=bindings();precision.verify()
    assert g.c.sha(precision.OUT/'manifest.json')==json.loads((OUT/'startup-binding.json').read_text())['sha256']
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    table=json.loads((OUT/'table-result.json').read_text());assert table['completed'] and table['cells']==plan['cells']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS molecular GR reconstruction bindings; consult finite gates, not time evolution',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

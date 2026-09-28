"""Test the existing molecular continuation against the failed material increment."""
import ast,inspect,json,sys
import numpy as np
import gr_material_thermo_increment as original
import eos_molecular_switch as switch

g=original.g;OUT=g.OUT/'gr-material-thermo-continuation'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


class Continuation:
    def __init__(self):self.eos=switch.EOS()
    def __call__(self,mode,r,t,X):
        assert mode==2
        return self.eos.sample(r,t,X,True)[0]['eos']


def overlay():
    source=inspect.getsource(original.run);tree=ast.parse(source)
    calls=[node for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
        and node.func.id=='save' and isinstance(node.args[0],ast.Constant) and node.args[0].value=='plan.json']
    assert len(calls)==1
    replacements={'assert not OUT.exists();OUT.mkdir()':'assert OUT.exists()',
        ast.get_source_segment(source,calls[0]):"plan=json.loads((OUT/'plan.json').read_text())",
        'eos=g.EOS()':'eos=Continuation()'}
    changed=source
    for old,new in replacements.items():assert changed.count(old)==1;changed=changed.replace(old,new)
    reverse=changed
    for old,new in reversed(list(replacements.items())):reverse=reverse.replace(new,old,1)
    assert reverse==source
    return source,changed,replacements


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify();switch.verify()
    failed=json.loads((original.OUT/'result.json').read_text())['failed_cases']
    assert len(failed)==1 and failed[0]['kind']=='known temperature cut'
    source,changed,replacements=overlay()
    (OUT/'original-run.py').write_text(source);(OUT/'continuation-run.py').write_text(changed)
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='d8689de',substitutions=replacements,
        intervention='Use the already built, previously independently audited molecular-switch library with retain_molecules=True. No rebuild or source/runtime replacement. Keep every path, quadrature, endpoint comparison and original numerical gate.',
        physical_scope='This removes the explicit temperature cutoff only inside the existing Taylor partition continuation comparison. That continuation has not been calibrated as a physical molecular partition function. Numerical continuity/increment recovery does not certify physical EOS or replace the old GR/time paths.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_material_thermo_continuation.py',original.OUT/'manifest.json',
        switch.OUT/'manifest.json',OUT/'continuation-run.py']})
    plan['runtime']=switch.read('manifest.json')['runtime'];save('plan.json',plan)
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));i=failed[0]['cell'];eos=switch.EOS();baseline=g.EOS()
    lr=state['lnd'][i];lt=np.log(1e6);X=state['X'][i];rows=[];values=[]
    for h in [0.,1e-8,1e-10,1e-12]:
        snapshots=[]
        for t in [lt-h,lt+h]:
            off,_=eos.sample(lr,t,X,False);on,_=eos.sample(lr,t,X,True)
            assert np.array_equal(off['eos'],baseline(2,lr,t,X))
            snapshots.extend([off['eos'],on['eos']])
        a=np.array(snapshots);difference=a[2].astype(np.longdouble)-a[3].astype(np.longdouble)
        rows.append(dict(cell=int(i),log_halfwidth=h,high_side_pressure_model_difference=float(difference[1]),
            high_side_energy_model_difference=float(difference[2]),high_side_entropy_model_difference=float(difference[3]),
            high_side_free_energy_model_difference=float(difference[2]-np.exp(np.longdouble(lt+h))*difference[3]),
            default_false_bitwise_original=True))
        values.append(a)
    np.savez_compressed(OUT/'switch-isolation.npz',values=np.array(values),log_halfwidths=np.array([0.,1e-8,1e-10,1e-12]))
    save('switch-isolation.json',dict(classification='Counterexample candidate',rows=rows,
        scope='Same-input switch intervention and finite one-sided probes. Original and retained outputs are distinct models; these differences are not actual latent heat or an observed phase transition.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['runtime'].items():assert g.c.sha(path)==digest,path
    _,changed,_=overlay();assert changed==(OUT/'continuation-run.py').read_text()
    namespace=dict(original.run.__globals__,OUT=OUT,save=save,verify=verify,Continuation=Continuation)
    exec(compile(changed,str(OUT/'continuation-run.py'),'exec'),namespace);namespace['run']()


def verify():
    from types import FunctionType
    FunctionType(original.verify.__code__,dict(original.verify.__globals__,OUT=OUT))()
    for path,digest in json.loads((OUT/'plan.json').read_text())['runtime'].items():assert g.c.sha(path)==digest,path


if __name__=='__main__':globals()[sys.argv[1]]()

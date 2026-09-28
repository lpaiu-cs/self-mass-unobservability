"""Split the exactly affine lapse-gradient force before material differencing."""
from pathlib import Path
from types import FunctionType
import json,resource,sys,time
import numpy as np
import solve_native_incident_lift as lift

base=lift.base;OUT=lift.OUT;MATERIAL=OUT/'material-affine';read,write,sha=lift.read,lift.write,lift.sha
LD=np.longdouble;old_initialize=lift.initialize


def initialize():
    old_initialize();Parent=base.Material
    class Material(Parent):
        split_ap=False
        def __init__(self,reference=128,steps=128):
            self.precise=(OUT/'material-precision-plan.json').exists()
            super().__init__(reference,steps);self.ap_maps={};self.ap_errors=[]
            if self.precise:self.model.join0=np.asarray(self.model.join0,LD)
        def raw(self,k,delta,field,eps,row=None):
            if self.precise:
                row=dict(self.point(k) if row is None else row)
                for key in ['Q','h','theta','seed']:row[key]=np.asarray(row[key],LD)
                # np.interp consumes the saved binary64 field samples; the
                # probe multiplication and geometry arithmetic remain extended.
                delta,field,eps=np.asarray(delta,LD),np.asarray(field,float),LD(eps)
            return super().raw(k,delta,field,eps,row)
        def fields(self,t):
            j,f,field,rates=super().fields(t)
            if self.split_ap:field=field.copy();field[4]=0.
            return j,f,field,rates
        def ap_map(self,k):
            if k in self.ap_maps:return self.ap_maps[k]
            row=self.point(k);delta=np.zeros_like(row['Q']);field=np.zeros((5,self.n))
            field[4]=np.maximum(abs(self.ap),self.a/self.R);values=[]
            for h in [1.,.5]:
                a=self.raw(k,delta,field,h);b=self.raw(k,delta,field,-h)
                assert np.array_equal(a[0],b[0]),'Acceleration entered a material flux'
                values.append((a[1].astype(LD)-b[1].astype(LD))/(2*LD(h)*field[4]))
            error=float(np.max(abs(values[1]-values[0]))/max(np.max(abs(values[1])),1e-290))
            assert error<1e-12,error
            self.ap_errors.append(error);self.ap_maps[k]=np.asarray(values[1],float)
            return self.ap_maps[k]
        def rhs(self,t,z,probe=1.):
            j,f,field,_=super().fields(t)
            force=field[4]*((1-f)*self.ap_map(j)+f*self.ap_map(j+1))
            self.split_ap=True
            try:r,l,dt=super().rhs(t,z,probe)
            finally:self.split_ap=False
            r[1]+=force;l[1]+=np.sum(force,dtype=LD)
            return r,l,dt
        run=FunctionType(Parent.run.__code__,dict(Parent.run.__globals__,OUT=MATERIAL),argdefs=Parent.run.__defaults__)
    base.Material=Material;base.MATERIAL=MATERIAL;base.PHOTON=lift.PHOTON


def prepare():
    assert not read(OUT/'material/pilot-64.json')['passed'];assert not MATERIAL.exists();MATERIAL.mkdir()
    write(OUT/'material-affine-plan.json',dict(classification='Counterexample candidate',
        failure=read(OUT/'material/pilot-64.json')['directional_relative'],
        repair='At fixed primitive/geometry state, the lapse-gradient acceleration ap enters only the momentum source affinely. Measure its diagonal coefficient once at each saved background knot, verify opposite directions and two probe sizes, and remove that channel before the existing native finite directional response. Restore its exact linear force in both rate and momentum ledger.',
        reason='A rapidly varying incident wave can make ap large relative to the hydrostatic acceleration while all primitive and dimensionless geometry changes stay small. Its old joint probe cap then suppresses the recoverable thermodynamic increment. The physical amplitude is unchanged.',
        controls='Keep8x Richardson for the remaining channels and original4/8/16 comparisons,0.2percent derivative gate,1e-8 conservation and owner gates,1percent actual donor gate and2percent time gate. No changed waveform,clock or physical state.',
        scope='Exact affine separation in the current finite material equation, not a full native EOS Jacobian certificate. The ap coefficient itself must agree below1e-12.',
        budgets=dict(check=45,pilot=60,production=400),total_production_allocation=1500,
        bindings={str(p):sha(p) for p in [Path(__file__),Path(lift.__file__),Path(base.__file__),OUT/'material/pilot-64.json',lift.PHOTON/'result.json']}))


def check():
    initialize();m=base.Material(128,64);d=np.load(OUT/'material/pilot-64.npz');rows=[]
    for i in [1,2]:
        t=d['t'][i];z=d['history_scaled'][i];rates=[m.rhs(t,z,p)[0] for p in [.5,1.,2.]]
        norm=np.maximum(np.sum(abs(rates[1]),axis=1),1.)
        errors=[(np.sum(abs(r-rates[1]),axis=1)/norm).tolist() for r in [rates[0],rates[2]]]
        rows.append(dict(time=float(t),probe_4_8_16=errors,passed=bool(np.max(errors)<.002)))
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,
        affine_map_relative=max(m.ap_errors),physical_branch_ratio=m.physical_branch_ratio)
    write(OUT/('material-precision-check.json' if m.precise else 'material-affine-check.json'),result)
    print(json.dumps(result),flush=True);assert result['passed'],result


check_json=check
check_extended=check
check_interpolation=check


def precision_plan():
    assert not read(OUT/'material-affine-check.json')['passed']
    write(OUT/'material-precision-plan.json',dict(classification='Counterexample candidate',
        evidence='Saved early-state probes do not share a converged original binary64 arithmetic range. The largest momentum discrepancy is in deep cells0..5; changing only the affine acceleration did not solve it. Preserve both failures.',
        repair='Carry Q,h,temperature seed,field,probe and deep shared-face state in longdouble through the existing recovery and Horner EOS evaluation, reusing the existing extended flux algebra. Do not change EOS coefficients, physical amplitude, state or tolerance.',
        gates='Original4/8/16 Richardson comparisons below0.2percent, owner and conservation1e-8, branch1percent; affine map1e-12. Actual prefixes must pass again before any production.',
        caps=dict(check=35,pilot=60,production=400),source_sha256=sha(__file__),
        bindings={str(p):sha(p) for p in [OUT/'material-affine-check.json',OUT/'material-affine-diagnostic.json',OUT/'diagnostic-material-affine-producer.py']}))


def diagnose():
    initialize();m=base.Material(128,64);d=np.load(OUT/'material/pilot-64.npz');rows=[];saved={}
    for i in [1,2]:
        t=d['t'][i];z=d['history_scaled'][i];rates={}
        for p in [1.,2.,4.,8.,16.,32.]:
            rates[p]=m.rhs(t,z,p)[0];saved[f'rate-{i}-{p}']=rates[p]
        errors={str(p):(np.sum(abs(rates[p]-rates[p/2]),axis=1)/np.maximum(np.sum(abs(rates[p]),axis=1),1.)).tolist() for p in [2.,4.,8.,16.,32.]}
        diff=abs(rates[2]-rates[1]);ids=np.argsort(diff[1])[-6:][::-1]
        rows.append(dict(time=float(t),comparisons=errors,largest_momentum_cells=ids.tolist(),
                         difference=diff[1,ids].tolist(),reference_momentum=rates[2][1,ids].tolist()))
    np.savez_compressed(OUT/'material-affine-diagnostic.npz',**saved)
    write(OUT/'material-affine-diagnostic.json',dict(classification='Counterexample candidate',rows=rows,
        physical_branch_ratio=m.physical_branch_ratio,trajectory_replayed=False,
        planned_scope='Two saved early states,6 arithmetic probe sizes only; physical eta remains1e-30. No acceptance by this diagnostic alone.'))


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','check','check_json','diagnose','precision_plan','check_extended','check_interpolation','pilot','production']
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));base.drive.native.deadline(dict(prepare=30,check=45,check_json=39,diagnose=30,precision_plan=30,check_extended=35,check_interpolation=30,pilot=60,production=400)[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action in ['pilot','production']:
            assert read(OUT/'material-precision-check.json')['passed'];base.initialize=initialize;base.PHOTON=lift.PHOTON;base.MATERIAL=MATERIAL
            if action=='production':assert read(MATERIAL/'pilot.json')['upper_remaining_seconds']<400
            base.material(action=='pilot')
            if action=='pilot':assert read(MATERIAL/'pilot.json')['upper_remaining_seconds']<400
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        p=OUT/f'material-affine-{action}-receipt.json';assert not p.exists()
        write(p,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))

"""Differentiate the actual donor branch at its true, sub-ulp amplitude.

Counterexample candidate: amplified probes must not silently switch a nonzero
background mass flux. At exactly zero flux retain the directional donor.
"""
from pathlib import Path
from types import FunctionType,MethodType
import inspect,textwrap,json,signal,sys,time
import numpy as np
import sympy as sp
import def_native_material_response as base

OUT=base.OUT/'branch-aware';write=base.write;sha=base.sha;AMP=base.AMP


class Material(base.Material):
    def __init__(self,reference):
        super().__init__(reference);self.reference_flux=None;self.branch_crossings=0;self.physical_branch_ratio=0.
        m=self.model;f=m.flow
        def choose(flux,section):
            if self.reference_flux is None:return flux>=0
            saved=self.reference_flux[0,section]
            self.branch_crossings+=int(np.count_nonzero((saved!=0)&((saved>=0)!=(flux>=0))))
            return np.where(saved==0,flux>=0,saved>=0)
        m.material_donor=lambda flux:choose(flux,slice(0,self.nb+1))
        f.material_donor=lambda flux:choose(flux,slice(self.nb,None))
        m.join_donor=lambda flux:bool(choose(np.array([flux]),slice(self.nb,self.nb+1))[0])
        atmosphere,deep=self.expanded
        atmosphere=base.replace(atmosphere,'flux[0]>=0','self.material_donor(flux[0])')
        assert deep.count('mass>=0')==2;deep=deep.replace('mass>=0','self.material_donor(mass)')
        ns=dict(base.flow.hydro_ns);exec(compile(atmosphere,__file__,'exec'),ns);self.atmosphere=MethodType(ns['rhs'],f)
        ns=dict(m.material_rhs.__func__.__globals__);exec(compile(deep,__file__,'exec'),ns);self.deep=MethodType(ns['material_rhs'],m)
        join=textwrap.dedent(inspect.getsource(m.join_flux));join=base.replace(join,'flux[0]>=0','self.join_donor(flux[0])')
        ns=dict(m.join_flux.__func__.__globals__);exec(compile(join,__file__,'exec'),ns);m.join_flux=MethodType(ns['join_flux'],m);f.join_flux=m.join_flux
        self.expanded=(atmosphere,deep);self.join_source=join

    def raw(self,k,delta,field,eps,row=None):
        row=self.point(k) if row is None else row;self.reference_flux=row.get('flux')
        result=super().raw(k,delta,field,eps,row)
        if self.reference_flux is not None and eps:
            saved=self.reference_flux[0];active=saved!=0
            ratio=AMP*np.max(abs((result[0][0,active]-saved[active])/eps/saved[active]))
            self.physical_branch_ratio=max(self.physical_branch_ratio,float(ratio));assert ratio<.01,('True amplitude leaves donor branch',ratio)
        return result

    run=FunctionType(base.Material.run.__code__,dict(vars(base),OUT=OUT),argdefs=base.Material.run.__defaults__)


check=FunctionType(base.check.__code__,dict(vars(base),OUT=OUT,Material=Material))
production=FunctionType(base.production.__code__,dict(vars(base),OUT=OUT,Material=Material))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Resolve the failed material directional derivative, then complete the same three material response paths; no changed physical flux or acceptance threshold.',
        root_cause='Amplified finite probes cross tiny nonzero deep mass-flux signs, changing donor thermochemical cells. The physical1e-26 response stays on the original branch. A step-size limit cannot be represented by unguarded finite probes across that branch.',
        repair='Use the saved nonzero mass-flux donor; at exactly zero flux choose the directional donor. Apply consistently to deep energy/H, atmosphere H and the shared H face. Require actual-amplitude flux variation/background ratio<.01. Keep original HLL flux, native primitive recovery and all other branches.',
        costs='First failed complete64 path cost40.82s including setup. Preserve it. Corrected three-path cap250s after measured prefixes and2x margin; additional check35s,pilot35s. No native bank, new horizon, spatial grid or clock path.',
        gates=dict(owner=1e-8,directional=.002,conservation=1e-8,time=.02,background=.02,small_relative=1e-6,physical_branch_ratio=.01),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),base.OUT/'production.json',base.OUT/'probe-diagnostic.json',base.OUT/'steps-64-reference-128.npz']}))
    F,dF,y,dy,h=sp.symbols('F dF y dy h');assert sp.expand((F+h*dF)*(y+h*dy)-F*y-h*(dF*y+F*dy))==h*h*dF*dy
    write(OUT/'branch-symbolic.json',dict(classification='Proven',passed=True,
        scope='For nonzero F, the donor branch is unchanged whenever abs(h*dF)<abs(F). On that branch delta(F*y)=deltaF*y+F*delta_y. At F=0 the one-sided directional donor is set by deltaF. No assertion about all limiter branches or a uniform EOS derivative bound.'))


def pilot():
    assert json.loads((OUT/'check.json').read_text())['passed'];assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,base.flow.old.optical.timeout);signal.alarm(35);rows=[]
    for steps,reference in [[64,128],[128,128],[128,64]]:
        m=Material(reference);r=m.run(steps,f'pilot-{steps}-{reference}',2);r.update(branch_crossings=m.branch_crossings,physical_branch_ratio=m.physical_branch_ratio);rows.append(r)
        if not r['passed']:break
    forecast=sum(r['seconds']*(r['steps']-2)/2 for r in rows);upper=2*forecast+10
    result=dict(classification='Counterexample candidate',paths=rows,forecast_seconds=forecast,upper_seconds=upper,eligible=len(rows)==3 and all(r['passed'] for r in rows) and upper<250,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);print(json.dumps(result));signal.alarm(0)
    if result['eligible']:
        write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,production_seconds=250,forecast_seconds=forecast,upper_seconds=upper,paths=[[64,128],[128,128],[128,64]],
            binds={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),OUT/'pilot.json',OUT/'check.json',OUT/'plan.json']}))


if __name__=='__main__':globals()[sys.argv[1]]()

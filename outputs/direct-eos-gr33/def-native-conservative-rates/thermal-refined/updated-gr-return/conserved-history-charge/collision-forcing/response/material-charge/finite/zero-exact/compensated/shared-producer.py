"""Subtract the fixed deep pressure flux before forming the finite defect."""
from pathlib import Path
from types import FunctionType,MethodType,SimpleNamespace
import inspect,json,signal,sys,textwrap,time
import numpy as np
import sympy as sp
import def_native_finite_collision_return as old

OUT=old.OUT/'compensated';GR=OUT/'gr';C=old.C;AMP=old.AMP
write=old.write;sha=old.sha;configure=old.configure


def replace(s,a,b):assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)


class Material(old.Material):
    def __init__(self,reference,steps=128):
        super().__init__(reference,steps);m=self.model
        common='momentum=self.area_gas*(self.face_p+ps+self.cx*self.face_rho*vf*vf)'
        centered='momentum=self.area_gas*(ps+self.cx*self.face_rho*vf*vf)'
        join='mass[-1],momentum[-1],energy[-1],neutral[-1]=self.mflux'
        subtract=join+'\n    momentum[-1]-=self.base_momentum_flux[-1]'
        # Change the common owner and its flux-export version together.
        physical=old.old.physical.branch.base.centered.method
        physical=replace(physical,common,centered);physical=replace(physical,join,subtract)
        physical=replace(physical,'np.diff(momentum-self.base_momentum_flux)','np.diff(momentum)')
        scope=dict(m.material_rhs.__func__.__globals__);exec(compile(physical,__file__,'exec'),scope)
        m.material_rhs=MethodType(scope['material_rhs'],m)
        deep=replace(self.expanded[1],common,centered);deep=replace(deep,join,subtract)
        deep=replace(deep,'np.diff(self.base_momentum_flux)+b.volume*self.f0["initial_support"]+geometry','b.volume*self.f0["initial_support"]+geometry')
        scope=dict(self.deep.__func__.__globals__);exec(compile(deep,__file__,'exec'),scope);self.deep=MethodType(scope['material_rhs'],m)
        self.centered_sources=(physical,deep)

    raw_source=inspect.getsource(old.old.physical.branch.base.Material.raw)
    raw_source=textwrap.dedent(raw_source)
    raw_source=replace(raw_source,'    shared=float(',
        '    original_shared=abs((deep[1,-1]+C*m.base_momentum_flux[-1])-aflux[1,0])/max(abs(aflux[1,0]),1.)\n    assert original_shared<1e-12\n    aflux[1,0]=deep[1,-1]\n    ag[0]+=C*m.base_momentum_flux[-1]\n    shared=float(')
    scope=dict(old.old.physical.branch.base.__dict__);exec(compile(raw_source,__file__,'exec'),scope)
    raw=scope['raw']
    run=FunctionType(old.Material.run.__code__,dict(old.Material.run.__globals__,OUT=OUT),argdefs=old.Material.run.__defaults__)


def prepare():
    assert not OUT.exists();OUT.mkdir();GR.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='3ca3359c5',
        claim='Remove the fixed deep pressure flux before finite subtraction, retain the same material equation, then finish the actual native collision return and GR charge.',
        algebra='For fixed face vector B, use Fnew=F-B and Gnew=G-div(B), giving exactly -div(Fnew)+Gnew=-div(F)+G. Deep B is the original area*reference pressure; exterior faces have B=0, and its one shared-face contribution is included identically on both sides.',
        limits='Finite actual-amplitude material return on the same17 stored backgrounds and accepted photon response. No new microscopic EOS inversions, full nonlinear radiation feedback, uniform EOS error bound, or final physical closure is claimed.',
        budgets=dict(check_seconds=30,pilot_seconds=60,production_seconds=780,source_seconds=90,GR_seconds=90,CPU_threads=1,virtual_GiB=3),
        resource_reassessment='Reuse the unstarted780s material production and90s source/GR caps. Permit at most30s owner/symbolic checks and60s measured two-path pilot for this arithmetic repair. No additional clocks, cells, support or horizon. Original failures remain frozen.',
        forecast='Reuse prior complete late-CFL raw-call counts and measured finite prefix cost. Require twice total forecast under780s and resume the new accepted prefixes.',
        gates=dict(owner=1e-8,balance=1e-8,finite_resolution=.002,time=.02,pressure_resolution=.002,quadrature=.002,independent_GR=1e-9),
        stop='Stop on an original numerical gate or budget. No amplitude reduction, fine-clock expansion, or threshold change.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),old.OUT/'pilot-64.json',old.old.photons.OUT/'result.json',old.OUT/'pilot-64.npz']}))
    a,b,c,d=sp.symbols('Fleft Fright Bleft Bright');g=sp.symbols('G')
    assert sp.expand(-((b-d)-(a-c))+g-(d-c)-(-(b-a)+g))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='A fixed face-flux subtraction and its matching cell-source subtraction leave every semidiscrete conservation rate unchanged; the shared face still cancels. No EOS or temporal error theorem.'))


def check():
    start=time.monotonic();signal.alarm(30);configure();m=Material(128,64);before=old.Material(128,64)
    saved=dict(np.load(old.OUT/'pilot-64.npz'));rows=[]
    for k in [0,8,16]:
        a=m.point(k);b=before.point(k);z=saved['delta_scaled']
        now=m.raw(k,z,np.zeros((5,m.n)),AMP);prior=before.raw(k,z,np.zeros((5,m.n)),AMP)
        rn=-np.diff(now[0],axis=1);rn[1]+=now[1]
        ro=-np.diff(prior[0],axis=1);ro[1]+=prior[1]
        err=(np.sum(abs(rn-ro),axis=1)/np.maximum(np.sum(abs(ro),axis=1),1.)).tolist()
        assert max(err)<1e-8 and a['owner_error']<1e-8
        rows.append(dict(k=k,rate_equivalence=err,owner_error=a['owner_error']))
    for name,source in zip(['physical','export'],m.centered_sources):(OUT/f'expanded-{name}.py').write_text(source)
    (OUT/'expanded-raw.py').write_text(Material.raw_source)
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start));signal.alarm(0)


def repair_prepare():
    write(OUT/'first-check-failure.json',dict(classification='Counterexample candidate',error='Shared compensated face mismatch3.800172907318631e-5',
        physical_steps=0,charged_check_seconds=30,reason='The same physical face was converted and its common pressure subtracted in two arithmetic orders; cancellation exposed their rounding difference.'))
    p=json.loads((OUT/'plan.json').read_text());p['shared_face_repair']='Check the original full face at1e-12, then use ONE existing deep-unit centered face value identically for both cells. Do not independently subtract the same pressure twice.'
    p['budgets']['total_check_seconds']=60
    p['resource_reassessment']='One additional30s owner check for the shared arithmetic repair; original check charged30s, no physical time steps repeated. All pilot/production/source/GR caps unchanged.'
    p['bindings'].update({str(q):sha(q) for q in [Path(__file__),OUT/'first-producer.py',OUT/'plan.json',OUT/'first-check-failure.json']})
    write(OUT/'corrected-plan.json',p)


worker=FunctionType(old.worker.__code__,dict(old.worker.__globals__,OUT=OUT,Material=Material,configure=configure))
pilot=FunctionType(old.pilot.__code__,dict(old.pilot.__globals__,OUT=OUT,worker=worker))
production=FunctionType(old.production.__code__,dict(old.production.__globals__,OUT=OUT,worker=worker))


def sources():
    # Same finite primitive pressure readout, now using the repaired flux owner.
    fn=FunctionType(old.sources.__code__,dict(old.sources.__globals__,OUT=OUT,GR=GR,Material=Material))
    fn()


charge=FunctionType(old.charge.__code__,dict(old.charge.__globals__,OUT=OUT,GR=GR,__file__=__file__))


if __name__=='__main__':
    signal.signal(signal.SIGALRM,old.old.repaired.forcing.history.flow.old.optical.timeout);action=sys.argv[1];started=time.monotonic()
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc),seconds=time.monotonic()-started));raise

"""Collocated source/readout SEM on the Phase93 native background.

Counterexample candidate: replace the spatial sampling/mass rule, never a
physical source or a failed gate. No further local interface patch is added.
"""
from pathlib import Path
import argparse
import inspect
import json
import resource
import signal
import textwrap
import time
import numpy as np
from scipy.sparse import diags
import def_native_coupled_readjustment as old

OUT=old.OUT.parent/'def-native-wave-collocation'
PRIOR=old.OUT
write=old.write

# Retain the entire Phase93 assembler. Replace its trial mesh in exactly one
# guarded location; no source, thermal mesh, heat map or boundary is replaced.
source=textwrap.dedent(inspect.getsource(old.Model.__init__))
anchor='self.local,self.polynomials=hierarchy.task.basis(degree);dx=np.diff(self.cells)'
assert source.count(anchor)==1
source=source.replace(anchor,"self.cells=np.unique(np.r_[self.edges,self.native,np.linspace(1,1.1,257)])\n    "+anchor)
namespace=dict(vars(old));exec(compile(source,__file__,'exec'),namespace)
assemble=namespace['__init__']


def symbolic():
    import sympy as s
    speed,E,F,P,e,p,q=s.symbols('v E F P e p q')
    boost=s.Matrix([[1,speed],[speed,1]])
    radiation=s.Matrix([[E,F],[F,P]])
    transformed=(boost*radiation*boost.T).applyfunc(lambda a:s.series(a,speed,0,2).removeO())
    flux=(transformed*s.Matrix([-speed,1])).applyfunc(lambda a:s.series(a,speed,0,2).removeO())
    assert flux==s.Matrix([F+speed*P,P+speed*F])
    material=(boost*s.Matrix([[e,q],[q,p]])*boost.T).applyfunc(lambda a:s.series(a,speed,0,2).removeO())
    material_flux=(material*s.Matrix([-speed,1])).applyfunc(lambda a:s.series(a,speed,0,2).removeO())
    assert material_flux.subs({q:F,p:P})==flux
    mu=s.symbols('mu',nonnegative=True)
    assert s.integrate(mu,(mu,0,1))==s.Rational(1,2)
    assert s.integrate(mu**2,(mu,0,1))==s.Rational(1,3)
    return dict(classification='Proven',passed=True,
        moving_surface='For a comoving radiation tensor(E,F,P), the stationary-frame energy/momentum fluxes through the moving surface are F+v*P and P+v*F to first order. Equal comoving heat flux and normal pressure match material and radiation regardless of the material rest density.',
        hemisphere='An isotropic outgoing hemisphere has F=c*E/2 and Pr=E/3=2F/(3c), so matching total material pressure to this radiation requires Pgas=0 or an explicit external stress carrier.',
        work='The emitted stationary-frame energy includes pressure work, not only comoving heat. In Phase93 units W=4*pi*R^3*A^4*N*a*P*zeta.',
        boundary_limit='A nonzero maintained gas pressure is an external stress condition, not an isolated material-vacuum interface. A small p*dV/heat ratio does not prove a charge error bound.',
        scope='Local linear Lorentz flux and moment identities. This is not yet a solved moving atmosphere, boosted ray history or final scalar charge.')


class Model(old.Model):
    def __init__(self,degree,lumped=True):
        start=time.monotonic();assemble(self,degree,resolved=True)
        self.mass_rule='GLL with exact first central cell' if lumped else 'original positive six-point composite mass'
        consistent=self.M
        if lumped:
            x=self.local
            weight=1/(degree*(degree+1)*np.polynomial.legendre.Legendre.basis(degree)(2*x-1)**2)
            widths=np.diff(self.cells);points=(self.cells[:-1,None]+widths[:,None]*x).ravel()
            weights=(widths[:,None]*weight).ravel()
            # r^2/r^4 spherical inertia vanishes at the centre. Integrate that
            # first element with positive interior quadrature, retaining both
            # centre degrees rather than imposing a new central constraint.
            weights[:degree+1]=0
            p=self.bg.sample(points);r=p['r'];b=1-2*p['m']/np.maximum(r,1e-100);c=p['N']*np.sqrt(b)
            W0=4*np.pi*r**4*np.exp(-8*p['phi']**2)*(p['e']+p['p'])/(b*c)
            W0[np.repeat(self.cells[:-1]>=1,degree+1)]=0
            W1=r*r/c;V,_=self.evaluation(points)
            mass=V[0].T@diags(weights*W0)@V[0]+V[1].T@diags(weights*W1)@V[1]
            gx,gw=np.polynomial.legendre.leggauss(6);r=self.cells[1]*(gx+1)/2;w=self.cells[1]*gw/2
            p=self.bg.sample(r);b=1-2*p['m']/r;c=p['N']*np.sqrt(b);V,_=self.evaluation(r)
            W0=4*np.pi*r**4*np.exp(-8*p['phi']**2)*(p['e']+p['p'])/(b*c);W1=r*r/c
            self.M=(mass+V[0].T@diags(w*W0)@V[0]+V[1].T@diags(w*W1)@V[1]).tocsc()
            assert np.all(self.M.diagonal()>0)
        self.mass_constant_controls=[]
        for field in range(2):
            q=np.zeros(self.size);indices=self.indices[::degree,field];q[indices[indices>=0]]=1
            baseline=q@(consistent@q);actual=q@(self.M@q)
            self.mass_constant_controls.append(float(abs(actual/baseline-1)))
        indices=np.searchsorted(self.cells,self.native)*degree
        assert np.array_equal(self.grid[indices],self.native)
        self.native_endpoint_indices=self.indices[indices]
        self.setup_seconds=time.monotonic()-start

    def evolve(self,*args,**kwargs):
        previous=old.OUT;old.OUT=OUT
        try:return super().evolve(*args,**kwargs)
        finally:old.OUT=previous


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Test and, if successful, fix the native-point wave readout by collocating every material face and native readout, with positive GLL inertia in the existing canonical weak equations.',
        cause_hypothesis='The rejected pointwise velocities arise from resolving short source-interface waves with cell-polynomial interpolation and consistent inertia. This is a hypothesis, not a confirmed cause.',
        changed='One global source/readout collocation rule, no additional local patch; GLL inertia except the regular central element. The old physical weak stiffness, extended thermal force and temperature/heat/photons are retained.',
        fixed='Phase93 native EOS, baryons, composition, horizon, physical closures and global max-norm gates. No field smoothing, RMS replacement or source retuning.',
        paths=['p4-8 pilot','p4-32','p4-64','p2-64','p4-64 consistent-mass contrast only if all primary gates pass'],
        gates=json.loads((PRIOR/'plan.json').read_text())['gates'],
        budget=dict(total_seconds=240,pilot_seconds=45,CPU_threads=1,memory_GB=5,native_calls=0,new_background_roots=0,automatic_expansion=False),
        decision='Use actual coupled paths and unchanged criteria. Stop the candidate on failed time/space or the measured budget; do not keep modifying this mesh or expand the horizon. The moving-interface identity is separate theorem progress, not a substitute for a converged charge.',
        references=['https://doi.org/10.1046/j.1365-246x.1999.00967.x','https://doi.org/10.1093/mnras/216.2.403'],
        reference_scope='The first paper supports the nodal GLL diagonal-mass principle, not our GR solver. The second studies radiating-star junctions; our hemispheric local flux identity is derived here rather than copied from its radial-null closure.',
        symbolic=symbolic(),bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(old.__file__),PRIOR/'inputs.npz',PRIOR/'coefficients.npz',PRIOR/'resolved-result.json',old.prior.OUT/'background.npz']}))
    (OUT/'reused-assembler.py').write_text(source,encoding='utf-8')


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert old.photons.digest(old.ROOT/p)==sha,p
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)));start=time.monotonic()
    horizon=json.loads((PRIOR/'plan.json').read_text())['horizon_seconds'];p=Model(4)
    pilot=p.evolve(horizon,8,'pilot',True)
    forecast=1.4*(3*p.setup_seconds+pilot['seconds']*168/8+10)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',forecast_seconds=forecast,
        setup_seconds=p.setup_seconds,evolution8_seconds=pilot['seconds'],dofs=p.size,
        mass_constant_controls=p.mass_constant_controls,memory_GB=pilot['memory_GB'],
        assumption='Measured p4 costs scaled to168 steps and three assemblies with40percent margin; p2 and the optional consistent-mass contrast are not measured.'))
    assert forecast<240,'Registered budget does not fit; no automatic expansion'
    rows=[p.evolve(horizon,n,f'p4-{n}',True) for n in [32,64]]
    a=np.load(OUT/'p4-32.npz');b=np.load(OUT/'p4-64.npz');compare={};passed=True
    for name in ['temperature','velocity','scalar']:
        err=float(np.max(abs(a[name]-b[name][::2]))/max(np.max(abs(b[name])),1e-100));limit=.02 if name=='temperature' else .03
        compare[name]=dict(time_relative=err,time_pass=err<limit,space_pass=False);passed &= err<limit
    if passed:
        del p
        other=Model(2);rows.append(other.evolve(horizon,64,'p2-64',True));a=np.load(OUT/'p2-64.npz')
        for name in compare:
            err=float(np.max(abs(a[name]-b[name]))/max(np.max(abs(b[name])),1e-100));limit=.02 if name=='temperature' else .03
            compare[name].update(space_relative=err,space_pass=err<limit);passed &= err<limit
    contrast=None
    if passed:
        assert time.monotonic()-start+forecast/2<240,'No time for required mass contrast'
        del other
        q=Model(4,lumped=False);q.evolve(horizon,64,'consistent-p4-64',True);a=np.load(OUT/'consistent-p4-64.npz');contrast={}
        for name in compare:
            err=float(np.max(abs(a[name]-b[name]))/max(np.max(abs(b[name])),1e-100));contrast[name]=err
            passed &= err<(.02 if name=='temperature' else .03)
    result=dict(classification='Counterexample candidate',passed=bool(passed),actual_coupled_evolution=True,
        comparisons=compare,consistent_mass_contrast=contrast,endpoint=rows[1]['history'][-1],seconds=time.monotonic()-start,
        memory_GB=max(r['memory_GB'] for r in rows),physical_moving_boundary_solved=False,final_dynamic_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print('COLLOCATION',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()

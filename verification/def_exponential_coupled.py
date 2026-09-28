"""Counterexample candidate: exponential scalar propagation with reciprocal work.

The prescribed boundary has an autonomous forcing reservoir. Its energy loss,
not an arbitrary scalar defect, supplies boundary work. Native matter is kept.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve
import sympy as sp
import def_nonstationary_interior as previous
import def_exact_scalar_carrier as carrier

s,e,ld=previous.s,previous.e,previous.ld
OUT=e.g.OUT/'def-exponential-coupled'


def solve(a,b):
    factor=lu_factor(np.asarray(a,float))
    x=lu_solve(factor,np.asarray(b,float)).astype(ld)
    for _ in range(4):x+=lu_solve(factor,np.asarray(b-a@x,float)).astype(ld)
    return x


def exponential(a):
    """Small dense extended-precision exponential; no native material calls."""
    scale=max(0,int(np.ceil(np.log2(max(float(abs(a).sum(1).max()),1)))))
    z=a/ld(2)**scale;answer=np.eye(len(a),dtype=ld);term=answer.copy()
    for j in range(1,80):
        term=term@z/j;answer+=term
        if abs(term).max()<ld('1e-24'):break
    else:raise RuntimeError('Exponential series cap')
    for _ in range(scale):answer=answer@answer
    return answer


def gradient(star,z):
    p=star.previous;n=star.n;size=n+9;R=star.rf[-1];sw=star.canonical_weight
    G=star.clock_matrix.copy();b=np.zeros(2*size,dtype=ld)
    c=star.scalar_area*z['kf']/star.distance/(4*np.pi*R)
    for j in range(1,n):
        v=np.zeros(2*size,dtype=ld);v[j-1]=-1/sw[j-1];v[j]=1/sw[j]
        G+=c[j]*np.outer(v,v)
    v=np.zeros(2*size,dtype=ld);v[n-1]=-1/sw[-1];v[n:n+9]=star.boundary_basis
    G+=c[-1]*np.outer(v,v)
    factor=-e.GRAV/R*z['H']*star.volume*star.beta/(p['a']+z['a'])
    G[np.arange(n),np.arange(n)]+=2*factor*z['a']*z['trace']/sw**2
    b[:n]=factor*(p['a']*p['trace']-z['a']*z['trace'])*p['psi']/sw
    G[size+np.arange(n),size+np.arange(n)]=z['k']
    return G,b


def setup(star,z,boundary,amplitude):
    n=star.n;size=n+9;R=star.rf[-1]
    star.canonical_weight=np.sqrt(star.scalar_volume/(4*np.pi*R**3))
    basis=np.zeros(9,dtype=ld);basis[[0,1,3,5,7]]=ld(amplitude)*np.array([35,-56,28,-8,1],dtype=ld)/128
    basis[0]+=boundary;star.boundary_basis=basis
    clock=np.zeros((2*size,2*size),dtype=ld)
    for j in range(1,5):
        a=n+2*j-1;b=a+1;omega=2*ld(str(np.pi))*j
        clock[a,size+b]=clock[size+b,a]=omega
        clock[b,size+a]=clock[size+a,b]=-omega
    star.clock_matrix=clock
    y=np.zeros(2*size,dtype=ld);y[:n]=star.canonical_weight*z['psi']
    y[size:size+n]=star.canonical_weight*R*z['Pi'];y[n+np.array([0,1,3,5,7])]=1
    z['canonical']=y;star.previous=z
    M,b=gradient(star,z);assert abs(b).max()==0
    Q=np.block([[np.zeros((size,size),dtype=ld),np.eye(size,dtype=ld)],
        [-np.eye(size,dtype=ld),np.zeros((size,size),dtype=ld)]])
    h=star.h*e.C/R;dim=len(M);aug=np.zeros((2*dim,2*dim),dtype=ld)
    aug[:dim,:dim]=h*Q@M;aug[:dim,dim:]=h*Q
    exp=exponential(aug);star.reference_matrix=M
    star.propagator=exp[:dim,:dim];star.phi_propagator=exp[:dim,dim:]
    star.reference_exponential_defect=float(abs(star.propagator.T@M@star.propagator-M).max())
    return y


def wave(star,z):
    old=star.previous['canonical'];G,b=gradient(star,z);D=G-star.reference_matrix
    B=star.phi_propagator;A=star.propagator
    # This form remains regular when A+I is singular (e.g. a half-period).
    matrix=np.eye(len(G),dtype=ld)-B@D/2
    rhs=(A+B@D/2)@old+B@b
    y=solve(matrix,rhs);mid=(old+y)/2;g=G@mid+b
    star.new_canonical=y;star.amplitude=star.boundary_basis@y[star.n:star.n+9]
    star.scalar_gradient_work=float(g@(y-old))
    n=star.n;size=n+9;R=star.rf[-1];sw=star.canonical_weight
    phi=y[:n]/sw;pi=y[size:size+n]/(R*sw)
    phir=np.r_[ld(0),np.diff(phi)/np.diff(star.r),(star.amplitude-phi[-1])/star.distance[-1]]
    error=float(abs(matrix@y-rhs).max()/max(abs(y).max(),ld('1e-30')))
    return phi,pi,phir,error


_evaluate=FunctionType(s.PulseStar.evaluate.__code__,dict(s.PulseStar.evaluate.__globals__,wave=wave))


class ExponentialStar(s.PulseStar):
    def evaluate(self,delta):
        z=_evaluate(self,delta)
        z['canonical']=self.new_canonical.copy()
        return z


def initialize(pool,h,amplitude):
    star,delta,z,boundary=previous.initialize(pool,h)
    star.__class__=ExponentialStar
    # Restore the actual initial metric rate instead of the old unused at=0.
    s.m.finish_fluxes(star,z)
    setup(star,z,boundary,amplitude)
    return star,delta,z,boundary


def symbolic():
    x0,x1,p0,p1,gx,gp,h=sp.symbols('x0 x1 p0 p1 gx gp h')
    assert sp.expand(gx*(h*gp)+gp*(-h*gx))==0
    M=sp.diag(2,3);Q=sp.Matrix([[0,1],[-1,0]])
    assert (Q*M).T*M+M*(Q*M)==sp.zeros(2)
    R,G,dm,clock,material=sp.symbols('R G dm clock material')
    assert sp.expand((dm/R+clock-G/R*material)*R/G-(dm/G+R/G*clock-material))==0
    return dict(classification='Proven',passed=True,
        scope='Skew work cancellation, reference quadratic generator identity, and dimensionless-to-physical mass/clock work scaling. Full discrete mass identity still requires the numerical coupled controls.')


def controls():
    star,_=s.initialize(s.old.imported.CachedOnly(),-4,ld('.001'),ld(1))
    z=dict(np.load(previous.OUT/'plus-12/initial.npz'));star.previous=z
    star.h=star.rf[-1]/e.C
    z['psi']=np.zeros(star.n,dtype=ld);z['Pi']=np.zeros(star.n,dtype=ld)
    z['k']=z['H']*z['b'];z['kf']=z['Hgrad']*z['bf']
    y=setup(star,z,ld(0),ld('1e-5'));phi,pi,phir,residual=wave(star,z)
    L,f,_=carrier.operator(star,z);exact=carrier.carrier(L,f,1e-5,1)
    relative=float(abs(phi-exact[:star.n]).max()/max(abs(exact[:star.n]).max(),1e-30))
    loss=-(star.new_canonical@star.clock_matrix@star.new_canonical-y@star.clock_matrix@y)/2
    physical=star.reference_matrix-star.clock_matrix
    gain=(star.new_canonical@physical@star.new_canonical-y@physical@y)/2
    balance=float(abs(gain-loss)/max(abs(loss),ld('1e-30')))
    row=dict(classification='Counterexample candidate',relative_carrier_error=relative,
        reservoir_energy_relative_error=balance,linear_solve_residual=residual,
        reference_exponential_defect=star.reference_exponential_defect,symbolic=symbolic())
    row['passed']=relative<2e-11 and balance<2e-14 and residual<2e-17
    assert row['passed'],row
    return row


def pilot():
    assert not OUT.exists();OUT.mkdir();began=time.monotonic()
    check=controls();e.write(OUT/'controls.json',check)
    files=[Path(__file__),Path(previous.__file__),Path(s.__file__),Path(s.prior.__file__),Path(s.old.__file__),
        previous.OUT/'manifest.json',previous.OUT/'plus-12/initial.npz']
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2ba5f1ea',
        bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        claim='First native reciprocal coupling of the exponential scalar carrier and boundary energy reservoir; no defect assigned to material heat.',
        gates=dict(native=1.,scalar=2e-17,mass_work_identity=2e-15),
        budget=dict(native_steps=2,hard_timeout_seconds=90,workers=4,automatic_expansion=False),
        scope='Same finite cavity; no free surface, radiative boundary or orbital/observational completion. Pilot is not a time convergence claim.'))
    history=[]
    try:
        with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
            R=ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm']);h=R/e.C/6
            star,delta,z,boundary=initialize(pool,h,ld('1e-5'))
            np.savez_compressed(OUT/'initial.npz',delta=delta,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            for j in [1,2]:
                before=z;star.previous=z;star.previous_boundary=star.amplitude;star.initial_guess=delta.copy()
                def log(row):
                    row['step']=j
                    with (OUT/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
                delta,z,value=s.stage(star,(delta,z),log)
                dy=z['canonical']-before['canonical'];mid=(z['canonical']+before['canonical'])/2
                work=-R/e.GRAV*(star.clock_matrix@mid)@dy
                dm=(z['dmf'][-1]-before['dmf'][-1])/e.GRAV
                res=np.sum(z['H']*star.volume*star.heat0*value[:,1])
                scale=np.sum(star.volume*(abs(z['Ephi'])+abs(before['Ephi'])+abs(z['dEm']-before['dEm'])))+abs(work)
                identity=float(abs(dm-work-res)/max(scale,ld('1e-100')))
                boundary_expected=boundary+ld('1e-5')*np.sin(ld(str(np.pi))*ld(j)/6)**8
                row=dict(step=j,native_norm=float(np.max(abs(value)/s.ATOL)),scalar=z['scalar_residual'],
                    mass_work_identity=identity,boundary_error=float(abs(star.amplitude-boundary_expected)),
                    clock_work_erg=float(work),native_calls=star.pool.evaluations,
                    scalar_gradient_work=star.scalar_gradient_work)
                history.append(row)
                np.savez_compressed(OUT/f'step-{j:03d}.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
                e.write(OUT/'progress.json',row);print('EXPONENTIAL',json.dumps(row),flush=True)
                assert identity<2e-15 and row['boundary_error']<2e-19,row
        result=dict(classification='Counterexample candidate',passed=True,history=history,
            seconds=time.monotonic()-began,full_dynamic_charge_solved=False)
        e.write(OUT/'result.json',result)
    except Exception as error:
        e.write(OUT/'failure.json',dict(error=repr(error),history=history,seconds=time.monotonic()-began));raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['controls','pilot']);action=p.parse_args().action
    print(json.dumps(controls())) if action=='controls' else pilot()

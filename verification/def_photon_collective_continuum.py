"""Continuous Maxwell dielectric spectrum; retain only unresolvable plasmons as poles.

The failed discrete-velocity calculation and its thresholds stay frozen. This
uses the same EOS, collision model, photon cells and coupled solver, replacing
only velocity atoms by an explicitly positive continuous density.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import signal
import time
import numpy as np
from scipy.special import dawsn
from scipy.optimize import brentq
import def_photon_collective as m

OUT=m.OUT/'continuum'


def response(y):
    y=np.asarray(y);small=abs(y)<20
    R=np.empty_like(y);derivative=np.empty_like(y)
    d=dawsn(y[small]/np.sqrt(2));R[small]=1-np.sqrt(2)*y[small]*d
    derivative[small]=-np.sqrt(2)*d-y[small]*R[small]
    z=y[~small];t=1/z**2
    R[~small]=-t*(1+t*(3+t*(15+t*(105+t*(945+t*(10395+t*135135))))))
    derivative[~small]=2*t/z*(1+t*(6+t*(45+t*(420+t*(4725+t*(62370+t*945945))))))
    imaginary=np.sqrt(np.pi/2)*y*np.exp(-y*y/2)
    return R+1j*imaginary,derivative


def dielectric(v,x,s):
    speeds=np.r_[1.,np.sqrt(m.base.m_e/s['ion_mass'])]
    strength=np.r_[1.,s['ion_strength']]
    y=np.asarray(v)[...,None]/speeds
    R,derivative=response(y)
    eps=1+(R*strength).sum(axis=-1)/x**2
    ep=(derivative*strength/speeds).sum(axis=-1)/x**2
    return R,eps,ep


def density(v,x,s):
    speeds=np.r_[1.,np.sqrt(m.base.m_e/s['ion_mass'])]
    R,eps,_=dielectric(v,x,s);ce=R[...,0]/x**2
    ci=(R[...,1:]*s['ion_strength']).sum(axis=-1)/x**2
    phi=np.exp(-np.asarray(v)[...,None]**2/(2*speeds**2))/(np.sqrt(2*np.pi)*speeds)
    return (phi[...,0]*abs(1+ci)**2+(phi[...,1:]*s['ion_strength']).sum(axis=-1)*abs(ce)**2)/abs(eps)**2


def measure(x,order,s):
    # Piecewise-linear POSITIVE density has exact polynomial cumulative moments.
    # No interpolation of independent CDFs and no renormalization of sum rules.
    edges=np.array([0,.001,.003,.01,.03,.1,.3,1.,2.,4.,8.,12.])
    v=np.unique(np.concatenate([np.linspace(a,b,order+1) for a,b in zip(edges[:-1],edges[1:])]))
    pole=0.;weight=0.;width=0.
    if dielectric(2.2,x,s)[1].real<0:
        center=brentq(lambda z:dielectric(z,x,s)[1].real,2.2,2/x+5,xtol=1e-13,rtol=1e-14)
        R,eps,ep=dielectric(center,x,s);width=float(eps.imag/ep)
        if width<1e-8:
            pole=center;weight=float(R[0].real**2/(x*x*center*ep))
            # For a Maxwell Landau width below 1e-8, retain a residue pole. The omitted core
            # and residue approximation are checked by independent sum rules
            # and dielectric responses; they are not an exact kinetic theorem.
            cut=.001*center
            if center-cut<12:v=np.unique(np.r_[v,center-cut,center,center+cut])
        else:
            # Log-distance knots resolve the wings as well as the narrow core.
            # Uniform arctan knots followed by linear interpolation in v create
            # huge false triangles at the last transformed interval.
            delta=np.geomspace(width*.001,2,int(np.ceil(np.log(2/(width*.001))*order/2))+1)
            extra=np.r_[center,center-delta,center+delta]
            v=np.unique(np.r_[v,extra[(extra>0)&(extra<12)]])
    f=density(v,x,s)
    if pole:f[abs(v-pole)<=.001*pole*(1+1e-12)]=0
    assert np.all(np.isfinite(f)) and np.all(f>=0)
    h=np.diff(v);slope=np.diff(f)/h;left=v[:-1]
    mass=f[:-1]*h+slope*h*h/2
    first=left*mass+f[:-1]*h*h/2+slope*h**3/3
    second=left*left*mass+left*(f[:-1]*h*h+2*slope*h**3/3)+f[:-1]*h**3/3+slope*h**4/4
    return dict(v=v,f=f,slope=slope,mass=np.r_[0.,np.cumsum(mass)],
        first=np.r_[0.,np.cumsum(first)],second=float(second.sum()+weight*pole*pole),
        pole=pole,weight=weight,width=width)


def cumulative(row,z):
    v=row['v'];f=row['f'];i=np.searchsorted(v,z,side='right')-1;i=np.clip(i,0,len(v)-2)
    t=np.clip(z-v[i],0,v[i+1]-v[i]);s=row['slope'][i]
    inc=f[i]*t+s*t*t/2
    M=row['mass'][i]+inc;J=row['first'][i]+v[i]*inc+f[i]*t*t/2+s*t**3/3
    take=z>=row['pole'];M=M+row['weight']*take;J=J+row['weight']*row['pole']*take
    return M,J


def overlap(table,node,a,b,c,d):
    if table is None:return m.overlap(None,node,a,b,c,d)
    row=table['rows'][node]
    (A,MA),(B,MB),(C,MC),(D,MD)=[cumulative(row,z) for z in [a,b,c,d]]
    return (MB-MA)-a*(B-A)+(b-a)*(C-B)+d*(D-C)-(MD-MC)


def spectrum(order,points):
    path=OUT/f'spectrum-{order}-{points}.npz'
    if path.exists():
        data=np.load(path);return dict(grid=data['grid'],rows=[{key:data[f'{j}-{key}'] for key in ['v','f','slope','mass','first','second','pole','weight','width']} for j in range(len(data['grid']))])
    s=np.load(m.OUT/'inventory.npz');grid=np.geomspace(1e-7,1e3,points);rows=[];errors=[]
    for x in grid:
        row=measure(x,order,s);rows.append(row)
        exact=(x*x+s['ion_strength'].sum())/(x*x+1+s['ion_strength'].sum())
        errors.append([abs(2*(row['mass'][-1]+row['weight'])/exact-1),abs(2*row['second']-1)])
    assert np.max(errors)<2e-4,('continuum sum rules',np.max(errors))
    arrays={'grid':grid}
    for j,row in enumerate(rows):arrays.update({f'{j}-{key}':value for key,value in row.items()})
    np.savez_compressed(path,**arrays)
    m.old.ex.write(path.with_suffix('.json'),dict(order=order,points=points,maximum_moment_errors=np.max(errors,axis=0).tolist(),fitted_normalization=False))
    return dict(grid=grid,rows=rows)


cell_density=FunctionType(m.cell_density.__code__,dict(vars(m),overlap=overlap))
raw_coefficients=FunctionType(m.coefficients.__code__,dict(vars(m),spectrum=spectrum,cell_density=cell_density),argdefs=m.coefficients.__defaults__)


def coefficients(bank,order,nangle,nfreq=1):
    i,j,a,info=raw_coefficients(bank,order,nangle,nfreq)
    info['model']='Positive piecewise-linear continuous Maxwell/RPA density plus residues for Landau width below 1e-8 sigma_e; exact cell-overlap moments; no opacity normalization.'
    return i,j,a,info


class Operator(m.Operator):
    __init__=FunctionType(m.previous.Operator.__init__.__code__,dict(vars(m.previous),coefficients=coefficients),argdefs=m.previous.Operator.__init__.__defaults__,closure=m.previous.Operator.__init__.__closure__)


def prepare():
    OUT.mkdir(exist_ok=True);assert not (OUT/'pilot-plan.json').exists()
    prior=json.loads((m.OUT/'result.json').read_text());assert prior['passed'] is False
    m.old.ex.write(OUT/'pilot-plan.json',dict(classification='Counterexample candidate',
        reason='Discrete velocity and angular differences exceed the 1e-6 response gate and are comparable to the proposed collective effect. Replace discrete continuum atoms by the positive Maxwell dielectric density. Do not densify the failed velocity matrix or change its verdict.',
        claim='Resolve the same collisionless RPA continuum and its very narrow plasmon contribution before any further coupled production.',
        controls='Static and second sum rules plus the independent positive-imaginary-frequency dielectric response at x=1e-7,.01,.1,.13,.2,.3,.5,1,10,1000. No normalization; fixed 32/64/128 subdivisions per interval.',
        gates=dict(moment_relative=2e-4,imaginary_response_relative=2e-4),hard_seconds=120,CPU_workers=1,GPU=False,new_EOS_calls=0,new_stellar_steps=0,
        decision='Only after the spectral controls pass, reuse the existing solver and compare the same response at unchanged 1e-6 gates. Measure and declare that production budget first.',
        bindings={p.relative_to(m.old.h.ROOT).as_posix():m.old.h.digest(p) for p in [Path(__file__),Path(m.__file__),m.OUT/'inventory.npz',m.OUT/'result.json']}))
    print('PREPARED CONTINUUM PILOT',flush=True)


def pilot(target='uniform-log-pilot.json'):
    assert not (OUT/target).exists();signal.alarm(120);start=time.monotonic()
    s=np.load(m.OUT/'inventory.npz');rows=[]
    for order in [32,64,128]:
        for x in [1e-7,.01,.1,.13,.2,.3,.5,1.,10.,1000.]:
            row=measure(x,order,s);exact=(x*x+s['ion_strength'].sum())/(x*x+1+s['ion_strength'].sum())
            errors=[float(abs(2*(row['mass'][-1]+row['weight'])/exact-1)),float(abs(2*row['second']-1))]
            # Positive eight-node integration of the linear density in each cell.
            g,w=m.leggauss(8);v=row['v'];h=np.diff(v);q=(v[:-1]+v[1:])[:,None]/2+h[:,None]*g/2
            f=row['f'][:-1,None]+row['slope'][:,None]*(q-v[:-1,None]);weight=h[:,None]*w*f
            residual=[]
            for y in [.001,.01,.1,1.,10.]:
                value=np.sum(weight*q*q/(q*q+y*y))+2*row['weight']*row['pole']**2/(row['pole']**2+y*y)
                residual.append(float(abs(value/m.continuum_response(x,y,s)-1)))
            rows.append(dict(order=order,x=x,moment_errors=errors,imaginary_errors=residual,pole=float(row['pole']),pole_weight=float(row['weight']),width=float(row['width']),cells=len(v)-1))
    result=dict(classification='Counterexample candidate',rows=rows,seconds=time.monotonic()-start,
        passed=bool(max(max(r['moment_errors']+r['imaginary_errors']) for r in rows if r['order']==128)<2e-4))
    m.old.ex.write(OUT/target,result);print('PILOT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot']);globals()[p.parse_args().action]()

"""Exact-in-time finite-grid direct scalar carrier on fixed matter/metric.

This is a reusable carrier propagator, not a replacement for coupled matter
work. The native coupled equations are not modified by this diagnostic.
"""
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import expm, eigh
import def_nonstationary_interior as n

s,e,ld=n.s,n.e,n.ld


def operator(star,z):
    tc=star.rf[-1]/e.C;k=z['H']*z['b']
    kf=z['Hgrad']*z['bf'];c=star.scalar_area*kf/star.distance
    d=np.diag(-(c[:-1]+c[1:]))+np.diag(c[1:-1],1)+np.diag(c[1:-1],-1)
    L=(tc*e.C)**2*k[:,None]/star.scalar_volume[:,None]*d
    potential=(tc*e.C)**2*k*z['H']*4*np.pi*e.GRAV*star.volume/star.scalar_volume*star.beta*z['trace']
    L+=np.diag(potential)
    f=np.zeros(star.n,dtype=ld);f[-1]=(tc*e.C)**2*k[-1]/star.scalar_volume[-1]*c[-1]
    mass=star.scalar_volume/k;mass/=mass.max()
    return np.asarray(L,float),np.asarray(f,float),np.asarray(mass,float)


def carrier(L,f,amplitude,tau):
    """Return phi and dphi/d(t/tc) through the end of a sin^8 pulse (0<=tau<=1)."""
    assert 0<=tau<=1 and len(L)==len(f)
    size=len(f);A=np.zeros((2*size+9,2*size+9))
    A[:size,size:2*size]=np.eye(size);A[size:2*size,:size]=L
    coefficients=np.array([35,-56,28,-8,1])/128*amplitude
    A[size:2*size,2*size]=f*coefficients[0]
    initial=np.zeros(len(A));initial[2*size]=1
    for j in range(1,5):
        col=2*size+2*j-1;omega=2*np.pi*j
        A[size:2*size,col]=f*coefficients[j]
        A[col,col+1]=-omega;A[col+1,col]=omega;initial[col]=1
    return (expm(A*tau)@initial)[:2*size]


def check():
    start=time.monotonic();plan=n.bindings();star,_=s.initialize(s.old.imported.CachedOnly(),-4,ld('.001'),ld(1))
    z=dict(np.load(n.OUT/'plus-12/initial.npz'));L,f,mass=operator(star,z)
    symmetric=np.sqrt(mass[:,None]/mass[None,:])*L
    assert np.max(abs(symmetric-symmetric.T))<1e-11
    eigenvectors=eigh(symmetric);eigenvalues,U=eigenvectors
    assert eigenvalues.max()<0
    amp=plan['drive_amplitude'];exact=carrier(L,f,amp,1)[:star.n]
    # Independent modal forced-oscillator integral for each Fourier carrier.
    frequencies=np.sqrt(-eigenvalues);forcing=U.T@(np.sqrt(mass)*f)
    modal=np.zeros(star.n)
    for j,coefficient in enumerate([35,-56,28,-8,1]):
        omega=2*np.pi*j;denominator=frequencies**2-omega**2
        assert abs(denominator).min()>1e-5
        modal+=amp*coefficient/128*forcing*(np.cos(omega)-np.cos(frequencies))/denominator
    independent=U@modal/np.sqrt(mass)
    assert np.max(abs(exact-independent))<1e-12*amp
    B=np.load(s.OUT/'geometry.npz')['baryons'];B=B/B.sum()
    norm=lambda x:float(np.sqrt(np.sum(B*x*x)))
    A=np.block([[np.zeros_like(L),np.eye(star.n)],[L,np.zeros_like(L)]])
    rhs=np.r_[np.zeros(star.n),f];rows=[]
    for steps in plan['time_steps']:
        dt=2/steps;y=np.zeros(2*star.n)
        advance=np.linalg.solve(np.eye(len(A))-dt*A/2,np.eye(len(A))+dt*A/2)
        drive=np.linalg.solve(np.eye(len(A))-dt*A/2,dt*rhs)
        for j in range(steps//2):
            g0=amp*np.sin(np.pi*j*dt)**8;g1=amp*np.sin(np.pi*(j+1)*dt)**8
            y=advance@y+drive*(g0+g1)/2
        p=np.load(n.OUT/f'plus-{steps}'/f'step-{steps//2:03d}.npz')['psi']
        m=np.load(n.OUT/f'minus-{steps}'/f'step-{steps//2:03d}.npz')['psi'];native=(p-m)/2
        rows.append(dict(steps=steps,midpoint_rms=norm(y[:star.n]),
            midpoint_relative_time_error=norm(y[:star.n]-exact)/norm(exact),
            native_minus_fixed_matter_relative=norm(native-y[:star.n])/max(norm(native),1e-100)))
    return dict(classification='Counterexample candidate',passed=True,rows=rows,exact_carrier_rms=norm(exact),
        expm_modal_relative_difference=norm(exact-independent)/norm(exact),
        highest_drive_phase_per_step=[float(8*np.pi*2/k) for k in plan['time_steps']],
        seconds=time.monotonic()-start,new_native_calls=0,new_fluid_steps=0,
        source_sha256=e.digest(Path(__file__)),
        scope='Exponential and independent modal integration of the fixed initial finite-grid scalar operator. Isolates direct-wave temporal error, and provides a reusable exact carrier. No modification or acceptance of failed coupled histories; time-dependent metric, native matter work and physical exterior remain to be coupled.')


if __name__=='__main__':
    target=n.OUT/'exact-carrier.json';assert not target.exists()
    result=check();e.write(target,result);print(json.dumps(result))

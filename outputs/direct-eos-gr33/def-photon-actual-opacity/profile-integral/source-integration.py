"""Integrate native LTE Voigt lines separately from the background spectrum.

No finer global spectrum: export line parameters at their own frequency, keep
the native H/He special profiles in the background, and integrate line areas.
"""
from collections import Counter
from pathlib import Path
from types import FunctionType
import argparse
import difflib
import json
import shutil
import signal
import subprocess
import time
import numpy as np
from scipy.special import wofz, erfc
import def_photon_actual_opacity as a
import def_photon_actual_opacity_validate as v
import def_photon_atomic_repair as r

OUT=a.OUT/'profile-integral';CACHE=a.CACHE/'profile-integral';BUILD=CACHE/'build'


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir();BUILD.mkdir()
    a.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Replace unresolved triangular interpolation of ordinary LTE atomic lines with profile-area integration on the EXISTING photon bins. Keep special H/He profiles in the native background.',
        decision='Only rerun the short coupled response if the fixed-bath source-grid and independent quadrature checks are below the original 1e-3 threshold. No global grid refinement or new stellar steps.',
        model='Ordinary LTE lines use full Voigt wings over the inherited finite photon band. Native opacity-dependent wing truncation is removed for these lines, a declared change of line representation, not an accuracy claim for the atomic model.',
        budget=dict(native_runs=3,native_seconds=120,integration_seconds=180,coupled_seconds=480,CPU_workers=1,memory_GB=4,new_stellar_steps=0),
        forecast='Prior same-size native spectra took 5.63/12.55/1.56s. Profile integration is unmeasured and capped at180s. Prior coupled response155.57s; capped480s. Stop on failed preflight; preserve original false verdict.',
        gates=dict(source_grid=.001,quadrature=.0001,time_order=1.8,time_error=.001),
        bindings={str(p):a.digest(p) for p in [Path(__file__),r.BUILD/'synspec54.f',r.BUILD/'PARAMS.FOR',a.OUT/'response.json',r.OUT/'bank.npz']}))
    for name in ['data','models','selected-full.19']:(CACHE/name).symlink_to(a.CACHE/name,target_is_directory=name!='selected-full.19')
    for p in r.BUILD.glob('*.FOR'):shutil.copyfile(p,BUILD/p.name)
    old=(r.BUILD/'synspec54.f').read_text();new=old
    start=new.index('      SUBROUTINE LINOP(');end=new.index('      SUBROUTINE LINOPW(')
    line=new[start:end]
    anchor='''         CALL PROFIL(IL,IAT,ID,AGAM)
         DOP1=DOPA1(IAT,ID)
         FR0=FREQ0(IL)'''
    assert line.count(anchor)==1
    line=line.replace(anchor,'''         FR0=FREQ0(IL)
C        The line width belongs to the line frequency, not a batch center.
         DSAVE=DOPA1(IAT,ID)
         DOP1=1.D0/(FR0*3.33564D-11*
     *      SQRT(1.651D8*TEMP(ID)/AMAS(IAT)+VTURB(ID)))
         DOPA1(IAT,ID)=DOP1
         CALL PROFIL(IL,IAT,ID,AGAM)
         DOPA1(IAT,ID)=DSAVE''')
    anchor='         if(ab0.le.0.and.lasdel) go to 100'
    assert line.count(anchor)==1
    line=line.replace(anchor,anchor+'''
C        Export ordinary LTE profiles; retain special profiles natively.
         IF(IMODE0.EQ.-3.AND.INNLT.EQ.0.AND.LPR) THEN
            WRITE(89,'(i8,1p9e26.17)') INDAT(IL),FR0,
     *         GF0(IL),EXCL0(IL),AB0/STIM(ID),DOP1,AGAM,
     *         FREQ(1),FREQ(2),TEMP(ID)
            GO TO 100
         END IF''')
    new=new[:start]+line+new[end:]
    (BUILD/'synspec54.f').write_text(new)
    (OUT/'profile.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='frequency-repair/synspec54.f',tofile='profile-integral/synspec54.f')))
    cmd=['gfortran','-std=legacy','-fno-automatic','-mcmodel=medium','-g','-fcheck=all','-fbacktrace','synspec54.f','-o','synspec54']
    t=time.monotonic();p=subprocess.run(cmd,cwd=BUILD,capture_output=True,text=True,timeout=90)
    (OUT/'build.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    a.write(OUT/'build.json',dict(seconds=time.monotonic()-t,binary_sha256=a.digest(BUILD/'synspec54'),source_sha256=a.digest(BUILD/'synspec54.f'),command=cmd))


native=FunctionType(a.run.__code__,dict(vars(a),OUT=OUT,CACHE=CACHE,BUILD=BUILD),argdefs=a.run.__defaults__)
spectrum=FunctionType(v.spectrum.__code__,dict(vars(v),OUT=OUT))
background=FunctionType(v.bank.__code__,dict(vars(v),OUT=OUT,spectrum=spectrum,
    populations=FunctionType(v.populations.__code__,dict(vars(v),OUT=OUT))))


def spectra():
    for name,start,stop,step in [('balanced-coarse',100.,20000.,2.),('balanced-fine',100.,20000.,.5),('infrared',20000.,1e7,20000.)]:
        native(name,metals=True,lines=True,fixed_state=True,wstart=start,wend=stop,step=step)
        shutil.copyfile(CACHE/name/'fort.89',OUT/name/'fort.89')
    background()
    (OUT/'bank.npz').rename(OUT/'background-bank.npz')


def lines(name):
    # Preserve coincident transitions: maximum multiplicity within a native
    # batch, never a simple set of line centers. Only remove repeated batches.
    batches={};params={}
    for band in [name,'infrared']:
        for row in np.loadtxt(OUT/band/'fort.89'):
            # Coincident centers/oscillator strengths may have distinct
            # damping constants; those are distinct profiles, not duplicates.
            key=tuple(float(f'{x:.12g}') for x in row[:7])
            if key in params:assert np.allclose(row[:7],params[key],rtol=2e-11,atol=0)
            params[key]=row[:7]
            batch=tuple(row[7:9]);batches.setdefault(batch,Counter())[key]+=1
    counts=Counter()
    for counter in batches.values():counts|=counter
    return np.array([list(params[key])+[count] for key,count in sorted(counts.items())])


def integrate(b,profiles,order,tail_budget=1e-4):
    """Positive cell quadrature; discarded wings have a uniform operator bound."""
    n=int((b['u']<=60).sum());edges=b['edges_u'][:n+1];T=float(b['T'])
    h=2*np.pi*a.previous.old.base.hbar;scale=h/(a.previous.old.base.k*T)
    center=profiles[:,1]*scale;width=scale/profiles[:,5]
    strength=profiles[:,4]*profiles[:,7];damping=profiles[:,6]
    weights=np.cbrt(strength*damping*width**2);weights/=weights.sum()
    factor=float(b['H'])*(1+float(b['Ci'][:n].sum()/b['Cm']))
    bounds=tail_budget/factor*weights
    radius=np.maximum(16,2*np.sqrt(strength*damping/(np.sqrt(np.pi)*bounds)))
    tail=strength*(damping/np.sqrt(np.pi)/((radius/2)**2+damping**2)
        +erfc(radius/2)/(np.sqrt(np.pi)*damping))
    total=np.zeros(n);g,w=np.polynomial.legendre.leggauss(order)
    capacity0=15*float(b['arad'])*T**3/np.pi**4
    stimulated_scale=4.79928144e-11/(T*scale)
    segments=0;begin=time.monotonic()
    for u0,du,U,ag,R in zip(center,width,strength,damping,radius):
        left=max(0.,u0-R*du);right=min(edges[-1],u0+R*du)
        if right<=left:continue
        first=max(0,np.searchsorted(edges,left,side='right')-1)
        last=min(n,np.searchsorted(edges,right))
        sd=du*max(1.,ag)
        ta,tb=np.arcsinh((np.array([left,right])-u0)/sd)
        # A unit interval in asinh coordinates resolves the core and wings
        # independently of the continuum source spacing and bin boundaries.
        knots=u0+sd*np.sinh(np.arange(np.ceil(ta),tb))
        mesh=np.unique(np.r_[left,edges[first+1:last],knots[(knots>left)&(knots<right)],right])
        t=np.arcsinh((mesh-u0)/sd);mid=(t[:-1]+t[1:])/2;half=np.diff(t)/2
        tq=mid[:,None]+half[:,None]*g
        uq=u0+sd*np.sinh(tq)
        x=(uq-u0)/du
        profile=wofz(x+1j*ag).real
        density=capacity0*uq**4*np.exp(-uq)/(-np.expm1(-uq))**2
        amount=(density*U*profile*(-np.expm1(-stimulated_scale*uq))
            *sd*np.cosh(tq)*half[:,None]*w).sum(axis=1)
        ids=np.searchsorted(edges,(mesh[:-1]+mesh[1:])/2,side='right')-1
        np.add.at(total,ids,amount);segments+=len(ids)
    return total,dict(seconds=time.monotonic()-begin,segments=segments,order=order,
        omitted_wing_response_bound=float(tail.sum()*factor),tail_budget=tail_budget,
        profiles=len(profiles),bound_scope='Exact Voigt model and exact positive bin integration; floating-point quadrature is checked separately. Streaming is skew and the common collision operator is dissipative.')


def bank():
    assert not (OUT/'bank.npz').exists();signal.alarm(180)
    b=dict(np.load(OUT/'background-bank.npz'));n=int((b['u']<=60).sum())
    coarse=np.load(OUT/'balanced-coarse-lines.npy');fine=np.load(OUT/'balanced-fine-lines.npy')
    assert coarse.shape==fine.shape and np.allclose(coarse,fine,rtol=2e-11,atol=0)
    a.write(OUT/'integration-plan.json',dict(classification='Counterexample candidate',
        method='Integrate 60140 native ordinary LTE profiles using positive Gauss4/8 in asinh frequency coordinates split at each inherited photon-cell edge. Full Voigt target, analytically bounded omitted wings rather than opacity-dependent native cutoff.',
        tail_bound_initial=1e-4,quadrature_gate=1e-4,source_plus_two_tails_gate=.001,
        seconds=180,new_stellar_steps=0,
        bindings={str(p):a.digest(p) for p in [Path(__file__),OUT/'background-bank.npz',OUT/'balanced-fine-lines.npy',OUT/'balanced-coarse-lines.npy']}))
    q4,info4=integrate(b,fine,4);print('QUADRATURE4',info4,flush=True)
    q8,info8=integrate(b,fine,8);print('QUADRATURE8',info8,flush=True)
    duration=float(b['H']/b['c']);C=b['Ci'][:n];Cm=float(b['Cm'])
    r4=b['rate_a_fine'][:n]+float(b['c'])*q4/C
    r8=b['rate_a_fine'][:n]+float(b['c'])*q8/C
    rc=b['rate_a_coarse'][:n]+float(b['c'])*q8/C
    def distance(x,y):return float(np.linalg.norm(np.sqrt(C)*(np.exp(-duration*x)-np.exp(-duration*y)))/np.sqrt(Cm))
    grid=distance(r8,rc);quad=distance(r4,r8)
    for label,rate in [('fine',r8),('coarse',rc),('quadrature4',r4)]:
        result=b['rate_a'].copy();result[:n]=rate;b['rate_a_'+label]=result
    np.savez_compressed(OUT/'bank.npz',**b)
    result=dict(classification='Counterexample candidate',passed=bool(grid+2*info8['omitted_wing_response_bound']<.001 and quad<1e-4),
        source_grid_fixed_bath_initial=grid,quadrature_fixed_bath_initial=quad,
        integration=[info4,info8],line_inventory_same=True,coupled_response_evolved=False,
        physical_opacity_certified=False,
        scope='Fixed-temperature material reservoir preflight only. These differences are not the coupled-response error or a certification of atomic populations, background special profiles, or subcell transport.')
    a.write(OUT/'preflight.json',result);print('PREFLIGHT',result,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','spectra','bank']);globals()[p.parse_args().action]()

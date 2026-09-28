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
import subprocess
import time
import numpy as np
from scipy.special import wofz
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
            key=tuple(row[:4]);p=row[4:7]
            if key in params:assert np.allclose(p,params[key],rtol=2e-13,atol=0),(key,p,params[key])
            params[key]=p
            batch=tuple(row[7:9]);batches.setdefault(batch,Counter())[key]+=1
    counts=Counter()
    for counter in batches.values():counts|=counter
    return np.array([list(key)+list(params[key])+[count] for key,count in sorted(counts.items())])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','spectra']);globals()[p.parse_args().action]()

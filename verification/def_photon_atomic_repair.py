"""Repair physical-frequency line windows and bind actual EOS nuclear masses.

Reuse the direct opacity runner and coupled validator without changing their
frozen failed results or running another stellar history.
"""
from pathlib import Path
from types import FunctionType
import argparse
import difflib
import json
import shutil
import subprocess
import time
import numpy as np
from scipy.constants import atomic_mass
import def_photon_actual_opacity as a
import def_photon_actual_opacity_validate as v

OUT=a.OUT/'frequency-repair';CACHE=a.CACHE/'frequency-repair';BUILD=CACHE/'build'


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir()
    paths=[Path(__file__),Path(a.__file__),Path(v.__file__),a.OUT/'response.json',
        a.OUT/'line-selection-counterexample.json',a.previous.OUT/'inventory.npz',
        a.CACHE/'build-balanced/synspec54.f',a.CACHE/'build-balanced/PARAMS.FOR']
    a.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Fix the physical-frequency selection of atomic line profiles and use the actual EOS nuclear mass convention. Reassess the existing coarse/fine spectra; no finer grid or longer evolution.',
        failure='Original response source-grid difference 0.0190532432 >0.001; raw verdict remains false.',
        budget=dict(native_runs=3,native_seconds_total=120,coupled_seconds=480,CPU_workers=1,memory_GB=4,new_stellar_steps=0),
        forecast='Previous direct coarse/fine/IR 5.69/12.59/1.61s; new bisection has unmeasured small cost. Prior coupled execution155.57s. Stop on bound or failed physics/algebra; do not expand the frequency step0.5 ceiling.',
        gates=dict(source_grid_initial=.001,time_initial=.001,time_order=1.8,EOS_stages=.001,native_LTE=1e-6),
        bindings={str(p):a.digest(p) for p in paths}))
    for name in ['data','models','selected-full.19']:(CACHE/name).symlink_to(a.CACHE/name,target_is_directory=name!='selected-full.19')


def build():
    assert not (OUT/'build.json').exists();BUILD.mkdir()
    for p in (a.CACHE/'build-balanced').glob('*.FOR'):shutil.copyfile(p,BUILD/p.name)
    old=(BUILD/'PARAMS.FOR').read_text()
    # D(1,I) contains relative atomic masses, so HMASS must be one atomic mass
    # unit in grams, not a hydrogen-atom mass. Use the same molar convention.
    mass_unit=1/a.previous.old.base.Avogadro
    new=old.replace('HMASS = 1.67333D-24',f'HMASS = {mass_unit:.17e}'.replace('e','D'))
    assert new!=old;(BUILD/'PARAMS.FOR').write_text(new)
    patch=list(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='balanced/PARAMS.FOR',tofile='repaired/PARAMS.FOR'))
    old=(a.CACHE/'build-balanced/synspec54.f').read_text();new=old
    s=np.load(a.previous.OUT/'inventory.npz');weights=s['neutral_masses']/atomic_mass
    active=s['ni'].sum(axis=1)>0;Z=a.previous.inventory_reader.g.d.CHARGES
    anchor='         AMAS(I)=D(1,I)';assert new.count(anchor)==1
    added='\n'.join(f'         IF(I.EQ.{int(z)}) AMAS(I)={float(w):.17e}'.replace('e','D') for z,w in zip(Z[active],weights[active]))
    new=new.replace(anchor,anchor+'\n'+added)
    start=new.index('      SUBROUTINE LINOP(');end=new.index('      SUBROUTINE LINOPW(')
    line=new[start:end]
    anchor='''         IJ1=int(MAX(float(IJCNTR(I))-XIJEXT,3.))
         IJ2=int(MIN(float(IJCNTR(I))+XIJEXT,float(NFREQS)))'''
    assert line.count(anchor)==1
    line=line.replace(anchor,'         CALL LINE_BOUNDS(FR0-EXT,FR0+EXT,IJ1,IJ2)')
    line=line.replace('IF(IJ1.GE.NFREQ.OR.IJ2.LE.2) GO TO 100','IF(IJ1.GT.IJ2) GO TO 100')
    new=new[:start]+line+new[end:]
    new+='''
      SUBROUTINE LINE_BOUNDS(LOW,HIGH,J1,J2)
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'SYNTHP.FOR'
      REAL*8 LOW,HIGH
C     FREQ(3:NFREQ) is descending, with no uniform-grid assumption.
      LFT=3
      RGT=NFREQ+1
   10 IF(LFT.LT.RGT) THEN
         MID=(LFT+RGT)/2
         IF(FREQ(MID).GT.HIGH) THEN
            LFT=MID+1
          ELSE
            RGT=MID
         END IF
         GO TO 10
      END IF
      J1=LFT
      LFT=3
      RGT=NFREQ+1
   20 IF(LFT.LT.RGT) THEN
         MID=(LFT+RGT)/2
         IF(FREQ(MID).GE.LOW) THEN
            LFT=MID+1
          ELSE
            RGT=MID
         END IF
         GO TO 20
      END IF
      J2=LFT-1
      END
'''
    # PARAMS makes L... logical and R... real; indices must be explicit.
    new=new.replace('      REAL*8 LOW,HIGH\nC     FREQ','      REAL*8 LOW,HIGH\n      INTEGER LFT,RGT,MID\nC     FREQ')
    (BUILD/'synspec54.f').write_text(new)
    patch+=list(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='balanced/synspec54.f',tofile='repaired/synspec54.f'))
    (OUT/'repair.patch').write_text(''.join(patch))
    cmd=['gfortran','-std=legacy','-fno-automatic','-mcmodel=medium','-g','-fcheck=all','-fbacktrace','synspec54.f','-o','synspec54']
    begin=time.monotonic();r=subprocess.run(cmd,cwd=BUILD,capture_output=True,text=True,timeout=90)
    (OUT/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    a.write(OUT/'build.json',dict(command=cmd,seconds=time.monotonic()-begin,
        binary_sha256=a.digest(BUILD/'synspec54'),patch_sha256=a.digest(OUT/'repair.patch'),
        mass_unit_g=mass_unit,nuclear_masses={str(int(z)):float(w) for z,w in zip(Z[active],weights[active])}))
    # One independent executable check for the bisection bounds on uneven grids.
    f=np.array([99.,90.,10.,9.99,9.,.1]);tests=0
    for lo in [-1.,.1,9.,9.995,50.,100.]:
        for hi in [lo,lo+.01,lo+200.]:
            first=np.searchsorted(-f,-hi,side='left');last=np.searchsorted(-f,-lo,side='right')
            assert np.array_equal(np.arange(first,last),np.flatnonzero((f>=lo)&(f<=hi)));tests+=1
    a.write(OUT/'bounds-check.json',dict(passed=True,cases=tests))


native=FunctionType(a.run.__code__,dict(vars(a),OUT=OUT,CACHE=CACHE,BUILD=BUILD),argdefs=a.run.__defaults__)
spectrum=FunctionType(v.spectrum.__code__,dict(vars(v),OUT=OUT))
populations=FunctionType(v.populations.__code__,dict(vars(v),OUT=OUT))
bank=FunctionType(v.bank.__code__,dict(vars(v),OUT=OUT,spectrum=spectrum,populations=populations))
response=FunctionType(v.run.__code__,dict(vars(v),OUT=OUT))


def spectra():
    for name,start,stop,step in [('balanced-coarse',100.,20000.,2.),('balanced-fine',100.,20000.,.5),('infrared',20000.,1e7,20000.)]:
        native(name,metals=True,lines=True,fixed_state=True,wstart=start,wend=stop,step=step)
    bank()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','build','spectra','response']);globals()[p.parse_args().action]()

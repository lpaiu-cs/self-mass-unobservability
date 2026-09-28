"""Diagnostic OP monochromatic reader; not an accepted photon collision model.

The full monochromatic record is used; no spectrum is inferred from means.
Atomic populations, opacity domains and scattering approximations stay explicit.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import json
import struct
import time
import numpy as np
from scipy.io import FortranFile
from scipy.integrate import simpson
import op_planck_connection as op
import def_ionic_structure_transport as conduction

ex=conduction.ex
base=conduction.base
h=conduction.h
OUT=conduction.OUT.parent/'def-photon-collision-response'
AREA=10**-16.55280
THOMSON_AU=2.37567e-8
FILES={}


@lru_cache(maxsize=24)
def mesh(z):
    path=op.SOURCE/'mono'/f'm{z:02}.mesh'
    with FortranFile(path,'r',header_dtype='<u4') as f:raw=f.read_record(np.uint8).tobytes()
    dv,n=struct.unpack('<fi',raw[:8]);u=np.frombuffer(raw[8:],dtype='<f4').astype(float)
    assert n==10000 and len(u)==n and np.all(np.diff(u)>0)
    FILES[str(path)]=h.digest(path)
    return dv,u


@lru_cache(maxsize=128)
def record(z,it,jn):
    path=op.SOURCE/'mono'/f'm{z:02}.{it:03}';dv,u=mesh(z)
    with FortranFile(path,'r',header_dtype='<u4') as f:
        head=f.read_record(np.uint8).tobytes()
        zz,tt,mass,umin,umax,ncoarse,ntot,dpack,jlo,jhi,jstep=struct.unpack('<iifffiifiii',head)
        assert zz==z and tt==it and ntot==len(u) and jn in range(jlo,jhi+1,jstep)
        for target in range(jlo,jn+1,jstep):
            raw=f.read_record(np.uint8).tobytes();j,epa,planck,ross,nlo,nhi=struct.unpack('<ifffii',raw[:24]);assert j==target
            ions=np.frombuffer(raw[24:],dtype='<f4').astype(float)
            n=int(f.read_ints('<i4')[0]);packed=f.read_record(np.uint8).tobytes()
        if n:
            a=np.frombuffer(packed,dtype=[('index','<i4'),('sigma','<f4')]);assert len(a)==n and a['index'][0]==1 and a['index'][-1]==ntot
            sigma=np.interp(np.arange(1,ntot+1),a['index'],a['sigma'])
        else:sigma=np.frombuffer(packed,dtype='<f4').astype(float);assert len(sigma)==ntot
    FILES[str(path)]=h.digest(path)
    return dict(z=z,it=it,jn=jn,epa=epa,planck=planck,ross=ross,ions=ions,nlo=nlo,nhi=nhi,
        u=u,dv=dv,sigma=sigma,packing_tolerance=dpack)


def probe():
    assert not OUT.exists();OUT.mkdir();rows=[]
    for z in [1,2,6]:
        r=record(z,170,58);u=r['u'];se=-np.expm1(-u);sigma=r['sigma'];weight=15/np.pi**4*u**3/np.expm1(u)
        scattering=r['epa']*THOMSON_AU
        options=[simpson(s*weight,x=u) for s in [sigma, sigma*se, sigma-scattering/se, sigma*se-scattering]]
        ross=1/(r['dv']*np.sum(1/sigma))
        rows.append(dict(z=z,it=170,jn=58,header_planck=r['planck'],header_ross=r['ross'],
            planck_variants=options,reconstructed_native_ross=ross,
            ross_relative=abs(ross/r['ross']-1),
            planck_relative=[abs(v/r['planck']-1) for v in options],
            minimum_stimulated_absorption=float(min(sigma*se-scattering)),
            packing_tolerance=r['packing_tolerance']))
        print('PROBE',z,'HEADERS',r['planck'],r['ross'],'PLANCK_VARIANTS',options,'ROSS',ross,'MIN_ABS',min(sigma*se-scattering),flush=True)
    FILES[str(Path(__file__).resolve())]=h.digest(Path(__file__))
    ex.write(OUT/'probe.json',dict(classification='Counterexample candidate',
        scope='Exploratory replay of three source-grid records, not a preregistered acceptance test or a stellar mixture calculation.',
        variants=['sigma','sigma*(1-exp(-u))','sigma-scattering/(1-exp(-u))','sigma*(1-exp(-u))-scattering'],
        checks=rows,source_sha256=FILES,photon_collision_model_constructed=False,
        conclusion='Tested reconstructions fail to recover all native Planck headers. No renormalization, mean-fitted kernel or physical surface closure is authorized by this probe. Frequency resolution, tails and conventions are not separately established as the unique cause.'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['probe']);globals()[p.parse_args().action]()

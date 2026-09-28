"""Resolve the three acoustic interfaces; retain the failed uniform-grid run."""
from pathlib import Path
import argparse
import inspect
import json
import resource
import signal
import textwrap
import time
import numpy as np
import def_photon_radial_gr as prior

FAILED=prior.OUT;OUT=FAILED.parent/'def-photon-radial-gr-resolved'
write=prior.write


def wave_scales():
    b,d,R,edges=prior.geometry();seconds=.008*float(b['H']/b['c'])
    rho,p,u,gamma=d['raw'][0,:5:1][[0,1,2,4]]
    c=float(b['c']);coordinate_clock=float(d['N'][0]/d['metric'][0])
    cs=c*np.sqrt(gamma*p/(rho*(c*c+u)+p))
    return cs*coordinate_clock*seconds/(100*R),c*coordinate_clock*seconds/(100*R)


class Background(prior.Background):
    def __init__(self,radiation,outer=2):
        super().__init__(radiation,outer)
        _,_,_,edges=prior.geometry();acoustic,light=wave_scales()
        # Fixed once from causality, not iterated in response to a verdict.
        extra=np.concatenate([np.linspace(edges[0]-1.25*light,edges[-1]+1.25*light,257),
            *[face+np.linspace(-1.25,1.25,33)*acoustic for face in edges]])
        self.grid=np.unique(np.r_[self.grid,extra]);self.surface_index=int(np.flatnonzero(self.grid==1)[0])
        self.nodes=self.sample(self.grid);self.mid=self.sample((self.grid[:-1]+self.grid[1:])/2)
        self.node_projection=radiation.projection(self.grid);self.mid_projection=radiation.projection(self.mid['r'])


# Preserve the frozen original implementation and replace only its readout
# quadrature, flux reconstruction and trial mesh. No new GR equation is added.
source=textwrap.dedent(inspect.getsource(prior.Mechanics.__init__))
old="""    gx,gw=np.polynomial.legendre.leggauss(12)
    r=np.concatenate([(lo+hi)/2+(hi-lo)*gx/2 for lo,hi in zip(self.edges[:-1],self.edges[1:])])
    dr=np.concatenate([(hi-lo)*gw/2 for lo,hi in zip(self.edges[:-1],self.edges[1:])])"""
new="""    inside=(m.points['r']>=self.edges[0])&(m.points['r']<self.edges[-1])
    r=m.points['r'][inside];dr=m.weights[inside];ids=(r>=self.edges[1]).astype(int)"""
assert old in source;source=source.replace(old,new)
old="""    self.volumes=weight.reshape(2,-1).sum(1)
    avg=sparse.csr_matrix((weight/np.repeat(self.volumes,len(gx)),(np.repeat(np.arange(2),len(gx)),np.arange(len(r)))),shape=(2,len(r)))"""
new="""    self.volumes=np.bincount(ids,weights=weight,minlength=2)
    avg=sparse.csr_matrix((weight/self.volumes[ids],(ids,np.arange(len(r)))),shape=(2,len(r)))"""
assert old in source;source=source.replace(old,new)
old="    shape=np.interp(r,self.edges,[0.,1.,0.],left=0,right=0)"
new="""    enclosed=self.sources(point)[2][:,:2].toarray()*(self.gr.C**4*R*(N*a)[:,None])/(self.gr.G*1e-7)
    shape=enclosed[:,0]/self.volumes[0]-enclosed[:,1]/self.volumes[1]
    shape[(r<=self.edges[0])|(r>=self.edges[-1])]=0."""
assert old in source;source=source.replace(old,new)
namespace=dict(vars(prior),Background=Background);exec(compile(source,'resolved-mechanics-source.py','exec'),namespace)


class Mechanics(prior.Mechanics):
    __init__=namespace['__init__']


class Coupled(prior.Coupled):
    def __init__(self,degree=4):
        tick=time.monotonic();self.m=Mechanics(degree)
        print('RESOLVED MECHANICS',degree,time.monotonic()-tick,flush=True)
        self.p=prior.Photons(self.m);print('RESOLVED PHOTONS',time.monotonic()-tick,flush=True)
        self.density_factor=prior.lu_factor(np.eye(2)-self.m.DS@self.p.Drho[:2])
        self.duration=.008*float(self.p.b['H']/self.p.b['c'])/self.m.tc
        self.maximum_port_residual=0.


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(120);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(6e9),int(6e9)))
    old=json.loads((FAILED/'result.json').read_text());assert old['time_passed'] and not old['spatial_passed']
    acoustic,light=wave_scales();_,_,R,edges=prior.geometry()
    plan=json.loads((FAILED/'plan.json').read_text())
    plan.update(reassessment='Uniform GR patch cells are about55 acoustic travel lengths wide on this short interval. The 65percent pressure-interface discrepancy is preserved. The old 12-point whole-volume compression readout also misses the acoustic layers; use every actual FEM quadrature row. Reconstruct the face momentum current by the same enclosed redshifted volume as the energy debit, instead of linear radius.',
        method='One fixed causal patch:32 subintervals over plus/minus1.25 acoustic travel lengths at each of3 pressure interfaces, and256 over the encompassing light cone. Keep original64 patch subintervals, full GR matrix, time32/64/128 and spatialp2/p4, photon table, duration and all original2percent gates. Do not refine again after failure.',
        scales=dict(acoustic_travel_cm=acoustic*100*R,light_travel_cm=light*100*R,
            original_patch_subinterval_cm=(edges[-1]-edges[0])*100*R/64),
        budget=dict(pilot_seconds=120,production_seconds=400,total_phase_seconds=900,CPU_threads=1,memory_GB=6,new_EOS_calls=0,automatic_expansion=False),
        bindings={str(p):prior.go.task.digest(p) for p in [Path(__file__),Path(prior.__file__),FAILED/'plan.json',FAILED/'result.json',
            prior.ph.OUT/'thermo.npz',prior.ph.OUT/'moments.npz',prior.ph.old.OUT/'bank.npz',prior.patch.OUT/'result.json']})
    write(OUT/'plan.json',plan);(OUT/'resolved-mechanics-source.py').write_text(source)
    write(OUT/'symbolic.json',prior.symbolic());prior.OUT=OUT;c=Coupled()
    write(OUT/'input.json',dict(classification='Counterexample candidate',**c.p.input,
        structural_checks=c.p.checks,coordinate_duration_seconds=c.duration*c.m.tc,GR_degrees_of_freedom=c.m.m.size))
    row=c.run(8,'pilot');elapsed=time.monotonic()-start
    forecast=1.4*(row['seconds']*480/8+2*(elapsed-row['seconds']))
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',seconds=elapsed,path_seconds=row['seconds'],
        forecast_seconds=forecast,maximum_production_seconds=400,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9))
    signal.alarm(0);print('RESOLVED PILOT',elapsed,forecast,flush=True)


def run():
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['forecast_seconds']<400
    source=inspect.getsource(prior.run).replace('signal.alarm(600)','signal.alarm(400)')
    namespace=dict(vars(prior),OUT=OUT,Coupled=Coupled)
    prior.OUT=OUT;exec(compile(source,'resolved-run-source.py','exec'),namespace);namespace['run']()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()

"""Regular radial quadrature: exact flat spherical discrete eigenmodes.

Counterexample candidate. New declared scalar quadrature, unchanged polarized
time update and numerical gates. The failed old boundary test stays failed.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import numpy as np
import def_scalar_hamiltonian as prior

ROOT,ld,PI=prior.ROOT,prior.ld,prior.PI
OUT=ROOT/'outputs/direct-eos-gr33/def-scalar-regular-radius'
digest=prior.digest


class Scalar(prior.Scalar):
    def __init__(self,n,gravity=True):
        super().__init__(n,gravity)
        # Midpoint quadrature of u=r*phi, p=r*Pi. The radial gradient
        # energy uses the exact piecewise-linear-u integral in each dual cell.
        self.volume=4*PI*self.r**2*np.diff(self.rf)
        self.w=self.volume/(4*PI)
        self.area=4*PI*np.r_[ld(0),self.r[:-1]*self.r[1:],self.r[-1]*self.rf[-1]]
        self.face_weight=self.area*self.distance
        self.fraction=np.full(n,ld('.5'))


evolve=FunctionType(prior.evolve.__code__,dict(vars(prior),Scalar=Scalar))
operator_check=FunctionType(prior.check.__code__,dict(vars(prior),Scalar=Scalar))


def check():
    result=operator_check();s=Scalar(32,False);old=prior.Scalar(32,False)
    phi=np.sinc(s.r);dx=ld(1)/s.n
    eigenvalue=4*np.sin(PI*dx/2)**2/dx**2
    new_error=np.max(abs(s.divergence(s.gradient(phi))+eigenvalue*phi))
    old_error=np.max(abs(old.divergence(old.gradient(phi))+eigenvalue*phi))
    # The imported pi/sinc arguments originate in binary64; the second
    # difference magnifies their boundary roundoff by O(dx^-2).
    assert new_error<ld('1e-12') and old_error>ld('.1'),(new_error,old_error)
    return dict(result,flat_eigen_identity_max_defect=float(new_error),
        predecessor_flat_eigen_defect=float(old_error),
        scope='Polarized mass identity and zero-momentum regularity retained; the complete flat radial operator including both boundaries has the exact cell-centered sine eigenmode.')


def save(name,value):return prior.save.__class__(prior.save.__code__,dict(vars(prior),OUT=OUT))(name,value)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,sha in p['bindings'].items():assert digest(ROOT/rel)==sha,rel
    return p


def prepare():
    assert not OUT.exists()
    old=json.loads((prior.OUT/'plan.json').read_text());failed=json.loads((prior.OUT/'result.json').read_text())
    assert failed['passed'] is False
    paths=[Path(__file__),Path(prior.__file__),prior.OUT/'plan.json',prior.OUT/'result.json',prior.OUT/'manifest.json']
    p=dict(old,bindings={p.relative_to(ROOT).as_posix():digest(p) for p in paths},
        source_sha256=digest(Path(__file__)),operator_check=check(),
        correction='Replace geometric scalar cell weights by midpoint r^2*dr and dual gradient weights by r_left*r_right*distance; the scalar field u=r*phi then has the exact flat cell-centered Laplacian, including the half boundary cells. Use midpoint mass interpolation. This is an explicit new spatial quadrature, not a threshold change.',
        budget=dict(cpu=1,blas_threads=1,gpu=False,full_timeout_seconds=30,maximum_runs=1,
            expected_wall_seconds=[3,15],basis='The same eight frozen-size paths with the predecessor operator took 2.82s. No new grids, durations or native EOS calls.'))
    OUT.mkdir();save('plan.json',p);save('pilot.json',dict(passed=True,reused_timing=failed['seconds']))
    print('PREPARED regular radial quadrature',json.dumps(p['operator_check']),flush=True)


run=FunctionType(prior.run.__code__,dict(vars(prior),Scalar=Scalar,evolve=evolve,OUT=OUT,bindings=bindings,save=save))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['check','prepare','run'])
    print(globals()[parser.parse_args().action]())

"""Conservative C1 enclosed-heat reconstruction in the actual native GR loop.

Counterexample candidate: a new spatial source approximation, not a changed
verdict for the rejected piecewise-constant source. Every face inventory stays.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags, vstack
import def_native_wave_collocation as prior

old=prior.old
OUT=prior.OUT.parent/'def-native-conservative-source'
write=old.write


class Model(prior.Model):
    def source_maps(self,r):
        if not hasattr(self,'secants'):
            # Difference first: a spatially constant enclosed inventory gives
            # zero local heat. Longdouble avoids adding a large common E offset
            # before taking this difference in the assembled source maps.
            rows=np.repeat(np.arange(self.n),2);cols=np.c_[np.arange(self.n)-1,np.arange(self.n)].ravel()
            valid=cols>=0
            difference=coo_matrix((np.tile([-1.,1.],self.n)[valid],(rows[valid],cols[valid])),shape=(self.n,self.n)).astype(np.longdouble).tocsr()
            self.secants=diags(1/self.volumes.astype(np.longdouble))@difference
            vl,vr=self.volumes[:-1],self.volumes[1:]
            inner=diags(vr/(vl+vr))@self.secants[:-1]+diags(vl/(vl+vr))@self.secants[1:]
            self.face_slopes=vstack([self.secants[0],inner,self.secants[-1]]).tocsr()
        ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,self.n-1)
        base=old.Model.face_map(self,r).astype(np.longdouble)
        x=np.asarray(base[np.arange(len(r)),ids]).ravel()
        vl=self.volumes[ids].astype(np.longdouble)
        dl=self.face_slopes[ids]-self.secants[ids];dr=self.face_slopes[ids+1]-self.secants[ids]
        enclosed=base+diags(vl*x*(1-x)**2)@dl+diags(vl*x*x*(x-1))@dr
        derivative=self.secants[ids]+diags(3*x*x-4*x+1)@dl+diags(3*x*x-2*x)@dr
        return enclosed.tocsc(),derivative.tocsc()

    def face_map(self,r):return self.source_maps(r)[0]

    def sources(self,p):
        r=p['r'];inside=r<=1;A4=np.exp(-8*p['phi']**2);b=1-2*p['m']/np.maximum(r,1e-100)
        enclosed,derivative=self.source_maps(r);geo=old.G*self.bg.R**2/old.C**4
        loss=diags(geo/A4*inside)@derivative
        # Gamma*P*rho_ref = -ad*heat_density is the native cp/cv identity.
        # Interpolate ad, then enforce that identity pointwise. Interpolating
        # its factors separately needlessly introduces pressure-source jumps.
        ad=np.interp(r,self.native,self.thermo[:,4]);ratio=np.zeros_like(r)
        np.divide(ad,p['gamma']*p['p'],out=ratio,where=inside)
        rr=diags(ratio)@loss
        J=diags(-old.G/(old.C**4*self.bg.R)*np.sqrt(b)/p['N']*inside)@enclosed
        return rr.tocsc(),loss.tocsc(),J.tocsc()

    def evolve(self,*args,**kwargs):
        previous=old.OUT;old.OUT=OUT
        try:return old.Model.evolve(self,*args,**kwargs)
        finally:old.OUT=previous


def controls():
    import sympy as s
    x,E0,E1,v,sl,sr=s.symbols('x E0 E1 v sl sr',nonzero=True)
    sec=(E1-E0)/v
    enclosed=E0+x*(E1-E0)+v*x*(1-x)**2*(sl-sec)+v*x*x*(x-1)*(sr-sec)
    assert s.simplify(enclosed.subs(x,0)-E0)==0 and s.simplify(enclosed.subs(x,1)-E1)==0
    derivative=s.diff(enclosed,x)/v
    assert s.simplify(derivative.subs(x,0)-sl)==0 and s.simplify(derivative.subs(x,1)-sr)==0
    assert s.simplify(s.integrate(derivative*v,(x,0,1))-(E1-E0))==0
    rp,rt,cp,cv,ad,gamma,P,rho=s.symbols('rp rt cp cv ad gamma P rho',nonzero=True)
    # ad=(P/rho)*chi_T/cvT, chi_T=-rho_T/rho_P,
    # Gamma=(cp/cv)/rho_P. Hence Gamma*P*(-rho_T/(rho*cpT))=ad.
    expression=gamma*P*(-rt/(rho*cp))-ad
    assert s.simplify(expression.subs({gamma:cp/(cv*rp),ad:-P*rt/(rho*rp*cv)}))==0
    m=Model.__new__(Model);m.n=3;m.edges=np.array([0.,.2,.6,1.]);m.bg=SimpleNamespace(R=1.)
    m.bg.sample=lambda r:dict(r=r,m=np.zeros_like(r),N=np.ones_like(r),phi=np.zeros_like(r))
    m.volumes=4*np.pi*np.diff(m.edges**3)/3
    gx,gw=np.polynomial.legendre.leggauss(4);fraction=(gx+1)/2
    radii=(m.edges[:-1,None]**3+np.diff(m.edges**3)[:,None]*fraction)**(1/3)
    face,_=m.source_maps(m.edges);_,grad=m.source_maps(radii.ravel())
    inventory=np.array([2.,-3.,5.],np.longdouble)
    boundary=face@inventory;density=(grad@inventory).reshape(3,4)
    increments=m.volumes/2*np.sum(density*gw,axis=1)
    error=float(max(abs(increments-np.diff(np.r_[np.longdouble(0),inventory]))))
    assert error<1e-12 and max(abs(boundary-np.r_[0.,inventory]))<1e-12
    return dict(classification='Proven',passed=True,
        identity='Cubic Hermite enclosed energy matches every thermal face and its common dE/dV slope. Its derivative integrates to the unchanged cell heat debit.',
        thermodynamic='Gamma*P*rho_ref=-ad*heat_density, from the same native cp/cv identities.',
        numerical_classification='Counterexample candidate',flat_cell_inventory_absolute=error,
        scope='Conservation and local thermodynamic algebra; not a physical subcell-profile certificate, whole response convergence or charge bound.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='88b332112',
        claim='Replace discontinuous pressure forcing by a conservative continuous subcell heat reconstruction in the actual native coupled evolution, testing whether this resolves the rejected pointwise GR velocity.',
        changed='One C1 cumulative heat profile in redshift proper volume. Face slopes are fixed volume-weighted adjacent secants; no fitted smoothing width. The native entropy-to-pressure identity is enforced between native points. All heat loss, enclosed metric source and heat lift use the same profile.',
        preserved='All material face energies, native EOS inputs, background, composition, heat law, emitted history, p2/p4 grid, original max-norm criteria and 3.434431ms horizon.',
        boundary='The same fixed emitting geometry and specified finite gas pressure. This run does not apply Phase94 moving-surface work or close the external stress reservoir.',
        scientific_boundary='This changes the subcell source approximation; it cannot relabel Phase93/94 failures or certify physical interpolation independence. Any eventual charge claim needs its source-profile error controlled.',
        paths=['p4-8 pilot','p4-32','p4-64','p2-64 only after temporal gates'],
        gates=json.loads((old.OUT/'plan.json').read_text())['gates'],
        budget=dict(total_seconds=240,CPU_threads=1,memory_GB=5,native_calls=0,background_roots=0,automatic_expansion=False),
        decision='Stop this fixed candidate on a failed gate or measured budget. Do not adjust slopes, source width, mesh, period or thresholds after seeing the response.',
        controls=controls(),bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(prior.__file__),Path(old.__file__),old.OUT/'inputs.npz',old.OUT/'coefficients.npz',old.prior.OUT/'background.npz',prior.OUT/'result.json']}))


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for name,h in plan['bindings'].items():assert old.photons.digest(old.ROOT/name)==h,name
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)));start=time.monotonic()
    horizon=json.loads((old.OUT/'plan.json').read_text())['horizon_seconds'];m=Model(4)
    native=m.bg.sample(m.native);direct=-m.raw[:,8]/(m.raw[:,0]*m.thermo[:,5])
    identity=m.thermo[:,4]/(native['gamma']*native['p']/(old.G*m.bg.R**2/old.C**4))
    thermo_error=float(max(abs(identity/direct-1)));assert thermo_error<1e-7
    pilot=m.evolve(horizon,8,'pilot',True)
    forecast=1.5*(3*m.setup_seconds+pilot['seconds']*168/8+10)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',forecast_seconds=forecast,
        setup_seconds=m.setup_seconds,evolution_seconds=pilot['seconds'],native_pressure_identity_relative=thermo_error,
        assumption='Measured p4 assembly/eight-step cost, 50percent margin. The p2 path is not measured.'))
    assert forecast<240,'Measured candidate does not fit; no automatic expansion'
    reports=[m.evolve(horizon,n,f'p4-{n}',True) for n in [32,64]]
    a=np.load(OUT/'p4-32.npz');b=np.load(OUT/'p4-64.npz');compare={};passed=True
    for field in ['temperature','velocity','scalar']:
        error=float(np.max(abs(a[field]-b[field][::2]))/max(np.max(abs(b[field])),1e-100));limit=.02 if field=='temperature' else .03
        compare[field]=dict(time_relative=error,time_pass=error<limit,space_pass=False);passed &= error<limit
    if passed:
        del m
        other=Model(2);reports.append(other.evolve(horizon,64,'p2-64',True));a=np.load(OUT/'p2-64.npz')
        for field in compare:
            error=float(np.max(abs(a[field]-b[field]))/max(np.max(abs(b[field])),1e-100));limit=.02 if field=='temperature' else .03
            compare[field].update(space_relative=error,space_pass=error<limit);passed &= error<limit
    result=dict(classification='Counterexample candidate',passed=bool(passed),comparisons=compare,
        actual_coupled_evolution=True,endpoint=reports[1]['history'][-1],seconds=time.monotonic()-start,
        memory_GB=max(r['memory_GB'] for r in reports),physical_source_profile_certified=False,
        moving_surface_solved=False,final_dynamic_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print('CONSERVATIVE SOURCE',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()

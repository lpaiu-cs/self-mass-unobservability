"""Causal hemispheric photons and their linear exterior GR/scalar response.

Counterexample candidate: the additive exterior response has dphi(surface)=0.
It is the forcing term for a future interior boundary match, not a whole star.
"""
from pathlib import Path
import argparse
import hashlib
import json
import signal
import time
import numpy as np
from scipy.integrate import simpson
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_radiative_exterior as vacuum

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'outputs/direct-eos-gr33'
OUT = BASE/'def-native-photon-exterior'
C = 2.99792458e10
G = 6.67430e-8


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')


def symbolic():
    import sympy as s
    r, m, v, f, fr, J, E, P = s.symbols('r m v f fr J E P', nonzero=True)
    b = 1-2*m/r
    mr = r*r*b*v*v/2
    vr = -(2/r+2*m/(r*r*b))*v
    Jr = 4*s.pi*r*r*E-r*v*v*J
    dm = r*r*b*v*f+J
    dl = dm/(r*b)
    dlr = sum(s.diff(dl,x)*y for x,y in [(r,1),(m,mr),(v,vr),(f,fr),(J,Jr)])
    dnr = dm/(r*r*b*b)+r*v*fr+4*s.pi*r*P/b
    target = 2*v*f/b+2*J/(r*r*b*b)-4*s.pi*r*(E-P)/b
    assert s.simplify(dnr-dlr-target) == 0
    mu = s.symbols('mu', nonnegative=True)
    energy = 2/(1+mu)
    radial = 2*(1+mu+mu**2)/(3*(1+mu))
    tangential = (energy-radial)/2
    assert s.simplify(-energy+radial+2*tangential) == 0
    return dict(classification='Proven', passed=True,
        mass='dm=r^2*b*Phi*dphi+J; J_prime+r*Phi^2*J=4*pi*r^2*E. E,P are Einstein-frame geometric photon stresses.',
        wave='(r^2/C)*dphi_tt-(r^2*C*dphi_r)_r-2*r^2*C*Phi^2*dphi/b = 2*C*Phi*J/b^2-4*pi*r^3*C*Phi*(E-P)/b; C=N*sqrt(b).',
        conservation='J=-q*sum(w*H(t-delay))/(N*a); w=2*mu_surface*dmu_surface; q=G*L_infinity/c^5. Body debit=-q*H(t), photon energy inside r=q*(H(t)-sum(w*H(t-delay))).',
        trace='-E+Pr+2*Pt=0 does not remove metric-mediated scalar forcing.',
        boundary='For a steady outgoing isotropic hemisphere E=2F/[c(1+mu_min)] and Pr=2F*(1+mu_min+mu_min^2)/[3c(1+mu_min)]; mu_min^2=1-N(r)^2/(N_surface^2*r^2).')


def history(t, integral=False):
    """Compact diagnostic luminosity, not a claimed stellar driving mechanism."""
    x = np.clip(t, 0, 2)
    if integral:
        return 3*x/8-np.sin(np.pi*x)/(2*np.pi)+np.sin(2*np.pi*x)/(16*np.pi)
    return np.sin(np.pi*x/2)**4


class Exterior:
    def __init__(self, name='fine', flat=False):
        d = np.load(BASE/'def-native-radiative-envelope'/f'{name}.npz')
        self.R = float(d['r'][-1])
        self.vb = float(d['v'][-1])*self.R
        self.mu = float(d['m'][-1])/self.R
        self.A = float(d['A'][-1])
        self.old_Nb = float(d['N'][-1])
        self.flat = flat
        if flat:
            self.M, self.K = 0., self.vb
        else:
            self.M, self.K, self.solution = vacuum.background(self.mu,self.vb,2e-13)
        self.Nb = float(self.metric(np.array([1.]))[1][0])
        # Preserve the saved local flux under the exterior clock normalization.
        self.L = float(d['Linfinity'])*(self.Nb/self.old_Nb)**2
        self.q = G*self.L/C**5
        self.scalar_scale = self.q*self.vb
        self.Pgas = float(d['Pgas'][-1])
        self.Prad = float(d['Prad'][-1])

    def metric(self, r):
        if self.flat:
            return np.zeros_like(r),np.ones_like(r),np.ones_like(r),np.ones_like(r),self.vb/r**2
        mass, c = self.solution.sol(1/r)
        b = 1-2*mass/r
        n = c/np.sqrt(b)
        return mass,n,b,c,self.K/(c*r*r)

    def rays(self, r, angles=64, order=32):
        x,w = np.polynomial.legendre.leggauss(angles)
        mu0 = (x+1)/2
        weights = w*mu0  # emitted power, NOT uniform angle inventory
        impact2 = (1-mu0**2)/self.Nb**2
        _,n,_,_,_ = self.metric(r)
        mu = np.sqrt(1-impact2[None,:]*(n/r)[:,None]**2)
        # Use the flat-ray longitudinal coordinate to regularize grazing rays.
        radial_span = ((r-1)*(r+1))[:,None]
        z1 = np.sqrt(radial_span+mu0**2)
        span = radial_span/(z1+mu0)
        y,v = np.polynomial.legendre.leggauss(order)
        z = mu0[None,:,None]+span[:,:,None]*(y+1)/2
        radius = np.sqrt(z*z+(1-mu0**2)[None,:,None])
        _,nn,_,cc,_ = self.metric(radius.ravel())
        nn,cc = nn.reshape(radius.shape),cc.reshape(radius.shape)
        mu_local = np.sqrt(1-impact2[None,:,None]*nn*nn/(radius*radius))
        integrand = z/(radius*cc*mu_local)
        delay = span/2*np.sum(integrand*v,axis=2)
        assert np.all(delay >= 0) and abs(weights.sum()-1)<1e-14
        return mu,weights,delay


class Response:
    def __init__(self, name, cells, angles, dt_fraction, end=6.):
        self.ext = e = Exterior(name)
        self.r = r = np.linspace(1,9,cells+1)
        self.dx = dx = 8/cells
        self.t = np.linspace(0,end,int(np.ceil(end/(dt_fraction*dx)))+1)
        self.dt = float(self.t[1])
        self.cells,self.angles = cells,angles
        # The same two-point positive FEM used in the canonical GR formulation.
        z = np.array([(1-1/np.sqrt(3))/2,(1+1/np.sqrt(3))/2])
        self.rq = rq = (r[:-1,None]+dx*z).ravel()
        _,n,b,c,v = e.metric(rq)
        self.mu,self.weights,self.delay = e.rays(rq,angles)
        self.n,self.b,self.c,self.v = n,b,c,v
        shape = np.array([1-z,z]).T
        self.shape = shape
        grad = np.array([-1.,1.])/dx
        p = rq*rq*c
        potential = -2*rq*rq*c*v*v/b
        inertia = rq*rq/c
        def matrix(coeff, gradient=False):
            local = dx/2*np.einsum('cq,qi,qj->cij',coeff.reshape(cells,2),
                np.broadcast_to(grad,(2,2)) if gradient else shape,
                np.broadcast_to(grad,(2,2)) if gradient else shape)
            diagonal = np.r_[local[0,0,0],local[:-1,1,1]+local[1:,0,0],local[-1,1,1]]
            return diags([local[:,1,0],diagonal,local[:,0,1]],[-1,0,1],format='csc')[1:-1,1:-1]
        self.M = matrix(inertia)
        self.K = matrix(p,True)+matrix(potential)
        self.solve = splu(self.M+self.dt**2/4*self.K).solve

    def source(self, t):
        h = history(t-self.delay)
        integ = history(t-self.delay,True)
        momE = (h/self.mu)@self.weights
        momP = (h*self.mu)@self.weights
        cumulative = integ@self.weights
        j = -cumulative*np.sqrt(self.b)/self.n
        ep = (momE-momP)/(4*np.pi*self.rq**2*self.n**2)
        load = 2*self.c*(self.v/self.ext.vb)*j/self.b**2-4*np.pi*self.rq**3*self.c*(self.v/self.ext.vb)*ep/self.b
        local = self.dx/2*np.einsum('cq,qi->ci',load.reshape(self.cells,2),self.shape)
        return (local[:-1,1]+local[1:,0])

    def run(self):
        start = time.monotonic()
        u = np.zeros(self.cells-1);vel=u.copy();acc=u.copy()
        force = self.source(0.)
        traces=[];work=0.;max_energy_error=0.;max_energy=0.
        boundary = []
        for t in self.t[1:]:
            nxt = self.source(t)
            predicted = u+self.dt*vel+self.dt**2/4*acc
            unew = self.solve(self.M@predicted+self.dt**2/4*nxt)
            anew = 4*(unew-predicted)/self.dt**2
            vnew = vel+self.dt*(acc+anew)/2
            work += (nxt+force)@(unew-u)/2
            energy = .5*(vnew@(self.M@vnew)+unew@(self.K@unew))
            max_energy=max(max_energy,abs(energy))
            max_energy_error=max(max_energy_error,abs(energy-work))
            u,vel,acc,force = unew,vnew,anew,nxt
            full = np.r_[0.,u,0.]
            # Second-order one-sided derivative of the actual solved field.
            slope=(4*full[1]-full[2])/(2*self.dx)
            traces.append([t,slope,2*np.interp(2,self.r,full),3*np.interp(3,self.r,full)])
            boundary.append(float(max(abs(full[self.r>=8.]))))
        traces=np.vstack([np.zeros(4),traces])
        energy_error=max_energy_error/max(max_energy,1e-30)
        assert energy_error<1e-8,energy_error
        return dict(traces=traces,r=self.r,field=np.r_[0.,u,0.],velocity=np.r_[0.,vel,0.]), dict(
            seconds=time.monotonic()-start,cells=self.cells,angles=self.angles,steps=len(self.t)-1,
            radius_cm=self.ext.R,luminosity_erg_s=self.ext.L,q=self.ext.q,scalar_scale=self.ext.scalar_scale,
            time_seconds=float(self.t[-1]*self.ext.R/C),surface_lapse_relative=self.ext.Nb/self.ext.old_Nb-1,
            remaining_gas_to_radiation_pressure=self.ext.Pgas/self.ext.Prad,
            maximum_normalized_response=np.max(abs(traces[:,1:]),axis=0).tolist(),
            maximum_physical_response=(np.max(abs(traces[:,1:]),axis=0)*abs(self.ext.scalar_scale)).tolist(),
            energy_work_relative_error=energy_error,maximum_field_at_r_ge_8=max(boundary),
            emitted_energy_erg=self.ext.L*self.ext.R/C*.75,
            body_mass_debit_g=-self.ext.L*self.ext.R/C**3*.75)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    bindings=[Path(__file__),Path(vacuum.__file__),Path(vacuum.q.exterior.__file__)]
    bindings += [BASE/'def-native-radiative-envelope'/p for p in ['fine.npz','thin.npz','result.json','audit.json']]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',symbolic=symbolic(),
        claim='Compute causal angular photon E,Pr and mass debit J from the native grey envelope and solve their coupled exterior scalar/metric response. Save the additive surface derivative needed by an interior boundary match.',
        input='Phase90 fine/thin native envelope endpoints. Compact diagnostic L(t)=L0*sin(pi*t/2)^4 for 0<t<2 in units R/c; not an astrophysically established drive.',
        boundary='Zero initial perturbations; dphi(1,t)=0 defines the additive exterior forcing. Outer r=9 is outside causal reach for t<=6. Homogeneous interior-dependent surface response remains to be matched.',
        approximation='First order in G*L/c^5 about exact vacuum scalar metric; outgoing isotropic hemisphere, collisionless exterior, fixed emitting radius. Finite gas pressure truncation retained as a physical limitation.',
        settings=[dict(name='coarse',envelope='fine',cells=320,angles=32,dt_fraction=.4),
                  dict(name='fine',envelope='fine',cells=640,angles=64,dt_fraction=.2),
                  dict(name='angular',envelope='fine',cells=640,angles=128,dt_fraction=.2),
                  dict(name='thin',envelope='thin',cells=640,angles=128,dt_fraction=.2)],
        gates=dict(trace_refinement_relative=.01,angular_relative=.002,envelope_cutoff_relative=.01,
                   energy_work_relative=1e-8,outer_field_relative=1e-6,conservation_relative=1e-10),
        budget=dict(total_seconds=240,pilot_cells=80,pilot_angles=16,pilot_end=1,production_paths=4,
                    native_EOS_calls=0,BLAS_threads=1,automatic_expansion=False),
        decision='If gates pass, use this exterior inhomogeneous response for the next whole-star boundary match. Never relabel the clamped-surface result as the final stellar charge or nonabsorption.',
        bindings={str(p.relative_to(ROOT)):digest(p) for p in bindings}))


def pilot():
    assert not (OUT/'pilot.json').exists()
    start=time.monotonic();signal.alarm(30)
    model=Response('fine',80,16,.4,end=1.)
    _,row=model.run()
    row['total_seconds']=time.monotonic()-start
    write(OUT/'pilot.json',row);signal.alarm(0)
    print(json.dumps(row),flush=True)


def compare(a,b):
    small=np.column_stack([np.interp(b[:,0],a[:,0],a[:,i]) for i in range(1,4)])
    return (np.max(abs(small-b[:,1:]),axis=0)/np.max(abs(b[:,1:]),axis=0)).tolist()


def run():
    assert not (OUT/'result.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert digest(ROOT/p)==sha,p
    budget=json.loads((OUT/'budget-review.json').read_text())
    assert budget['proceed']
    signal.alarm(210);start=time.monotonic();rows={};data={}
    for setting in plan['settings']:
        model=Response(setting['envelope'],setting['cells'],setting['angles'],setting['dt_fraction'])
        value,row=model.run();name=setting['name'];rows[name]=row;data[name]=value
        np.savez(OUT/f'{name}.npz',**value)
        write(OUT/'progress.json',rows)
        print('EXTERIOR',name,json.dumps(row),flush=True)
    errors=dict(refinement=compare(data['coarse']['traces'],data['fine']['traces']),
                angular=compare(data['fine']['traces'],data['angular']['traces']),
                envelope_cutoff=compare(data['angular']['traces'],data['thin']['traces']))
    passed=max(errors['refinement'])<.01 and max(errors['angular'])<.002 and max(errors['envelope_cutoff'])<.01
    passed=passed and all(r['maximum_field_at_r_ge_8']<1e-6*max(r['maximum_normalized_response']) for r in rows.values())
    result=dict(classification='Counterexample candidate',passed=bool(passed),rows=rows,errors=errors,
        seconds=time.monotonic()-start,causal_exterior_photon_metric_scalar_solved=True,
        interior_surface_value_matched=False,whole_star_inventory_luminosity_matched=False,
        final_stellar_charge_solved=False,astrophysical_drive_established=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0)
    print('FINAL',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','pilot','run'])
    globals()[parser.parse_args().action]()

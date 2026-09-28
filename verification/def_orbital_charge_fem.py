"""Physical Kepler harmonics on the accepted GR FEM, in the adiabatic sector.

Constant scalar lifting supplies the omitted boundary degree of freedom.
The mechanical contrast uses a momentum-only forcing and a direct impedance
difference, never subtraction of near-unit scalar solutions.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import jv
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_gr_interface_patch as patch
import def_radiative_exterior as exterior

go=patch.go;write=go.write
OUT=patch.OUT.parent/'def-orbital-charge-fem'
BENCH=exterior.s.OUT/'companion-benchmark.json'


def controls():
    d,c,z=sp.symbols('d c z',nonzero=True)
    assert sp.simplify(d/(c+z)-d/c+d*z/(c*(c+z)))==0
    x=sp.symbols('x',positive=True)
    j=sp.sin(x)/x
    assert sp.simplify(x*sp.diff(j,x)-(sp.I*x-1)*j-sp.exp(-sp.I*x))==0
    # Independent eccentric-anomaly integration, including dM=(1-e cos E)dE.
    b=json.loads(BENCH.read_text());e=b['parameters']['eccentricity'];rows=[]
    for n in [1,2,3]:
        exact=2*jv(n,n*e)
        direct=quad(lambda E:np.cos(n*(E-e*np.sin(E))),0,2*np.pi,epsabs=1e-13)[0]/np.pi
        assert abs(exact-direct)<2e-14
        rows.append(dict(harmonic=n,coefficient=float(exact),quadrature_error=abs(exact-direct)))
    return dict(classification='Proven',identities_passed=True,
        identities='a_free-a_clamped=-drive*dZ/(D_clamped*(D_clamped+dZ)); flat j0 affine outgoing boundary equals exp(-ikr). Kepler a/D=1+2 sum J_n(ne) cos(nM).',
        numerical_classification='Counterexample candidate',kepler_checks=rows)


class Problem:
    def __init__(self,degree=4,outer=2,quadrature=6):
        patch.install();self.p=go.Problem(degree=degree,outer=outer,quadrature=quadrature)
        m=self.m=self.p.model;self.outer=outer
        a=m.data[:,:4].reshape(-1,2,2).astype(np.longdouble)
        b=m.data[:,4:8].reshape(-1,2,2).astype(np.longdouble)
        c=m.data[:,8:12].reshape(-1,2,2).astype(np.longdouble)
        w=m.data[:,12:16].reshape(-1,2,2).astype(np.longdouble)
        weights=m.weights.astype(np.longdouble);v=[q.astype(np.longdouble) for q in m.V]
        cov=[q.astype(np.longdouble) for q in m.cov];lift=-a[:,:,1]
        self.kl=sum(cov[i].T@(weights*sum(b[:,i,j]*lift[:,j] for j in range(2)))+v[i].T@(weights*c[:,i,1]) for i in range(2))
        self.ml=sum(v[i].T@(weights*w[:,i,1]) for i in range(2))
        self.kll=np.sum(weights*(np.einsum('ni,nij,nj->n',lift,b,lift)+c[:,1,1]),dtype=np.longdouble)
        self.mll=np.sum(weights*w[:,1,1],dtype=np.longdouble)
        self.fluid=m.indices[:m.surface_index+1,0]
        self.scalar=m.indices[:-1,1]
        point=m.bg.sample(np.array([float(outer)]))
        self.mu=float(point['m'][0]/outer);self.flux=float(point['v'][0]*outer)
        self.F=float(point['N'][0]*np.sqrt(1-2*self.mu))
        geom=m.original.radiation.geometry;self.R=float(geom.R)
        self.delay=quad(lambda r:float(geom.metric(np.array([r]))[1][0]/geom.metric(np.array([r]))[0][0]),1,outer,epsabs=1e-12)[0]
        # The accepted free-surface module owns the exact same ADM normalization.
        self.mass=json.loads((patch.base.h.OUT/'absolute-shoot/result-0.001.json').read_text())['ADM_geom_m']
        period=json.loads(BENCH.read_text())['leading_drive']['period_seconds']
        self.omega=2*np.pi*self.R/(patch.base.h.gr.C*period)
        self.error=0.

    def solve(self,n):
        start=time.monotonic();m=self.m;p=self.p;omega=np.longdouble(n)*np.longdouble(self.omega)
        matrix=(p.K-float(omega**2)*p.M).tocsc()
        rhs=-(self.kl-omega**2*self.ml).astype(np.clongdouble)
        def factor(ids):
            A=matrix[ids,:][:,ids].tocsc();scale=np.sqrt(abs(A.diagonal()));D=diags(1/scale)
            lu=splu((D@A@D).astype(complex).tocsc(),permc_spec='NATURAL')
            K=p.Kx[ids,:][:,ids];M=p.Mx[ids,:][:,ids]
            def linear(b):
                u=(lu.solve(np.asarray(b/scale,complex))/scale).astype(np.clongdouble)
                for _ in range(3):u+=(lu.solve(np.asarray((b-K@u+omega**2*(M@u))/scale,complex))/scale).astype(np.clongdouble)
                err=float(np.max(abs(b-K@u+omega**2*(M@u))/(abs(b)+abs(K)@abs(u)+omega**2*(abs(M)@abs(u))+1e-100)))
                self.error=max(self.error,err);assert err<1e-9,err
                return u
            return linear
        clamped=np.zeros(m.size,np.clongdouble);clamped[self.scalar]=factor(self.scalar)(rhs[self.scalar])
        force=np.zeros(m.size,np.clongdouble)
        force[self.fluid]=(rhs-p.Kx@clamped+omega**2*(p.Mx@clamped))[self.fluid]
        delta=factor(np.arange(m.size))(force)
        rc=-rhs@clamped+self.kll-omega**2*self.mll
        dr=-rhs@delta
        wave=exterior.outgoing(self.mu,self.flux,float(omega*self.outer))
        h=wave['h'];Z=wave['impedance'];phase=np.exp(-1j*complex(omega)*(1+self.delay))
        drive=phase/(self.F*h);Dc=rc/(self.outer*self.F)-Z;dZ=dr/(self.outer*self.F)
        ac=drive/Dc;af=drive/(Dc+dZ);da=-drive*dZ/(Dc*(Dc+dZ))
        outgoing=self.outer*da/h*np.exp(-1j*complex(omega)*self.delay)
        charge=-self.R/self.mass*outgoing
        # Independent endpoint derivative checks the weak reaction sign/units.
        V,D=m.evaluation(np.array([float(self.outer)]))
        derivative=self.outer*(D[1]@clamped)[0]
        reaction_check=float(abs(derivative-rc/(self.outer*self.F))/max(abs(rc/(self.outer*self.F)),1e-30))
        delta_native=af*(m.nativeV[0]@delta)+da*(m.nativeV[0]@clamped)
        pair=lambda q:[float(q.real),float(q.imag)]
        return dict(harmonic=n,omega_R_over_c=float(omega),charge=pair(charge),outgoing=pair(outgoing),
            Z_clamped=pair(rc/(self.outer*self.F)),Z_difference=pair(dZ),
            boundary_amplitude=pair(af),clamped_boundary_amplitude=pair(ac),
            derivative_reaction_relative=reaction_check,linear_residual=self.error,
            radiative_current_error=wave['current_relative_error'],
            native_displacement_RMS=float(np.sqrt(np.sum(m.original.weights*abs(delta_native)**2))),
            seconds=time.monotonic()-start)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Resolve outgoing material-motion response of the existing finite-temperature adiabatic GR benchmark and apply actual declared Kepler harmonics.',
        scope='Frozen mechanical background, linear l=0 regular incident scalar. No thermal stationarity assumption is certified. Companion alpha=beta*phi0 remains a leading weak-body input. Clamped comparison requires support force and is not the same-inventory static comparator.',
        method='Accepted hierarchical GR stiffness/inertia; constant scalar lifting; variational boundary reaction; direct fluid-row contrast; outgoing impedance and common surface-referenced tortoise phase. Static and harmonics1,2,3, degree2/4; outer3 and quadrature8 conditional.',
        gates=dict(spatial_relative=.02,outer_relative=.002,quadrature_relative=.002,linear_residual=1e-9,current=2e-8,endpoint_reaction=.02),
        decision='Pilot p4 harmonic1 reused in production. Stop after gate failure, no automatic degree/frequency expansion. Static-subtracted differences require their own spatial comparison before interpretation.',
        budget=dict(pilot_seconds=120,total_compute_seconds=600,CPU_threads=1,memory_GB=4,new_EOS_calls=0,time_steps=0),
        bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(patch.__file__),Path(go.__file__),Path(go.space.__file__),Path(go.space.task.__file__),Path(exterior.__file__),BENCH,patch.OUT/'p4-result.json']}))
    write(OUT/'symbolic.json',controls());signal.alarm(120);start=time.monotonic();p=Problem();setup=time.monotonic()-start
    row=p.solve(1);forecast=1.5*(4*setup+16*row['seconds']+20)
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',setup_seconds=setup,row=row,forecast_seconds=forecast,
        seconds=time.monotonic()-start,R_m=p.R,ADM_geom_m=p.mass,dofs=p.m.size,
        assumptions='Four assemblies and sixteen harmonics scaled from measured p4; smaller degrees/outer3/quadrature8 are unmeasured. 50 percent margin plus20s.'))
    print('PILOT',json.dumps(row),'FORECAST',forecast,flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert go.task.digest(Path(p))==h,p
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['forecast_seconds']<600
    assert not (OUT/'result.json').exists();signal.alarm(int(600-pilot['seconds']))
    resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));started=time.monotonic();cases={}
    for label,degree,outer,nq in [('p4',4,2,6),('p2',2,2,6),('outer',4,3,6),('quadrature',4,2,8)]:
        p=Problem(degree,outer,nq)
        rows=[pilot['row'] if label=='p4' and n==1 else p.solve(n) for n in range(4)]
        cases[label]=rows;write(OUT/(label+'.json'),dict(classification='Counterexample candidate',rows=rows))
        print('CASE',label,[r['charge'] for r in rows],flush=True)
        del p
        if label=='p2':
            differences=[abs(complex(*a['charge'])-complex(*b['charge']))/abs(complex(*b['charge'])) for a,b in zip(rows,cases['p4'])]
            if max(differences)>=.02:break
    benchmark=json.loads(BENCH.read_text());ecc=benchmark['parameters']['eccentricity'];par=benchmark['parameters'];inp=benchmark['inputs']
    amplitude=-benchmark['leading_drive']['alpha_companion']*inp['material_mass_geom_cm']/(par['semimajor_axis_over_star_radius']*inp['radius_cm'])
    static=complex(*cases['p4'][0]['charge']);rows=[]
    for n in range(4):
        response=complex(*cases['p4'][n]['charge']);drive=2*amplitude*jv(n,n*ecc) if n else amplitude
        row=dict(harmonic=n,charge_gain=[response.real,response.imag],drive_amplitude=drive,
            delta_alpha_over_phi0=abs(drive*response/par['background_phi']),
            frozen_static_subtracted_over_phi0=abs(drive*(response-static)/par['background_phi']))
        for name in ['p2','outer','quadrature']:
            if name in cases:
                value=complex(*cases[name][n]['charge']);row[name+'_relative']=abs(value-response)/max(abs(response),1e-100)
                if n:row[name+'_static_subtracted_relative']=abs((value-complex(*cases[name][0]['charge']))-(response-static))/max(abs(response-static),1e-100)
        rows.append(row)
    passed=len(cases)==4 and all(r.get('p2_relative',1)<.02 and r.get('outer_relative',1)<.002 and r.get('quadrature_relative',1)<.002 for r in rows)
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,cases=cases,seconds=time.monotonic()-started,
        actual_declared_kepler_harmonics_applied=True,thermal_background_evolved=False,
        full_dynamic_charge_solved=False,same_inventory_static_comparator_solved=False,observational_nuisance_applied=False)
    write(OUT/'result.json',result);print('RESULT',passed,rows,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()

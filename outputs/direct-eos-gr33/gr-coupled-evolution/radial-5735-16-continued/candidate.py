"""Counterexample candidate: executable spherical fluid/heat/metric evolution.

The native EOS is evaluated on every Runge--Kutta stage.  The material time
matrix is the differential primitive inverse; composition is advected, so no
uncomputed composition derivative is dropped from the material EOS derivative.
The outer ghost state is a computational outflow boundary, not an atmosphere.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import ast
import hashlib
import json
import subprocess
import time

import numpy as np
import sympy as sp
import direct_eos_gr as g
import opacity_tables

ROOT = g.ROOT
ld = np.longdouble
C = ld(g.c.gr.C) * 100
GRAV = ld(g.c.gr.G) * 1000 / C**4
TAU = ld('0.0004113000088929766')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def worker_init():
    global EOS, OPACITY
    EOS = g.EOS()
    OPACITY = opacity_tables.Opacity()


def material(row):
    lr, lt, x = row
    a = EOS(2, float(lr), float(lt), np.asarray(x, float))
    nuclei = x / g.c.A
    hydrogen, helium = x[g.c.Z == 1].sum(), x[g.c.Z == 2].sum()
    par = [(nuclei @ g.c.Z) / nuclei.sum(), hydrogen,
           np.clip(1-hydrogen-helium,0,.1), lr / np.log(10.), lt / np.log(10.)]
    op = OPACITY(np.asarray(par, float))
    return np.r_[a, op]


def material_key(row):
    """Cache only when every binary64 argument actually sent to EOS/opacity agrees."""
    lr,lt,x=row
    nuclei=x/g.c.A
    hydrogen,helium=x[g.c.Z==1].sum(),x[g.c.Z==2].sum()
    par=np.asarray([(nuclei@g.c.Z)/nuclei.sum(),hydrogen,np.clip(1-hydrogen-helium,0,.1),
                    lr/np.log(10.),lt/np.log(10.)],float)
    return (float(lr),float(lt),*np.asarray(x,float),*par)


def symbolic():
    v, e, p, q, er, et, pr, pt = sp.symbols('v e p q er et pr pt')
    W2 = 1/(1-v*v)
    E = (e+p*v*v+2*q*v)*W2
    S = ((e+p)*v+q*(1+v*v))*W2
    assert sp.cancel(sp.diff(E, v)-2*W2**2*(v*(e+p)+q*(1+v*v))) == 0
    assert sp.cancel(sp.diff(S, v)-W2**2*((1+v*v)*(e+p)+4*q*v)) == 0
    assert sp.cancel(sp.diff(E, e)*er+sp.diff(E, p)*pr-W2*(er+pr*v*v)) == 0
    assert sp.cancel(sp.diff(S, e)*et+sp.diff(S, p)*pt-W2*v*(et+pt)) == 0
    # Project the covariant energy/momentum balances before eliminating baryon
    # compression.  This checks the finite-v heat terms, including 2 Q a_s.
    a,N,Nr,at,r,c=sp.symbols('a N Nr at r c', nonzero=True)
    ev,pv,qv,vv,ex,px,qx,vx=sp.symbols('ev pv qv vv ex px qx vx')
    variables=[e,p,q,v]
    dt=lambda f:sum(sp.diff(f,x)*d for x,d in zip(variables,[ev,pv,qv,vv]))
    dr=lambda f:sum(sp.diff(f,x)*d for x,d in zip(variables,[ex,px,qx,vx]))
    R=((e+p)*v*v+p*(1-v*v)+2*q*v)*W2
    eb=dt(E)+at/a*(E+R)+c*N/a*dr(S)+2*c*Nr/a*S+2*c*N/(a*r)*S
    sb=dt(S)+2*at/a*S+c*N/a*dr(R)+c*Nr/a*(R+E)+2*c*N/(a*r)*(R-p)
    theta=W2*v*vv+c*N/a*W2*vx+at/a+c*v*Nr/a+2*c*N*v/(a*r)
    acceleration=W2*(vv+c*N*v/a*vx)+c*Nr/a+v*at/a
    internal=ev+c*N*v/a*ex+(e+p)*theta+v*qv+c*N/a*qx+2*q*acceleration+2*c*N*q/(a*r)
    assert sp.cancel(eb-v*sb-internal)==0
    # Material composition is constant: D_u X=0.  With baryon conservation,
    # theta + D_u ln(tau/(K*T^2)) = -kr D_u ln(rho)-kt D_u ln(T).
    n, t, kn, kt = sp.symbols('n t kn kt')
    assert sp.expand(-n-kn*n-(kt+2)*t+(kn+1)*n+(kt+2)*t) == 0
    return dict(classification='Proven', passed=True,
                scope='Conditional algebra for the declared smooth EOS and constant proper tau; no continuum error bound.')


class Star:
    def __init__(self, cells, pool):
        self.pool = pool
        self.material_cache = {}
        self.source = g.OUT/'initial-state-17-4.npz'
        s = dict(np.load(self.source))
        count = len(s['lnd'])
        edges = np.unique(np.linspace(0, count, min(cells, count)+1).astype(int))
        # Reverse the original outer-to-inner ordering.  Coarsening is explicit;
        # it is not a mass refit or a claim to preserve all original cell masses.
        faces = s['radius_faces_m'][::-1].astype(ld)*100
        self.rf = faces[edges]
        assert self.rf[0] == 0 and np.all(np.diff(self.rf) > 0)
        r_all = s['r_mid_m'][::-1].astype(ld)*100
        desired = np.cbrt((self.rf[1:]**3+self.rf[:-1]**3)/2)
        chosen = np.array([lo+np.argmin(abs(r_all[lo:hi]-r))
                           for lo, hi, r in zip(edges[:-1], edges[1:], desired)])
        self.indices = count-1-chosen
        self.r = r_all[chosen]
        self.volume = 4*ld(str(np.pi))/3*np.diff(self.rf**3)
        self.area_difference_over_volume = 4*np.pi*np.diff(self.rf**2)/self.volume
        self.n = len(chosen)
        self.base = np.column_stack([s['lnd'][self.indices], s['lnT'][self.indices],
                                     np.zeros((self.n, 2)), s['X'][self.indices]]).astype(ld)
        aux = np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos'][self.indices]
        rho = np.exp(self.base[:, 0])
        rest = (self.base[:, 4:]/g.c.A) @ g.c.W * C*C
        self.qscale = rho*(rest+aux[:, 2])+aux[:, 1]
        lum = np.r_[0., np.load(g.OUT/'gr-transport/diagnostics.npz')['interior_Linf'], 0.]
        Q = (lum[self.indices]+lum[self.indices+1])/2/(4*np.pi*self.r**2*np.exp(2*s['nu'][self.indices])*C)
        self.base[:, 3] = Q/self.qscale
        self.original_baryons = ld(s['dm'].sum())  # source dm is already grams
        self.original_mass = ld(s['mass_faces_geom'][0])*100

    def gradient(self, value, odd=False):
        answer = np.gradient(value, self.r, axis=0, edge_order=2)
        # Centre parity uses the regular extension through r=0.
        answer[0] = value[0]/self.r[0] if odd else 0
        return answer

    def faces(self, value, odd=False):
        # One shared flux at every internal face; no cell-local duplicate flux.
        t = (self.rf[1:-1]-self.r[:-1])/np.diff(self.r)
        result = np.empty(self.n+1, dtype=ld)
        result[1:-1] = value[:-1]+t*(value[1:]-value[:-1])
        result[0] = 0 if odd else value[0]
        result[-1] = value[-1]  # specified zero-gradient computational ghost
        return result

    def divergence(self, flux):
        return np.diff(4*np.pi*self.rf**2*flux)/self.volume

    def state(self, y):
        assert np.all(np.isfinite(y)) and np.max(abs(y[:, 2])) < .2
        x = y[:, 4:]
        assert x.min() >= -1e-15 and np.max(abs(x.sum(1)-1)) < 1e-10
        # Native Fortran state is isolated in worker processes.
        inputs=list(zip(y[:,0],y[:,1],x))
        keys=[material_key(row) for row in inputs]
        missing={key:row for key,row in zip(keys,inputs) if key not in self.material_cache}
        values=self.pool.map(material,missing.values(),chunksize=4)
        self.material_cache.update(zip(missing,values))
        rows=[self.material_cache[key] for key in keys]
        if len(self.material_cache)>8*self.n:
            self.material_cache={key:row for key,row in zip(keys,rows)}
        aux = np.asarray(rows, dtype=ld)
        rho, T = np.exp(y[:, 0]), np.exp(y[:, 1])
        P, u, cvT = aux[:, 1], aux[:, 2], aux[:, 10]
        assert np.all(cvT > 0) and np.all(P > 0)
        rest = (x/g.c.A) @ g.c.W * C*C
        eps = rho*(rest+u)
        w = eps+P
        v, Q = y[:, 2], y[:, 3]*self.qscale
        W = 1/np.sqrt(1-v*v)
        D = rho*W
        E = (eps+P*v*v+2*Q*v)*W*W
        S = (w*v+Q*(1+v*v))*W*W
        R = (eps*v*v+P+2*Q*v)*W*W
        mf = np.r_[ld(0), np.cumsum(GRAV*E*self.volume)]
        m = mf[:-1]+GRAV*E*4*np.pi/3*(self.r**3-self.rf[:-1]**3)
        f = 1-2*m/self.r
        assert f.min() > 0 and 1-2*mf[-1]/self.rf[-1] > 0
        a = 1/np.sqrt(f)
        nur = a*a*(m/self.r**2+4*np.pi*GRAV*self.r*R)
        # Outer normalization is a coordinate-time gauge.  It does not match
        # a finite-pressure, heat-carrying material boundary to Schwarzschild.
        nu = np.empty(self.n, dtype=ld)
        nu[-1] = .5*np.log1p(-2*mf[-1]/self.rf[-1])-nur[-1]*(self.rf[-1]-self.r[-1])
        nu[:-1] = nu[-1]-np.cumsum((.5*(nur[1:]+nur[:-1])*np.diff(self.r))[::-1])[::-1]
        N = np.exp(nu)
        at = -4*np.pi*GRAV*C*self.r*N*a*a*S
        K = 16*ld('5.670400e-5')*T**3/(3*rho*aux[:, 21])
        return dict(aux=aux, rho=rho, T=T, P=P, u=u, rest=rest, w=w, v=v, Q=Q,
                    W=W, D=D, E=E, S=S, R=R, mf=mf, m=m, a=a, N=N, nur=nur, at=at, K=K)

    def rhs(self, delta):
        y = self.base+delta
        z = self.state(y)
        rho,T,P,w,v,Q,W,D,E,S,R,a,N,nur,at,K = [z[k] for k in
            ['rho','T','P','w','v','Q','W','D','E','S','R','a','N','nur','at','K']]
        aux = z['aux']
        speed = C*N*v/a
        fB = self.faces(N*D*v, odd=True)
        fE = self.faces(N*S/a, odd=True)
        fS = self.faces(N*R)
        Bdot = -C*self.divergence(fB)-at*D
        Edot = -C*self.divergence(fE)
        Sdot = -C*self.divergence(fS)-at*S-N*nur*C*E+C*N*P*self.area_difference_over_volume
        # Material derivatives remove composition changes from EOS time rows.
        right = np.column_stack([Bdot+a*speed*self.gradient(D),
                                 np.zeros(self.n),
                                 Sdot+a*speed*self.gradient(S, odd=True), np.zeros(self.n)]).astype(ld)
        er = rho*(z['rest']+z['u']+aux[:, 9])
        et = rho*aux[:, 10]
        pr, pt = P*aux[:, 5], P*aux[:, 6]
        matrix = np.zeros((self.n, 4, 4), dtype=ld)
        matrix[:, 0, 0] = a*D
        matrix[:, 0, 2] = a*rho*W**3*v
        # Project total energy into the material frame analytically.  Independent
        # discrete rest-energy and baryon fluxes otherwise create false heat.
        # This is a thermodynamic MOL formulation: total-energy conservation is
        # an independently measured budget, not imposed by a post-step rescale.
        matrix[:,1] = np.column_stack([rho*aux[:,9]-P, et, 2*Q*W*W, v*self.qscale])
        right[:,1] = -C*N/(a*W*W)*self.gradient(Q,odd=True)-2*Q*(C*N*nur/a+v*at/a+C*N/(a*self.r))
        matrix[:, 2] = a[:,None]*np.column_stack([W*W*v*(er+pr), W*W*v*(et+pt),
                           W**4*((1+v*v)*w+4*Q*v), (1+v*v)*W*W*self.qscale])
        h = K*T/C**2
        kr, kt = -aux[:,22], 5-aux[:,23]
        matrix[:,3] = (W/N)[:,None]*np.column_stack([-TAU*Q*kr/2,
                            h*v-TAU*Q*kt/2, h*W*W, TAU*self.qscale])
        right[:,3] = -Q-h*(C*self.gradient(y[:,1])/(a*W)+C*W*nur/a+W*v*at/(N*a))
        scales = np.max(abs(matrix),axis=2)
        A = np.asarray(matrix/scales[:,:,None],float)
        b = np.asarray(right/scales,float)
        rates = np.linalg.solve(A,b[:,:,None])[:,:,0].astype(ld)
        residual = float(np.max(abs(np.einsum('nij,nj->ni', A, rates)-b)/(1+abs(b))))
        assert residual < 1e-10, residual
        out = np.zeros_like(y)
        for j in range(4):
            if j == 3:
                # Q, rather than Q/w0(r), is the advected physical field.
                spatial = self.gradient(Q,odd=True)/self.qscale
            else:
                spatial = self.gradient(y[:,j],odd=(j==2))
            out[:,j] = rates[:,j]-speed*spatial
        # Upwind passive composition; no reaction/frozen-composition shortcut.
        left = np.vstack([y[0,4:]*0, np.diff(y[:,4:],axis=0)/np.diff(self.r)[:,None]])
        right_x = np.vstack([np.diff(y[:,4:],axis=0)/np.diff(self.r)[:,None],y[-1,4:]*0])
        out[:,4:] = -speed[:,None]*np.where(speed[:,None]>=0,left,right_x)
        z['time_matrix_residual'] = residual
        z['boundary_rates'] = np.array([-C*4*np.pi*self.rf[-1]**2*fB[-1],
                                       -C*4*np.pi*self.rf[-1]**2*fE[-1]])
        return out,z


def run(args):
    folder = g.OUT/'gr-coupled-evolution'/args.name
    assert not folder.exists(), 'Use a new name; keep failed paths immutable.'
    folder.mkdir(parents=True)
    candidate=folder/'candidate.py'
    candidate.write_bytes(Path(__file__).read_bytes())
    inputs = [candidate, g.OUT/'initial-state-17-4.npz',
              g.OUT/'gr-microphysics/auxiliaries.npz', g.OUT/'gr-transport/diagnostics.npz',
              ROOT/'verification/direct_eos_gr.py',ROOT/'verification/direct_ion_eos.py',
              ROOT/'verification/opacity_tables.py',ROOT/'verification/opacity_cubic.py',
              opacity_tables.o.OUT/'tables-enriched-tables.npz']
    plan = dict(classification='Counterexample candidate', checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        bindings={str(p.relative_to(ROOT)):digest(p) for p in inputs}, requested_cells=args.cells,
        steps=args.steps,proper_tau_seconds=float(TAU),duration_seconds=float(TAU*args.duration_tau),
        method='SSPRK2 on primitive increments; solve the full material baryon/internal-energy/momentum/causal-heat time matrix at each stage. The energy projection is symbolically equivalent before space discretization. Fresh native EOS, composition advection, common baryon/momentum face fluxes, radial energy-mass constraint and polar lapse at every stage. Total-energy budget is measured independently; it is not exactly conserved by this primitive MOL scheme.',
        boundary='Regular centre; zero-gradient computational outer ghost. Finite photospheric pressure retained. No insulating Q=0 enforcement after t=0 and no physical atmosphere/exterior match.',
        limitations=['Space discretization and primitive chain-rule defects are measured, not certified.',
                     'Coarse grids change initial quadrature masses; never renormalize or claim the original mass match.',
                     'No reactions on this short heat/fluid path; nuclear evolution, physical EOS, calibrated heat time and observations remain open.'],
        symbolic=symbolic())
    resume=None
    if args.resume:
        parent=g.OUT/'gr-coupled-evolution'/args.resume
        receipt=json.loads((parent/'checkpoint-stop.json').read_text())
        assert receipt['stopped']
        for rel,h in receipt['sha256'].items(): assert digest(parent/rel)==h,rel
        before=json.loads((parent/'plan.json').read_text())
        for key in ['requested_cells','steps','proper_tau_seconds','duration_seconds']:
            assert plan[key]==before[key],key
        # Changes to orchestration cannot silently change an ongoing equation.
        def science(path):
            tree=ast.parse(path.read_text())
            return [ast.dump(node,include_attributes=False) for node in tree.body
                    if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name in
                    ['Star','material','material_key','worker_init','symbolic']]
        assert science(candidate)==science(parent/'candidate.py')
        rows=json.loads((parent/'progress.json').read_text())['rows']
        start=rows[-1]['step']
        resume=dict(parent=parent,start=start,history=rows[:-1],state=dict(np.load(parent/f'step-{start:04d}.npz')),
                    initial=np.load(parent/'step-0000.npz')['totals'])
        plan['resume']=dict(parent=args.resume,step=start,receipt_sha256=digest(parent/'checkpoint-stop.json'),
                            science_AST_identical=True,checkpoint_replay_required_bitwise=True)
        for step in range(start):
            (folder/f'step-{step:04d}.npz').write_bytes((parent/f'step-{step:04d}.npz').read_bytes())
    write(folder/'plan.json',plan)
    started=time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers,initializer=worker_init) as pool:
        star=Star(args.cells,pool)
        delta=np.zeros_like(star.base) if resume is None else resume['state']['delta'].copy()
        dt=TAU*args.duration_tau/args.steps
        accumulated=np.zeros(2,dtype=ld) if resume is None else resume['state']['boundary_exchange'].copy()
        history=[] if resume is None else resume['history']
        initial=None if resume is None else resume['initial'].copy()
        start=0 if resume is None else resume['start']
        for step in range(start,args.steps+1):
            k1,z=star.rhs(delta)
            totals=np.array([np.sum(z['a']*z['D']*star.volume),np.sum(z['E']*star.volume)],dtype=ld)
            if resume is not None and step==start:
                for key in ['m','a','N','Q']:
                    assert np.array_equal(z[key],resume['state'][key]),('Checkpoint replay differs',key)
                assert np.array_equal(totals,resume['state']['totals'])
                write(folder/'checkpoint-replay.json',dict(passed=True,bitwise_fields=['m','a','N','Q','totals'],step=start))
                print('COUPLED CHECKPOINT REPLAY BITWISE',start,flush=True)
            if initial is None: initial=totals.copy()
            row=dict(step=step,time_seconds=float(step*dt),maximum_changes=np.max(abs(delta[:,:4]),axis=0).astype(float).tolist(),
                relative_budget_defect=np.asarray((totals-initial-accumulated)/initial,float).tolist(),
                maximum_time_matrix_residual=z['time_matrix_residual'],maximum_velocity=float(np.max(abs(z['v']))),
                outer_mass_cm=float(z['mf'][-1]),elapsed_seconds=time.monotonic()-started)
            history.append(row)
            np.savez_compressed(folder/f'step-{step:04d}.npz',delta=delta,m=z['m'],a=z['a'],N=z['N'],Q=z['Q'],
                                totals=totals,boundary_exchange=accumulated,indices=star.indices,radius_cm=star.r)
            write(folder/'progress.json',dict(classification='Counterexample candidate',rows=history))
            print('COUPLED STEP',args.name,step,'/',args.steps,json.dumps(row),flush=True)
            if step==args.steps: break
            k2,z2=star.rhs(delta+dt*k1)
            delta += dt*(k1+k2)/2
            accumulated += dt*(z['boundary_rates']+z2['boundary_rates'])/2
        changes=np.max(abs(delta[:,:4]),axis=0)
        assert np.all(changes>0), changes
        result=dict(classification='Counterexample candidate',completed=True,actual_nonlinear_finite_time_path=True,
                    cells=star.n,steps=args.steps,initial_baryon_relative_to_original=float(initial[0]/star.original_baryons-1),
                    initial_mass_relative_to_original=float(GRAV*initial[1]/star.original_mass-1),last=history[-1],
                    physical_EOS_certified=False,full_physical_GR_stellar_evolution=False,observational_closure=False)
        write(folder/'result.json',result)
    write(folder/'manifest.json',dict(sha256={p.name:digest(p) for p in folder.iterdir() if p.is_file()}))
    print('COUPLED FINITE PATH COMPLETE',args.name,flush=True)


def comparison_plan():
    folder=g.OUT/'gr-coupled-evolution'
    path=folder/'time-comparison-plan.json'
    assert not path.exists()
    write(path,dict(classification='Counterexample candidate',runs=['radial-5735-4','radial-5735-8','radial-5735-16'],
        fields=['delta_lnrho','delta_lnT','v_over_c','delta_Q_over_initial_enthalpy'],
        gate='All four maximum endpoint differences must decrease on halving dt; lnT and Q must have observed order at least 1.5. Report differences directly, without claiming a rigorous continuous-time error bound.',
        minimum_thermal_heat_order=1.5,
        restrictions='Same 5735 spatial cells and end time. No mass rescaling, omitted failed state or modified EOS/heat parameter. The 8-step runner has no memoization and an incorrect grams-to-kilograms factor only in its final original-baryon diagnostic; recompute that diagnostic here from the actual gram-valued input. The 4/16 runners memoize only identical native-EOS and opacity arguments. Evolution equations are unchanged.',
        physical_EOS_certified=False,observational_closure=False))


def compare():
    folder=g.OUT/'gr-coupled-evolution'
    plan=json.loads((folder/'time-comparison-plan.json').read_text())
    assert not (folder/'time-comparison.json').exists()
    arrays=[];records=[];bindings={};resolved=[]
    for name in plan['runs']:
        runpath=folder/name
        if (runpath/'continuation.json').exists():
            link=json.loads((runpath/'continuation.json').read_text())
            bindings[str((runpath/'continuation.json').relative_to(ROOT))]=digest(runpath/'continuation.json')
            continued=folder/link['name']
            assert json.loads((continued/'checkpoint-replay.json').read_text())['passed']
            assert json.loads((continued/'plan.json').read_text())['resume']['receipt_sha256']==digest(runpath/'checkpoint-stop.json')
            runpath=continued
        resolved.append(runpath)
        rp=json.loads((runpath/'plan.json').read_text())
        result=json.loads((runpath/'result.json').read_text())
        assert result['completed'] and result['cells']==5735
        for rel,h in rp['bindings'].items(): assert digest(ROOT/rel)==h,rel
        for rel,h in json.loads((runpath/'manifest.json').read_text())['sha256'].items(): assert digest(runpath/rel)==h,rel
        initial=np.load(runpath/'step-0000.npz')
        final=np.load(runpath/f"step-{result['steps']:04d}.npz")
        arrays.append(final['delta'])
        state=np.load(g.OUT/'initial-state-17-4.npz')
        aux=np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
        volume=-4*np.pi/3*np.diff((state['radius_faces_m'].astype(ld)*100)**3)
        capacity=np.sum(np.exp(state['lnd'].astype(ld))*aux[:,10]*volume)
        error=final['totals']-initial['totals']-final['boundary_exchange']
        record=dict(name=name,steps=result['steps'],duration_seconds=rp['duration_seconds'],
            corrected_initial_baryon_relative_to_original=float(initial['totals'][0]/np.sum(state['dm'].astype(ld))-1),
            initial_mass_relative_to_original=result['initial_mass_relative_to_original'],
            relative_budget_defect=np.asarray(error/initial['totals'],float).tolist(),
            energy_budget_defect_over_initial_thermal_capacity=float(error[1]/capacity),
            maximum_changes=np.max(abs(final['delta'][:,:4]),axis=0).astype(float).tolist(),
            maximum_composition_change=float(np.max(abs(final['delta'][:,4:]))),
            maximum_metric_changes={k:float(np.max(abs(final[k]-initial[k]))) for k in ['m','a','N']})
        records.append(record)
        bindings[str((runpath/'manifest.json').relative_to(ROOT))]=digest(runpath/'manifest.json')
    assert [r['steps'] for r in records]==[4,8,16]
    assert len(set(r['duration_seconds'] for r in records))==1
    for k in ['indices','radius_cm','totals','m','a','N']:
        baseline=np.load(resolved[0]/'step-0000.npz')[k]
        assert all(np.array_equal(np.load(path/'step-0000.npz')[k],baseline) for path in resolved),k
    errors=np.array([np.max(abs(b[:,:4]-a[:,:4]),axis=0) for a,b in zip(arrays[:-1],arrays[1:])])
    assert np.all(errors>0)
    orders=np.log2(errors[0]/errors[1])
    passed=bool(np.all(errors[1]<errors[0]) and np.min(orders[[1,3]])>=plan['minimum_thermal_heat_order'])
    write(folder/'time-comparison.json',dict(classification='Counterexample candidate',passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(),observed_orders=orders.astype(float).tolist(),runs=records,
        physical_EOS_certified=False,rigorous_time_error_bound=False,physical_exterior_match=False,observational_closure=False))
    for name in ['time-comparison-plan.json','time-comparison.json']:
        bindings[str((folder/name).relative_to(ROOT))]=digest(folder/name)
    write(folder/'time-comparison-manifest.json',dict(sha256=bindings))
    print('COUPLED TIME COMPARISON',passed,'orders',orders,'differences',errors,flush=True)
    assert passed,'Keep the failed time comparison; do not relabel it convergence.'


def verify():
    folder=g.OUT/'gr-coupled-evolution'
    for rel,h in json.loads((folder/'time-comparison-manifest.json').read_text())['sha256'].items():
        assert digest(ROOT/rel)==h,rel
        if rel.endswith('/manifest.json'):
            parent=(ROOT/rel).parent
            for name,binding in json.loads((ROOT/rel).read_text())['sha256'].items(): assert digest(parent/name)==binding,name
    result=json.loads((folder/'time-comparison.json').read_text())
    assert result['passed'] and not result['rigorous_time_error_bound'] and not result['physical_exterior_match']
    assert all(r['maximum_metric_changes'][k]>0 for r in result['runs'] for k in ['m','a','N'])
    assert symbolic()['passed']
    print('PASS full-grid coupled paths and specified finite time refinement; physical closure remains open',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('name')
    parser.add_argument('--cells',type=int,default=96)
    parser.add_argument('--steps',type=int,default=8)
    parser.add_argument('--duration-tau',type=float,default=1.)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--resume')
    args=parser.parse_args()
    if args.name in ['comparison_plan','compare','verify']: globals()[args.name]()
    else: run(args)

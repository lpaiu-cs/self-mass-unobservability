"""Counterexample candidate: native conservative two-carrier GR evolution.

Baryon, coordinate energy, momentum and isotope balances use shared faces.
Variable-step BDF2 follows an initial backward-Euler step. The heat equations
retain the full time/space derivative of their native entropy coefficient.
The LTE two-carrier model and reflecting wall remain conditional physics.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import subprocess
import time

import numpy as np
import sympy as sp
from scipy.linalg import solve_banded
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

import gr_covariant_two_carrier as parent

e, ld = parent.e, parent.ld
OUT = e.g.OUT/'gr-conservative-two-carrier'
SCALE = parent.parent.SCALE
ATOL = np.array([1e-13, 1e-9, 1e-14, 1e-17, 1e-17]+[1e-13]*26, dtype=ld)


def weights(h, previous_h):
    if previous_h is None:
        return np.array([1, -1, 0], dtype=ld)
    ratio = h/previous_h
    assert 0 < ratio <= ld('2.000000000001'), ratio
    return np.array([(1+2*ratio)/(1+ratio), -(1+ratio), ratio*ratio/(1+ratio)], dtype=ld)


def upwind(value, speed):
    take_left = speed[1:-1] >= 0
    if value.ndim == 2:
        take_left = take_left[:, None]
    return np.concatenate([value[:1], np.where(take_left, value[:-1], value[1:]), value[-1:]])


def fluxes(star, z):
    velocity = star.faces(z['N']*z['v']/z['a'], odd=True)
    baryon = velocity*upwind(z['B'], velocity)
    energy = velocity*upwind(z['E'], velocity)+star.faces(z['N']/z['a']*(z['P']*z['v']+z['Q']), odd=True)
    momentum = velocity*upwind(z['AS'], velocity)+star.faces(z['N']*(z['P']+z['Q']*z['v']))
    momentum[-1] = z['N'][-1]*z['R'][-1]  # Retain the parent's even wall stress.
    return baryon, energy, momentum


class ConservativeStar(parent.CovariantStar):
    def evaluate(self, delta):
        z = super().state(self.base+delta)
        zero = self.reference
        rho = zero['rho']*np.exp(delta[:, 0])
        T = zero['T']*np.exp(delta[:, 1])
        drest = (delta[:, 5:]/e.g.c.A)@e.g.c.W*e.C**2
        rest = zero['rest']+drest
        P, u, v, Q = z['P'], z['u'], delta[:, 2], z['Q']
        W = 1/np.sqrt(1-v*v)
        eps0 = zero['rho']*(zero['rest']+zero['u'])
        deps = zero['rho']*((zero['rest']+zero['u'])*np.expm1(delta[:, 0])
            +np.exp(delta[:, 0])*(drest+u-zero['u']))
        dE = (deps+(P+eps0)*v*v+2*Q*v)/(1-v*v)
        dmf = np.r_[ld(0), np.cumsum(e.GRAV*dE*self.volume)]
        dm = dmf[:-1]+e.GRAV*dE*self.partial_volume
        factor = np.sqrt(1-2*zero['a']**2*dm/self.r)
        assert np.all(factor > 0)
        a = zero['a']/factor
        da = zero['a']*(2*zero['a']**2*dm/self.r)/(factor*(1+factor))
        mf, m = zero['mf']+dmf, zero['m']+dm
        eps, w = rho*(rest+u), rho*(rest+u)+P
        S = (w*v+Q*(1+v*v))*W*W
        R = (eps*v*v+P+2*Q*v)*W*W
        nur = a*a*(m/self.r**2+4*np.pi*e.GRAV*self.r*R)
        nu = np.empty(self.n, dtype=ld)
        assert 1-2*mf[-1]/self.rf[-1] > 0
        nu[-1] = .5*np.log1p(-2*mf[-1]/self.rf[-1])-nur[-1]*(self.rf[-1]-self.r[-1])
        nu[:-1] = nu[-1]-np.cumsum((.5*(nur[1:]+nur[:-1])*np.diff(self.r))[::-1])[::-1]
        z.update(rho=rho, T=T, rest=rest, w=w, v=v, W=W, D=rho*W, E=eps0+dE,
                 S=S, R=R, mf=mf, m=m, a=a, da=da, N=np.exp(nu), nur=nur)
        transport = 16*ld('5.670400e-5')*T**3/(3*rho)
        z['Krad'], z['Kcond'] = transport/z['aux'][:, 24], transport/z['aux'][:, 27]
        z['K'] = z['Krad']+z['Kcond']
        z['tau_rad'] = 1/(e.C*rho*z['aux'][:, 24])
        z['beta_delta'] = np.column_stack([-np.log(z['Kcond']/zero['Kcond'])-2*delta[:, 1], -5*delta[:, 1]])
        self.finish_moments(delta, z, dE)
        return z

    def finish_moments(self, delta, z, dE):
        zero = self.reference
        rho_ratio, W, v = np.exp(delta[:, 0]), z['W'], z['v']
        dB = self.B0*(np.expm1(delta[:, 0])+rho_ratio*v*v*W*W/(W+1)
            +rho_ratio*W*z['da']/zero['a'])
        dS = z['w']*v*W*W+delta[:, 3]*self.qscale*(1+v*v)*W*W+2*zero['Q']*v*v*W*W
        dAS = z['a']*dS+z['da']*zero['Q']
        z['B'], z['AS'] = self.B0+dB, self.AS0+dAS
        z['dU'] = np.column_stack([dB, dE, dAS])
        z['dBX'] = dB[:, None]*self.base[:, 5:]+z['B'][:, None]*delta[:, 5:]
        z['fluxes'] = fluxes(self, z)
        # Differentiate the SAME finite-volume mass constraint used by E.
        mass_face_rate = -e.GRAV*e.C*4*np.pi*self.rf**2*z['fluxes'][1]
        mdot = mass_face_rate[:-1]+self.fraction*np.diff(mass_face_rate)
        z['at'] = z['a']**3*mdot/self.r
        self.current = z


def initialize(pool):
    star = parent.initialize(pool)
    initial = star.state(star.base)
    star.reference = {k:v.copy() if isinstance(v, np.ndarray) else v for k,v in initial.items()}
    star.B0 = initial['a']*initial['D']
    star.AS0 = initial['a']*initial['S']
    star.heat0 = initial['rho']*initial['aux'][:, 10]
    star.momentum_scale = initial['a']*initial['w']
    star.partial_volume = 4*np.pi/3*(star.r**3-star.rf[:-1]**3)
    star.fraction = star.partial_volume/star.volume
    star.beta0 = np.column_stack([np.log(e.TAU/initial['Kcond'])-2*star.base[:, 1],
        np.log(initial['tau_rad']/initial['Krad'])-2*star.base[:, 1]])
    star.__class__ = ConservativeStar
    return star


def residual(star, delta, previous, older, h, coefficients):
    z = star.evaluate(delta)
    c0,c1,c2 = coefficients
    dt = (c0*delta+c1*previous[0]+c2*older[0])/h
    dU = c0*z['dU']+c1*previous[1]['dU']+c2*older[1]['dU']
    fB,fE,fS = z['fluxes']
    result = np.zeros_like(delta)
    result[:, 0] = (dU[:, 0]+h*e.C*star.divergence(fB))/star.B0
    result[:, 1] = (dU[:, 1]+h*e.C*star.divergence(fE))/star.heat0
    force = -z['at']*z['S']-z['N']*z['nur']*e.C*z['E']+e.C*z['N']*z['P']*star.area_difference_over_volume
    result[:, 2] = (dU[:, 2]+h*e.C*star.divergence(fS)-h*force)/star.momentum_scale
    W,N,a,v,T = [z[k] for k in ['W','N','a','v','T']]
    speed = e.C*N*v/a
    vr, tr = star.gradient(v, odd=True), star.gradient(star.base[:, 1]+delta[:, 1])
    theta = W/N*(W*W*v*dt[:, 2]+e.C*N/a*W*W*vr+z['at']/a
        +e.C*N*v*z['nur']/a+2*e.C*N*v/(a*star.r))
    drive = (W*v/N*(dt[:, 1]+speed*tr)+e.C*tr/(a*W)
        +W**3/N*(dt[:, 2]+speed*vr)+e.C*W*z['nur']/a+W*v*z['at']/(N*a))
    beta_dt = (c0*z['beta_delta']+c1*previous[1]['beta_delta']+c2*older[1]['beta_delta'])/h
    beta_r = star.gradient(star.beta0+z['beta_delta'])
    qr, qrr = star.gradient(z['Q'], odd=True), star.gradient(z['Qrad'], odd=True)
    for column,q,qdot,qgrad,K,tau,index in [
        (3,z['Qcond'],star.qscale*(dt[:, 3]-dt[:, 4]),qr-qrr,z['Kcond'],e.TAU,0),
        (4,z['Qrad'],star.qscale*dt[:, 4],qrr,z['Krad'],z['tau_rad'],1)]:
        law = (tau*W/N*(qdot+speed*qgrad)+q+K*T/e.C**2*drive
            +tau*q/2*(theta+W/N*(beta_dt[:, index]+speed*beta_r[:, index])))
        result[:, column] = h*N*law/(tau*W*star.qscale)
    species_flux = fB[:, None]*upwind(star.base[:, 5:]+delta[:, 5:], fB)
    inventory = c0*z['dBX']+c1*previous[1]['dBX']+c2*older[1]['dBX']
    result[:, 5:] = (inventory+h*e.C*np.diff(4*np.pi*star.rf[:, None]**2*species_flux, axis=0)/star.volume[:, None])/star.B0[:, None]
    assert np.all(np.isfinite(result))
    return result, z


def species(star, candidate, previous, older, h, coefficients, metric):
    _,c1,c2 = coefficients
    v = candidate[:, 2]
    B = metric['a']*star.reference['rho']*np.exp(candidate[:, 0])/np.sqrt(1-v*v)
    velocity = star.faces(metric['N']*v/metric['a'], odd=True)
    flux = velocity*upwind(B, velocity)
    effective = -c1*previous[1]['B']-c2*older[1]['B']
    assert np.all(effective > 0)
    left = h*e.C*4*np.pi*star.rf[:-1]**2*np.maximum(flux[:-1], 0)/star.volume
    right = h*e.C*4*np.pi*star.rf[1:]**2*np.maximum(-flux[1:], 0)/star.volume
    band = np.zeros((3, star.n))
    band[1] = 1+(left+right)/effective
    band[0, 1:], band[2, :-1] = -right[:-1]/effective[:-1], -left[1:]/effective[1:]
    changes = np.diff(star.base[:, 5:], axis=0)
    source = left[:, None]*np.vstack([np.zeros(26), -changes])+right[:, None]*np.vstack([changes, np.zeros(26)])
    source += -c1*previous[1]['B'][:, None]*previous[0][:, 5:]-c2*older[1]['B'][:, None]*older[0][:, 5:]
    return solve_banded((1,1), band, np.asarray(source/effective[:, None], float)).astype(ld)


class LocalTangent(ConservativeStar):
    def evaluate(self, delta):
        change, aux = delta-self.anchor, self.linearization['aux'].copy()
        aux[:, 1] *= np.exp(aux[:, 5]*change[:, 0]+aux[:, 6]*change[:, 1])
        aux[:, 2] += aux[:, 9]*change[:, 0]+aux[:, 10]*change[:, 1]
        for j in [24,27]:
            aux[:, j] *= np.exp(aux[:, j+1]*change[:, 0]+aux[:, j+2]*change[:, 1])
        aux[:, 21] = aux[:, 24]*aux[:, 27]/(aux[:, 24]+aux[:, 27])
        y = self.base+delta
        self.material_cache = {e.material_key(row):value for row,value in zip(zip(y[:, 0],y[:, 1],y[:, 5:]),aux)}
        z = super().evaluate(delta)
        for key in ['m','mf','a','da','N','nur']:
            z[key] = self.linearization[key]
        self.finish_moments(delta, z, z['dU'][:, 1])
        return z


def tangent(star, delta, z):
    model = LocalTangent.__new__(LocalTangent)
    model.__dict__ = dict(star.__dict__)
    model.pool = parent.parent.old.CachedOnly()
    model.anchor, model.linearization = delta.copy(), z
    return model


def jacobian(model, delta, previous, older, h, coefficients):
    rows, columns, values = [], [], []
    step = ld('.001')
    for field in range(5):
        for color in range(5):
            selected = np.arange(color, model.n, 5)
            change = np.zeros_like(delta)
            change[selected, field] = step*SCALE[field]
            plus = residual(model, delta+change, previous, older, h, coefficients)[0][:, :5]
            minus = residual(model, delta-change, previous, older, h, coefficients)[0][:, :5]
            response = (plus-minus)/(2*step*SCALE)
            for offset in range(-2,3):
                target = selected+offset
                valid = (target >= 0)&(target < model.n)
                target, source = target[valid], selected[valid]
                rows.extend((5*target[:, None]+np.arange(5)).ravel())
                columns.extend(np.repeat(5*source+field,5))
                values.extend(np.asarray(response[target],float).ravel())
    return coo_matrix((values,(rows,columns)),shape=(5*model.n,5*model.n)).tocsc()


def stage(star, previous, older, h, coefficients, log):
    delta, history, factor = previous[0].copy(), [], None
    value,z = residual(star,delta,previous,older,h,coefficients)
    for iteration in range(24):
        star.last_delta, star.last_state = delta.copy(), z
        norm = float(np.max(abs(value)/ATOL))
        log(dict(iteration=iteration,residual_norm=norm,maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist()))
        if norm <= 1:
            return delta,z
        if iteration == 23:
            break
        if factor is None or iteration%4 == 0:
            factor = splu(jacobian(tangent(star,delta,z),delta,previous,older,h,coefficients))
        x,f = (delta[:, :5]/SCALE).ravel(),(value[:, :5]/SCALE).ravel()
        mixed,raw = parent.parent.direction(x,f,factor,history)
        accepted = None
        trials = []
        for proposal in ([mixed,raw] if history else [raw]):
            correction = proposal.reshape(star.n,5).astype(ld)*SCALE
            fraction = min(1.,.1/max(float(abs(correction[:, :2]).max()),1e-300),.01/max(float(abs(correction[:, 2]).max()),1e-300))
            for backtrack in range(8):
                candidate = delta.copy()
                candidate[:, :5] += fraction*correction
                try:
                    candidate[:, 5:] = species(star,candidate,previous,older,h,coefficients,z)
                    trial,state = residual(star,candidate,previous,older,h,coefficients)
                    score = float(np.max(abs(trial)/ATOL))
                    trials.append(dict(fraction=fraction,norm=score))
                    if score <= 1 or score < norm*(1-1e-4*fraction):
                        accepted = candidate,trial,state
                        break
                except (AssertionError,np.linalg.LinAlgError,FloatingPointError) as error:
                    trials.append(dict(fraction=fraction,invalid_trial=repr(error)))
                fraction /= 2
            if accepted is not None:
                break
        if accepted is None:
            raise RuntimeError(('Conservative native line search failed',norm,trials))
        history = (history+[(x.copy(),f.copy())])[-6:]
        delta,value,z = accepted
    raise RuntimeError(('Conservative native BDF stage did not converge',norm))


def symbolic():
    v,eps,P,Q = sp.symbols('v eps P Q')
    E,S,R = [(eps+P*v*v+2*Q*v)/(1-v*v),((eps+P)*v+Q*(1+v*v))/(1-v*v),(eps*v*v+P+2*Q*v)/(1-v*v)]
    assert sp.cancel(S-v*E-P*v-Q) == 0
    assert sp.cancel(R-v*S-P-Q*v) == 0
    r,h,t = sp.symbols('r h t',positive=True)
    c0,c1,c2 = (1+2*r)/(1+r),-(1+r),r*r/(1+r)
    for power in range(3):
        assert sp.simplify(c0*t**power+c1*(t-h)**power+c2*(t-h-h/r)**power-h*sp.diff(t**power,t)) == 0
    assert sp.simplify(-c2/c0+r*r/(1+2*r)) == 0
    # Species residual minus X times the baryon residual preserves sum X=1
    # throughout the nonlinear solve and is conservative when baryons close.
    X,Xp,Xo,B,Bp,Bo,F,FX = sp.symbols('X Xp Xo B Bp Bo F FX')
    assert sp.expand((c0*B*X+c1*Bp*Xp+c2*Bo*Xo+h*FX)-X*(c0*B+c1*Bp+c2*Bo+h*F)
        -(c1*Bp*(Xp-X)+c2*Bo*(Xo-X)+h*(FX-X*F))) == 0
    lr,lt,lv,lq,q,tau,W,N,H,force,kr,kt = sp.symbols('lr lt lv lq q tau W N H force kr kt')
    full = tau*W/N*lq+q+H*(W*v/N*lt+W**3/N*lv+force)+tau*q*W/(2*N)*(-lr-(kr-1)*lr-kt*lt)
    reduced = W/N*(-tau*q*kr/2*lr+(H*v-tau*q*kt/2)*lt+H*W*W*lv+tau*lq)+q+H*force
    assert sp.simplify(full-reduced) == 0
    assert parent.parent.symbolic()['passed']
    return dict(classification='Proven',passed=True,scope='Continuum flux decomposition, BDF2 degree-two consistency and species/baryon residual identity. Not a continuum, nonlinear stability, EOS or finite entropy certificate.')


def check():
    symbolic()
    errors = []
    for count in [16,32,64]:
        h = ld('.0001')/count
        value, older = ld(1), ld(1)
        for step in range(count):
            c0,c1,c2 = weights(h,None if step == 0 else h)
            value,older = (-c1*value-c2*older)/(c0+h*10000),value
        errors.append(float(abs(value-np.exp(ld(-1)))))
    orders = np.log2(np.array(errors[:-1])/errors[1:])
    assert np.all(orders > 1.8), orders
    # A contact with different rest coefficients must carry the SAME isotope
    # rest energy. The old central energy flux fails this positive-speed case.
    composition = np.array([[1.,0.],[0.,1.],[0.,1.]])
    coefficient = np.array([8.,11.])
    velocity = np.array([0.,.01,.01,0.])
    B = np.ones(3)
    class Contact:
        def faces(self,value,odd=False):
            result = np.r_[value[0],.5*(value[:-1]+value[1:]),value[-1]]
            if odd:
                result[[0,-1]] = 0
            return result
    lorentz = 1/np.sqrt(1-.01**2)
    energy = lorentz*(composition@coefficient)
    contact = dict(N=B,a=B,v=B*.01,B=B,E=energy,AS=energy*.01,P=B*0,Q=B*0,R=energy*.01**2)
    baryon_flux,energy_flux,_ = fluxes(Contact(),contact)
    species_flux = baryon_flux[:,None]*upwind(composition,velocity)
    assert np.allclose(energy_flux,lorentz*(species_flux@coefficient),rtol=0,atol=3e-17)
    assert abs(Contact().faces(energy*.01,odd=True)[1]-energy_flux[1]) > .01
    star = initialize(None)
    delta = np.zeros_like(star.base)
    initial = star.evaluate(delta)
    assert np.all(initial['dU'] == 0) and np.all(initial['dBX'] == 0)
    assert parent.cones(initial)['sampled_cone_inside_light_cone']
    for key in ['m','mf','a','N','Q','aux']:
        assert np.array_equal(initial[key],star.reference[key]),key
    previous = (delta,initial)
    h = ld('.00010282500222324415')
    coefficients = weights(h,None)
    model = tangent(star,delta,initial)
    J = jacobian(model,delta,previous,previous,h,coefficients)
    rng = np.random.default_rng(33)
    direction = rng.uniform(-1,1,size=(star.n,5)).astype(ld)
    change = np.zeros_like(delta)
    change[:, :5] = ld('.0001')*SCALE*direction
    finite = (residual(model,delta+change,previous,previous,h,coefficients)[0][:,:5]
        -residual(model,delta-change,previous,previous,h,coefficients)[0][:,:5])/(2*ld('.0001')*SCALE)
    predicted = J.dot(np.asarray(direction,float).ravel()).reshape(star.n,5)
    errors = np.max(abs(finite-predicted),axis=0)/np.maximum(np.max(abs(finite),axis=0),1e-30)
    assert np.max(errors) < 1e-3, errors
    result = dict(classification='Counterexample candidate',passed=True,symbolic=symbolic(),
        contact_rest_flux_consistent=True,initial_moments_and_geometry_preserved=True,
        tangent_direction_errors=errors.astype(float).tolist(),linear_BDF2_orders=orders.tolist(),actual_native_time_step=False)
    print('CONSERVATIVE TWO CARRIER CHECK',json.dumps(result),flush=True)
    return result


def prepare():
    assert not OUT.exists(), 'Preserve every previously bound method and failure.'
    checked = check()
    plan = json.loads((parent.OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest, rel
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        candidate='Native conservative baryon/energy/momentum/isotopes and two heat entropy laws.',
        method='Variable-step BDF2, initial backward Euler; original coordinate time nodes. Common upwind advection for aD,E,aS,X with centered pressure/heat fluxes. The native EOS appears in conservative endpoint increments. Differentiate the same finite-volume mass flux for a_t. Heat equations use the full native beta=tau/(K*T^2) change and explicit expansion, including composition dependence.',
        nonlinear_absolute_tolerances=ATOL.astype(float).tolist(),nonlinear_relative_tolerance=0,
        maximum_stage_iterations=24,implementation_check=checked,
        finite_conservation_gates=dict(local_energy_over_initial_heat_capacity=1e-6,
            local_baryon_relative=1e-9,isotope_inventory_over_total_initial_baryons=1e-9),
        time_refinement_gate='All five endpoint differences decrease; lnT,Qtotal,Qrad observed order at least 1.5. Require every full path, native residual and finite conservation gate. No prefix replacement.',
        correction='Use conservative a*S momentum with its -a_t*S source; equivalently the primitive source contains -2*a_t*S.',
        parent_plan_sha256=e.digest(parent.OUT/'plan.json'),
        limitations='Same conditional LTE two-carrier physics, prescribed conductive time, opacity conventions and reflecting wall. First-order upwind material advection and centered pressure/heat terms; not a well-balanced TOV reconstruction, Riemann solver, physical exterior, nonlinear entropy or continuum error certificate. Exact-arithmetic telescoping does not certify native EOS or the physical solution. No nuclear reactions or observational closure.')
    for path in [Path(__file__),parent.OUT/'plan.json',e.ROOT/'verification/gr_material_conservation.py']:
        plan['bindings'][path.relative_to(e.ROOT).as_posix()] = e.digest(path)
    OUT.mkdir()
    e.write(OUT/'plan.json',plan)
    print('PREPARED CONSERVATIVE NATIVE BDF2 EVOLUTION',flush=True)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():
        assert e.digest(e.ROOT/rel) == digest,rel
    for rel,digest in plan['runtime_sha256'].items():
        assert e.digest(Path(rel)) == digest,rel
    assert e.TAU == ld(plan['conduction_time'])
    return plan


def budget(star,z,energy_flux,baryon_flux):
    energy = (z['dU'][:,1]*star.volume+np.diff(energy_flux))/(star.heat0*star.volume)
    baryon = (z['dU'][:,0]*star.volume+np.diff(baryon_flux))/(star.B0*star.volume)
    inventory = np.sum(z['dBX']*star.volume[:,None],axis=0)/np.sum(star.B0*star.volume)
    return dict(maximum_local_energy_defect=float(abs(energy).max()),
        maximum_local_baryon_defect=float(abs(baryon).max()),
        maximum_isotope_inventory_defect=float(abs(inventory).max()),
        isotope_inventory_defects=inventory.astype(float).tolist(),
        total_energy_defect_over_initial_heat_capacity=float(np.sum(z['dU'][:,1]*star.volume)/np.sum(star.heat0*star.volume)),
        relative_total_baryon_defect=float(np.sum(z['dU'][:,0]*star.volume)/np.sum(star.B0*star.volume))),energy,baryon


def budget_passed(value,plan):
    gate = plan['finite_conservation_gates']
    return (value['maximum_local_energy_defect'] <= gate['local_energy_over_initial_heat_capacity']
        and value['maximum_local_baryon_defect'] <= gate['local_baryon_relative']
        and value['maximum_isotope_inventory_defect'] <= gate['isotope_inventory_over_total_initial_baryons'])


def save_state(path,delta,z,t,energy_flux,baryon_flux,energy_defect,baryon_defect):
    np.savez_compressed(path,delta=delta,**{k:z[k] for k in ['m','mf','a','N','Q','aux','dU']},
        time_seconds=t,integrated_energy_flux=energy_flux,integrated_baryon_flux=baryon_flux,
        normalized_energy_defect=energy_defect,normalized_baryon_defect=baryon_defect)


def run(refinement,workers):
    plan = bindings()
    assert refinement in plan['refinements'] and workers > 0
    folder = OUT/f'path-{refinement}'
    assert not folder.exists(), 'An existing path is never restarted or overwritten.'
    folder.mkdir()
    times = parent.wall.prior.time_nodes(plan,refinement)
    began = time.monotonic()
    with ProcessPoolExecutor(max_workers=workers,initializer=e.worker_init) as pool:
        star = initialize(pool)
        delta = np.zeros_like(star.base)
        z = star.evaluate(delta)
        previous = older = (delta,z)
        energy_flux = np.zeros(star.n+1,dtype=ld)
        baryon_flux = energy_flux.copy()
        older_energy,older_baryon = energy_flux.copy(),baryon_flux.copy()
        for step,t in enumerate(times):
            if step:
                h = t-times[step-1]
                coefficients = weights(h,None if step == 1 else times[step-1]-times[step-2])
                def log(record):
                    with (folder/'iterations.jsonl').open('a') as stream:
                        stream.write(json.dumps(dict(step=step,**record))+'\n')
                    print('CONSERVATIVE NATIVE STAGE',refinement,step,record['iteration'],record['residual_norm'],flush=True)
                try:
                    delta,z = stage(star,previous,older,h,coefficients,log)
                except Exception as error:
                    if hasattr(star,'last_delta'):
                        np.savez_compressed(folder/'failed-iterate.npz',delta=star.last_delta,aux=star.last_state['aux'],time_seconds=t)
                    e.write(folder/'failure.json',dict(classification='Counterexample candidate',completed=False,step=step,reason=repr(error)))
                    raise
                c0,c1,c2 = coefficients
                area = 4*np.pi*star.rf**2
                energy_flux,older_energy = (-c1*energy_flux-c2*older_energy+h*e.C*area*z['fluxes'][1])/c0,energy_flux
                baryon_flux,older_baryon = (-c1*baryon_flux-c2*older_baryon+h*e.C*area*z['fluxes'][0])/c0,baryon_flux
            value,energy_defect,baryon_defect = budget(star,z,energy_flux,baryon_flux)
            cone = parent.cones(z)
            with (folder/'cones.jsonl').open('a') as stream:
                stream.write(json.dumps(dict(step=step,**cone))+'\n')
            if not cone['sampled_cone_inside_light_cone'] or not budget_passed(value,plan):
                save_state(folder/'inadmissible-state.npz',delta,z,t,energy_flux,baryon_flux,energy_defect,baryon_defect)
                e.write(folder/'failure.json',dict(classification='Counterexample candidate',completed=False,step=step,
                    reason='Converged state failed a frozen cone or conservation gate.',budget=value,cone=cone))
                raise RuntimeError('Conservative state outside the frozen admissible domain')
            save_state(folder/f'step-{step:04d}.npz',delta,z,t,energy_flux,baryon_flux,energy_defect,baryon_defect)
            row = dict(classification='Counterexample candidate',step=step,time_seconds=float(t),
                elapsed_seconds=time.monotonic()-began,maximum_changes=np.max(abs(delta[:,:5]),axis=0).astype(float).tolist(),**value,**cone)
            e.write(folder/'progress.json',row)
            print('CONSERVATIVE TWO CARRIER STEP',refinement,json.dumps(row),flush=True)
            if step:
                older,previous = previous,(delta,z)
        e.write(folder/'result.json',dict(classification='Counterexample candidate',completed=True,steps=len(times)-1,
            duration_seconds=float(times[-1]),physical_EOS_certified=False,physical_exterior_match=False,observational_closure=False))
    e.write(folder/'manifest.json',dict(plan_sha256=e.digest(OUT/'plan.json'),sha256={p.name:e.digest(p) for p in folder.iterdir() if p.is_file()}))


def completed(refinement):
    plan = bindings()
    folder = OUT/f'path-{refinement}'
    result = json.loads((folder/'result.json').read_text())
    manifest = json.loads((folder/'manifest.json').read_text())
    assert result['completed'] and not (folder/'failure.json').exists()
    assert manifest['plan_sha256'] == e.digest(OUT/'plan.json')
    for rel,digest in manifest['sha256'].items():
        assert e.digest(folder/rel) == digest,rel
    times = parent.wall.prior.time_nodes(plan,refinement)
    assert result['steps'] == len(times)-1 and result['duration_seconds'] == float(times[-1])
    iterations = [json.loads(line) for line in (folder/'iterations.jsonl').read_text().splitlines()]
    assert {r['step'] for r in iterations} == set(range(1,len(times)))
    for step in range(1,len(times)):
        rows = [r for r in iterations if r['step'] == step]
        assert [r['iteration'] for r in rows] == list(range(len(rows))) and len(rows) <= 24
        assert rows[-1]['residual_norm'] <= 1 and len(rows[-1]['maximum_absolute_residual']) == 31
    star = initialize(None)
    previous = older = None
    energy_flux = baryon_flux = np.zeros(star.n+1,dtype=ld)
    older_energy,older_baryon = energy_flux.copy(),baryon_flux.copy()
    maximum_energy = maximum_baryon = maximum_isotope = maximum_speed = 0.
    for step,t in enumerate(times):
        cp = np.load(folder/f'step-{step:04d}.npz')
        assert cp['time_seconds'] == t
        delta = cp['delta'].copy()
        y = star.base+delta
        star.material_cache = {e.material_key(row):aux for row,aux in zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])}
        z = star.evaluate(delta)
        for key in ['m','mf','a','N','Q','aux','dU']:
            assert np.array_equal(z[key],cp[key]),(step,key)
        value,energy_defect,baryon_defect = budget(star,z,cp['integrated_energy_flux'],cp['integrated_baryon_flux'])
        assert np.array_equal(energy_defect,cp['normalized_energy_defect'])
        assert np.array_equal(baryon_defect,cp['normalized_baryon_defect'])
        assert budget_passed(value,plan)
        assert np.all(cp['integrated_energy_flux'][[0,-1]] == 0) and np.all(cp['integrated_baryon_flux'][[0,-1]] == 0)
        cone = parent.cones(z)
        assert cone['sampled_cone_inside_light_cone']
        maximum_speed = max(maximum_speed,cone['maximum_local_rest_characteristic_speed_over_c'])
        maximum_energy = max(maximum_energy,value['maximum_local_energy_defect'])
        maximum_baryon = max(maximum_baryon,value['maximum_local_baryon_defect'])
        maximum_isotope = max(maximum_isotope,value['maximum_isotope_inventory_defect'])
        if step:
            h = t-times[step-1]
            coefficients = weights(h,None if step == 1 else times[step-1]-times[step-2])
            native,_ = residual(star,delta,previous,older,h,coefficients)
            assert np.max(abs(native)/ATOL) <= 1,(step,float(np.max(abs(native)/ATOL)))
            c0,c1,c2 = coefficients
            area = 4*np.pi*star.rf**2
            energy_flux,older_energy = (-c1*energy_flux-c2*older_energy+h*e.C*area*z['fluxes'][1])/c0,energy_flux
            baryon_flux,older_baryon = (-c1*baryon_flux-c2*older_baryon+h*e.C*area*z['fluxes'][0])/c0,baryon_flux
            older,previous = previous,(delta,z)
        else:
            assert np.all(delta == 0)
            previous = older = (delta,z)
        assert np.array_equal(energy_flux,cp['integrated_energy_flux'])
        assert np.array_equal(baryon_flux,cp['integrated_baryon_flux'])
    return delta,dict(refinement=refinement,steps=len(times)-1,duration_seconds=float(times[-1]),
        maximum_path_speed_over_c=maximum_speed,maximum_path_local_energy_defect=maximum_energy,
        maximum_path_local_baryon_defect=maximum_baryon,maximum_path_isotope_inventory_defect=maximum_isotope)


def compare():
    assert not (OUT/'time-refinement.json').exists(), 'Preserve a failed full-time verdict.'
    paths = [completed(r) for r in [1,2,4]]
    errors = np.array([np.max(abs(b[0][:,:5]-a[0][:,:5]),axis=0) for a,b in zip(paths[:-1],paths[1:])])
    passed,orders = parent.audit.gate(errors)
    result = dict(classification='Counterexample candidate',passed=passed,
        endpoint_maximum_differences=errors.astype(float).tolist(),observed_orders=orders.tolist(),paths=[p[1] for p in paths],
        rigorous_time_error_bound=False,physical_EOS_certified=False,physical_exterior_match=False,observational_closure=False)
    e.write(OUT/'time-refinement.json',result)
    files = [OUT/'plan.json',OUT/'time-refinement.json']+[OUT/f'path-{r}/manifest.json' for r in [1,2,4]]
    e.write(OUT/'time-refinement-manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files}))
    print('CONSERVATIVE FULL TIME COMPARISON',json.dumps(result),flush=True)
    assert passed,'The frozen full-time convergence gate failed.'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command',choices=['check','prepare','run','completed','compare'])
    parser.add_argument('--refinement',type=int,choices=[1,2,4],default=1)
    parser.add_argument('--workers',type=int,default=4)
    args = parser.parse_args()
    if args.command == 'run':
        run(args.refinement,args.workers)
    elif args.command == 'completed':
        print(json.dumps(completed(args.refinement)[1],indent=2))
    else:
        globals()[args.command]()

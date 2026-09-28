"""Independent initial-data controls and ownership of unresolved source errors."""
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_initial_constraints as task

OUT=task.OUT;write=task.write;LD=task.LD


def run():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,task.previous.old.optical.timeout);signal.alarm(35)
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',seconds=35,native_calls=100,
        controls=['analytic uniform-density polar slice','first-law and Einstein identities','exact-source Gauss panel flux residual','native isentropic anchor recovery','represented-cell baryon inventory versus actual fluid mass'],
        source_sha256=task.sha(__file__),producer_sha256=task.sha(task.__file__),
        accept='Initial mathematical constraint solution and actual finite-volume source matching are separate gates. A per-cell inventory mismatch must be fixed before reusing this state as a new physical flow.'))
    r,e,p=sp.symbols('r e p',positive=True);k=8*sp.pi*e/3;mass=4*sp.pi*e*r**3/3;b=1-k*r*r
    nu=-sp.log(b)/4;assert sp.simplify(sp.diff(nu,r)-mass/(r*r*b))==0
    m,Phi,A,En,Pr,F,N,c,G=sp.symbols('m Phi A En Pr F N c G',positive=True)
    bb=1-2*m/r;mt=-4*sp.pi*r*r*N*sp.sqrt(bb)*A**4*G*F/c**4
    Krr=4*sp.pi*r*A**4*G*F/(c**5*sp.sqrt(bb))
    assert sp.simplify(mt+c*N*r*bb*Krr)==0
    rho,w,gamma,E0,p0=sp.symbols('rho w gamma E0 p0',positive=True)
    Eg=E0*w+p0*(w**gamma-w)/(gamma-1);Pg=p0*w**gamma
    assert sp.simplify(w*sp.diff(Eg,w)-Eg-Pg)==0
    q=task.Quadrature(np.linspace(0,1,23),12);rr=q.r;en=LD('.001');actual,faces=q.integrate(4*LD(np.pi)*en*rr**2)
    exact=4*LD(np.pi)*en*rr**3/3;uniform=float(max(abs(actual-exact).flat)/max(abs(exact).flat))
    assert uniform<1e-11
    data=task.Data();fine=np.load(OUT/'balanced-20.npz');z=np.load(OUT/'balanced-initial-state.npz');q=task.Quadrature(fine['edges'],20)
    # The primitive flux balance uses all panel contributions. No derivative
    # of a nearly constant rounded phi profile enters this independent check.
    increments=q.h.astype(LD)*(fine['panel_scalar_source']@q.w.astype(LD))
    scalar=float(max(abs(np.diff(fine['scalar_flux_faces'])-increments))/max(abs(increments)))
    assert scalar<1e-11
    # Isentropic native calls use the same frozen ionic inventory and entropy
    # as the actual initial local EOS, not a fitted Gamma pressure surrogate.
    bulk=data.model.bulk;d=bulk.d;native=task.previous.chem.old.Native(cap=100);controls=[]
    for j in range(bulk.n):
        task.previous.chem.setup(native,d,j);target=d['raw'][j,3];lt=np.log(d['T'][j]);x=float(z['density_log_ratio'][j]);y=d['y0'][j]
        for _ in range(4):
            state=native.state(x,lt,y);raw=state['raw'];delta=(raw[3]-target)*np.exp(lt)/raw[10]
            if abs(delta)<2e-12:break
            lt-=delta
        else:raise AssertionError(('Native entropy recovery',j,delta))
        proper_u=float(z['gas_energy'][j]/z['density'][j]-data.model.cx*task.C**2)
        controls.append(dict(cell=j,logrho=x,temperature_K=float(np.exp(lt)),entropy_over_cv=float(abs(delta)),
            pressure=float(abs(z['gas_pressure'][j]/raw[1]-1)),internal_energy=float(abs(proper_u/raw[2]-1))))
    native_error=max(max(c['pressure'],c['internal_energy']) for c in controls)
    assert native_error<.002,native_error
    # Preserve the material inventory of the ACTUAL solver, not merely of an
    # arbitrary continuous interpolation through its centers.
    source=data.source(q.r.ravel());phi0=data.scalar(q.r.ravel())[0].reshape(q.r.shape).astype(LD)
    old=data.bg.sample(q.r.ravel()/data.R);bb=1-2*old['m'].reshape(q.r.shape).astype(LD)*LD(data.R)/q.r
    baryon=4*LD(np.pi)*q.r**2*np.exp(-6*phi0**2)*source['rho'].reshape(q.r.shape).astype(LD)/np.sqrt(bb)
    _,enclosed=q.integrate(baryon)
    jordan_edges=np.r_[bulk.d['edges'][:-1],data.model.m.rf]
    edges=data.geometry.metric(jordan_edges-data.model.m.RJ)[3]*data.R
    # These faces need not be old source knots. Integrate the stored positive
    # source directly on each physical face interval with20 Gauss nodes.
    qb=task.Quadrature(edges,20);rb=qb.r;ph=data.scalar(rb.ravel())[0].reshape(rb.shape).astype(LD)
    old=data.bg.sample(rb.ravel()/data.R);bb=1-2*old['m'].reshape(rb.shape).astype(LD)*LD(data.R)/rb
    rho=data.source(rb.ravel())['rho'].reshape(rb.shape).astype(LD)
    reconstructed=qb.h*(4*LD(np.pi)*rb**2*np.exp(-6*ph**2)*rho/np.sqrt(bb)@qb.w)
    actual=np.r_[data.model.mass0,data.model.flow.initial[0]*data.model.flow.eos.rho0*4*np.pi*data.model.m.RJ**2*data.model.m.vol].astype(LD)
    active=actual>0;relative=np.zeros_like(actual);relative[active]=reconstructed[active]/actual[active]-1
    scale=max(abs(actual));absolute=float(max(abs(reconstructed-actual))/scale)
    np.savez_compressed(OUT/'inventory-audit.npz',actual_g=actual,reconstructed_g=reconstructed,relative=relative,edges_E=edges)
    row=dict(classification='Counterexample candidate',mathematical_projection_passed=bool(uniform<1e-11 and scalar<1e-11 and native_error<.002),
        uniform_star_relative=uniform,scalar_flux_panel_relative=scalar,native_constitutive_relative=native_error,native_calls=native.ion.calls,native_anchors=controls,
        actual_cell_inventory_max_relative=float(max(abs(relative))),actual_inventory_absolute_over_largest_cell=absolute,
        actual_deep_total_g=float(actual[:bulk.n].sum()),reconstructed_deep_total_g=float(reconstructed[:bulk.n].sum()),
        exact_actual_cell_inventory_matched=bool(max(abs(relative))<1e-10),production_ready=bool(max(abs(relative))<1e-10),
        momentary_balance_is_not_stationary_radiating_star=True,full_native_continuum_EOS=False,final_charge_solved=False,seconds=time.monotonic()-start)
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Polar momentum identity, uniform pressureless Cauchy mass/lapse identity and the anchored Gamma first law. These are algebraic statements, not empirical EOS or flow certification.'))
    write(OUT/'audit.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

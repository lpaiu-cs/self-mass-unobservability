"""Material metric work and stable rest-energy subtraction for evolving GR cells."""
import json,sys
import numpy as np
import sympy as s
import mpmath as mp
import gr_material_conservation as material

g=material.g;OUT=g.OUT/'gr-material-rest-energy'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();material.verify()
    paths=[g.ROOT/'verification/gr_material_rest_energy.py',material.OUT/'manifest.json',
        g.OUT/'gr-heat-initial-constraints/initial-data.npz',g.OUT/'initial-state-17-4.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='b84cb63',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        tangent_probe_seconds=[1e-6,1e-12],arithmetic_relative_tolerance=1e-12,
        assumptions='Smooth polar-areal spherical GR, c=G=1 for algebra; material baryon measure dB conserved, composition/rest constant C within each cell. Positive f=1-2m/r at both compared states.',
        scope='Exact conditional change-of-energy-variable identity plus arithmetic tests of saved initial metric tangents. The probes are not solutions of the coupled GR equations and do not certify continuum truncation or physical EOS.'))
    r,m,dr,dm,N,v,P,Q,B,C=s.symbols('r m dr dm N v P Q B C',positive=True)
    f=1-2*m/r;df=2*(m*dr-r*dm)/(r*(r+dr))
    assert s.factor((1-2*(m+dm)/(r+dr))-f-df)==0
    a=1/s.sqrt(f);rdot=N*v/a;mdot=-4*s.pi*r*r*N/a*(P*v+Q)
    inverse_a_rate=s.diff(s.sqrt(f),r)*rdot+s.diff(s.sqrt(f),m)*mdot
    expected=N*(v*(m/r**2+4*s.pi*r*P)+4*s.pi*r*Q)
    assert s.simplify(inverse_a_rate-expected)==0
    # Work from the independently derived material mass-flux law.
    assert s.simplify(inverse_a_rate.subs(v,0)-4*s.pi*r*N*Q)==0
    F=s.symbols('F0:5');k=s.symbols('k0:4');b=s.symbols('b0:4');c=s.symbols('c0:4')
    ubar=[F[i+1]-F[i]-c[i]*b[i]*k[i] for i in range(4)]
    assert s.expand(sum(ubar[i]+c[i]*b[i]*k[i] for i in range(4))-(F[4]-F[0]))==0
    x,y=s.symbols('x y',positive=True)
    assert s.simplify((s.sqrt(y)-s.sqrt(x))*(s.sqrt(y)+s.sqrt(x))-(y-x))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        rest_coordinate_weight='B*kappa=integral[D dVcoord]=integral[(1/a) dB], with dB=a*D*dVcoord. This identity applies to arbitrary smooth material profiles, including nonzero W; an extra 1/W would be incorrect.',
        metric_weight_rate='B*kappa_dot=integral N*[v*(m/r^2+4*pi*r*P)+4*pi*r*Q] dB. Both R_dot=N*v/a and material m_dot=-4*pi*r^2*N/a*(P*v+Q) are used.',
        stable_shell_energy='For fixed B and C in one cell, Ubar=U-C*B*kappa, dUbar/dt=F_M(in)-F_M(out)-C*B*kappa_dot. At finite changes, Delta Ubar=Delta U-C*integral Delta(1/a) dB. Preserve increments instead of subtracting two large total rest energies.',
        stable_metric_increment='delta_f=2*(m*delta_r-r*delta_m)/(r*(r+delta_r)); delta(1/a)=delta_f/(sqrt(f+delta_f)+sqrt(f)). This avoids both total-mass subtraction and subtraction of two near-unity square roots.',
        time_zero_velocity='At instantaneous v=0, kappa_dot still has integral 4*pi*r*N*Q dB/B. Setting kappa constant suppresses an actual metric-energy contribution when Q is nonzero.',
        species_limit='If C changes due to composition, include the full change of C/a under the material integral; fixed-C formula is not a nuclear burning energy law.',
        limitation='Shared fluxes telescope in U=Ubar+C*B*kappa. These identities alone do not certify a numerical metric constraint, primitive closure, atmosphere or time path.'))
    plan=json.loads((OUT/'plan.json').read_text());initial=dict(np.load(g.OUT/'gr-heat-initial-constraints/initial-data.npz'))
    radius=initial['r_cm'];mass=initial['m_geom_cm'];mdot=initial['dmass_geom_cm_dt']
    mp.mp.dps=90;rows=[]
    for dt in plan['tangent_probe_seconds']:
        delta=mdot*np.longdouble(dt);r=radius.astype(np.longdouble);m=mass.astype(np.longdouble)
        base=1-2*m/r;increment=-2*delta/r
        stable=increment/(np.sqrt(base+increment)+np.sqrt(base))
        naive=np.sqrt(1-2*(mass+np.asarray(delta,float))/radius)-np.sqrt(1-2*mass/radius)
        errors=[]
        for rr,mm,dd,actual in zip(radius,mass,delta,stable):
            rmp=mp.mpf(float(rr));mmp=mp.mpf(float(mm))
            # Convert the stored longdouble exactly through its integer ratio.
            numerator,denominator=dd.as_integer_ratio();dmp=mp.mpf(numerator)/denominator
            baseline=1-2*mmp/rmp;change=-2*dmp/rmp
            reference=change/(mp.sqrt(baseline+change)+mp.sqrt(baseline))
            numerator,denominator=actual.as_integer_ratio();got=mp.mpf(numerator)/denominator
            errors.append(float(abs(got/reference-1)) if reference else float(abs(got)))
        row=dict(seconds=dt,states=len(r),maximum_relative_high_precision_difference=max(errors),
            nonzero_stable_changes=int(np.count_nonzero(stable)),naive_zero_changes=int(np.count_nonzero((naive==0)&(stable!=0))))
        row['passed']=row['maximum_relative_high_precision_difference']<=plan['arithmetic_relative_tolerance'];rows.append(row)
        np.savez_compressed(OUT/f'tangent-{dt}.npz',mass_increment_geom_cm=delta,stable_inverse_metric_increment=stable,
            naive_inverse_metric_increment=naive,high_precision_relative_difference=np.array(errors))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_passed=all(r['passed'] for r in rows),full_GR_evolution=False,physical_EOS_certified=False,
        arbitrary_precision_comparison_is_interval_proof=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('MATERIAL REST ENERGY',rows,flush=True)


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS material rest-energy change and stored arithmetic bindings',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

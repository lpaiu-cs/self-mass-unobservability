"""Replay Phase94 verdict and check the moving-surface mass jump independently."""
import json
import time
import numpy as np
import sympy as s
import def_native_wave_collocation as task


def junction():
    r,A,e,p,b,Phi,psi,z,Ji,Jo=s.symbols('r A e p b Phi psi z Ji Jo')
    inside=r*r*b*Phi*psi-4*s.pi*r**3*A*(e+p)*z+Ji
    outside=r*r*b*Phi*psi+Jo
    # Equality of mass at the displaced surface, including the background
    # material/vacuum derivative jump. A4 is denoted A in this identity.
    defect=s.expand(outside-inside-4*s.pi*r**3*A*e*z)
    required=s.solve(defect,Jo)[0]
    assert s.simplify(required-Ji+4*s.pi*r**3*A*p*z)==0
    R,G,c,N,a,P,E=s.symbols('R G c N a P E',nonzero=True)
    W=4*s.pi*R**3*A*N*a*P*z
    conversion=G/(c**4*R*N*a)
    assert s.simplify(-conversion*W+4*s.pi*A*G*R**2*P*z/c**4)==0
    assert task.symbolic()['passed']
    return dict(classification='Proven',passed=True,
        mass_jump='J_out-J_in=-4*pi*r^3*A4*p*zeta for a moving material/vacuum reference interface.',
        energy='With J=-G*E/(c^4*R*N*a) at the surface, E_out=E_heat+4*pi*R^3*A4*N*a*P*zeta.',
        pressure_split='Only Prad work belongs to outgoing photons; Pgas work requires the specified external stress reservoir.',
        scope='First-order identities at the reference surface. This calculation does not solve its motion, radiation geometry or the external gas reservoir.')


def main():
    start=time.monotonic();out=task.OUT
    plan=json.loads((out/'plan.json').read_text());result=json.loads((out/'result.json').read_text())
    for path,sha in plan['bindings'].items():assert task.old.photons.digest(task.old.ROOT/path)==sha,path
    a=np.load(out/'p4-32.npz');b=np.load(out/'p4-64.npz');c=np.load(out/'p2-64.npz')
    bg=task.old.Background();p=bg.sample(b['radius']);rows={};accepted=True
    for field in ['temperature','velocity','scalar']:
        norm=np.max(abs(b[field]));time_error=float(np.max(abs(a[field]-b[field][::2]))/norm)
        difference=abs(c[field]-b[field]);idx=np.unravel_index(np.argmax(difference),difference.shape)
        space_error=float(difference[idx]/norm);gate=.02 if field=='temperature' else .03
        expected=result['comparisons'][field]
        assert abs(time_error-expected['time_relative'])<1e-14
        assert abs(space_error-expected['space_relative'])<1e-14
        assert expected['time_pass']==(time_error<gate) and expected['space_pass']==(space_error<gate)
        accepted &= time_error<gate and space_error<gate
        j,i=idx
        cs=task.old.C*np.sqrt(p['gamma'][i]*p['p'][i]/(p['e'][i]+p['p'][i]))
        distance=min(b['radius'][i]-b['edges'][i],b['edges'][i+1]-b['radius'][i])*bg.R
        rows[field]=dict(time_relative=time_error,space_relative=space_error,
            worst_time_s=float(b['emission_times'][j]*bg.tc),native_array_index=int(i),
            radius_fraction=float(b['radius'][i]),p2=float(c[field][idx]),p4=float(b[field][idx]),
            reference_maximum=float(norm),local_acoustic_crossing_estimate_s=float(distance/cs))
    assert bool(accepted)==result['passed'] and not accepted
    assert not (out/'consistent-p4-64.npz').exists()
    energy={}
    for label in ['pilot','p4-32','p4-64','p2-64']:
        x=np.load(out/(label+'.npz'));report=json.loads((out/(label+'.json')).read_text())
        assert all(np.all(np.isfinite(x[k])) for k in x.files)
        assert all(np.max(abs(x[k][0]))==0 for k in rows)
        assert report['max_linear_residual']<1e-9 and report['max_heat_identity']<1e-9
        E=x['E'];balance=float(abs(np.sum(-np.diff(np.r_[np.longdouble(0),E]))+E[-1])/max(abs(E)))
        assert balance<1e-12;energy[label]=balance
    env=np.load(task.old.prior.OUT/'final-envelope.npz');surface=np.flatnonzero(b['grid']==1)[0]
    z=float(b['q'][b['indices'][surface,0]]);Na=float(env['N'][-1])/np.sqrt(float(env['b'][-1]))
    pref=4*np.pi*bg.R**3*float(env['A'][-1])**4*Na*z
    # Independent local flux check, using the stored envelope flux and pressure.
    radiation_pressure=float(env['Prad'][-1]);hemisphere_pressure=2*float(env['F'][-1])/(3*task.old.C)
    match=abs(hemisphere_pressure/radiation_pressure-1);assert match<1e-7
    report=dict(classification='Counterexample candidate',artifact_checks_passed=True,
        whole_response_passed=False,comparisons=rows,heat_balance=energy,symbolic=junction(),
        source_front_hypothesis=dict(classification='Conjectural',
            statement='Unresolved waves from discontinuous cellwise heat sources remain a candidate cause. A local sound-crossing estimate is not a rigorous domain-of-dependence certificate in the coupled GR/heat system.',confirmed=False),
        surface=dict(displacement_fraction=z,gas_pressure_dyn_cm2=float(env['Pgas'][-1]),
            radiation_pressure_dyn_cm2=radiation_pressure,hemisphere_pressure_relative=match,
            omitted_radiation_work_erg=pref*radiation_pressure,
            omitted_external_gas_work_erg=pref*float(env['Pgas'][-1]),
            work_over_emitted_heat=pref*float(env['Ptotal'][-1])/float(b['E'][-1]),
            isolated_matter_photon_stress_match=False,work_applied_to_evolution=False),
        seconds=time.monotonic()-start,full_goal_complete=False)
    task.write(out/'audit.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':main()

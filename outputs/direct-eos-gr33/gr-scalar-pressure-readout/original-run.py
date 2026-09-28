def run():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_scalar_reference.py',g.ROOT/'verification/common_eos.py',
        g.OUT/'initial-state-17-4.npz',g.OUT/'gr-microphysics/auxiliaries.npz',
        g.OUT/'gr-increment-structure/manifest.json',g.c.OUT/'scalar-connection.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='5090e6f',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        beta=-4,tolerances=[1e-10,2e-12],relative_coefficient_tolerance=1e-8,
        equation='Use the existing linearized DEF scalar equation on the current initial GR background and its declared PCHIP coefficients. (x^2*a*phi_prime)_prime=x^2*b*phi, including the exact Schwarzschild exterior logarithmic tail. a=N*sqrt(f), b=4*pi*beta*R^2*N*(epsilon-3P)/sqrt(f).',
        interpretation='Response coefficient alpha_A/phi_infinity at infinitesimal nonzero boundary scalar. No finite scalar amplitude is fitted or assumed. Independently integrate the Riccati flux/value equation and compare two ODE tolerances.',
        limits='Current fixed reference only. No actual orbital drive, relaxation response, finite scalar backreaction, physical EOS or atmosphere, continuous error enclosure, observed force or likelihood claim. A static response coefficient is not a dynamic-chi observable.'))
    x=s.symbols('x',positive=True);phi=s.Function('phi')(x);a=s.Function('a')(x);b=s.Function('b')(x)
    flux=x*x*a*s.diff(phi,x);w=flux/phi
    residual=s.diff(w,x)-(x*x*b-w*w/(x*x*a))
    equation=s.solve(s.Eq(s.diff(flux,x),x*x*b*phi),s.diff(phi,x,2))[0]
    assert s.simplify(residual.subs(s.diff(phi,x,2),equation))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        identity='Where phi is nonzero, w=(x^2*a*phi_prime)/phi obeys w_prime=x^2*b-w^2/(x^2*a). At x=1, A=-w/[1-w*log(1-2mu)/(2mu)] is the exterior dimensionless 1/x tail. The linear response coefficient is -R*A/M.',
        limitation='Riccati can fail at field nodes. The current weak-field reference is checked for no sign change on the saved sample, not globally certified node-free.'))
    control=riccati(dict(p=lambda x:1.,v=lambda x:-.3**2),2e-12)
    expected=.3/np.tan(.3)-1
    control_error=float(abs(control.y[0,-1]/expected-1));assert control_error<1e-9
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));native=np.load(g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
    assert np.allclose(native[:,0],np.exp(state['lnd']),rtol=1e-12,atol=0)
    assert np.array_equal(native[:,2],state['u_W'])
    plan=json.loads((OUT/'plan.json').read_text());rows=[]
    for tolerance in plan['tolerances']:
        model=g.c.scalar_model('current-native-initial',state,tolerance)
        nonlinear=riccati(model,tolerance);ratio=nonlinear.y[0,-1];mu=model['mu']
        amplitude=-ratio/(1-ratio*np.log1p(-2*mu)/(2*mu))
        coefficient=-model['chi']/model['M'];comparison=-model['R']*amplitude/model['M']
        field,derivative=model['field'](model['x'][1:]);assert np.all(field>0)
        row=dict(tolerance=tolerance,tail_chi_m=float(model['chi']),mass_geom_m=model['M'],radius_m=model['R'],
            alpha_over_phi_infinity=float(coefficient),Riccati_alpha_over_phi_infinity=float(comparison),
            relative_two_equation_difference=float(abs(comparison/coefficient-1)),
            minimum_sampled_normalized_field=float(field.min()),maximum_sampled_normalized_field=float(field.max()))
        row['passed']=row['relative_two_equation_difference']<=plan['relative_coefficient_tolerance'];rows.append(row)
        np.savez_compressed(OUT/f'field-{tolerance}.npz',x=model['x'],p_coefficients=model['p'].c,
            b_coefficients=model['v'].c,field=field,field_derivative=derivative,
            Riccati_x=nonlinear.t,Riccati_flux_ratio=nonlinear.y[0])
    difference=abs(rows[-1]['alpha_over_phi_infinity']/rows[0]['alpha_over_phi_infinity']-1)
    save('result.json',dict(classification='Counterexample candidate',completed=True,beta=-4,rows=rows,
        finite_tolerance_relative_difference=difference,
        all_passed=all(r['passed'] for r in rows) and difference<=plan['relative_coefficient_tolerance'],
        analytic_uniform_potential_control_relative_error=control_error,
        nonzero_linear_response_coefficient=bool(rows[-1]['alpha_over_phi_infinity']!=0),
        specified_finite_background_amplitude=False,nonlinear_scalar_backreaction=False,
        dynamic_orbital_response=False,physical_EOS_certified=False,observation_inference_completed=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify();print('CURRENT GR SCALAR REFERENCE',rows,flush=True)

def run():
    mp.mp.dps=80;mp.iv.dps=80
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    control=controls();save('controls.json',control)
    data=dict(np.load(regular.original.OUT/'field-2e-12.npz'));rat=regular.rational
    x=list(map(rat,data['x']));field=list(map(rat,data['field']));derivative=list(map(rat,data['field_derivative']))
    assert len(field)==len(x)-1==len(derivative)
    a0=rat(data['p_coefficients'][-1,0]);b0=rat(data['b_coefficients'][-1,0])
    scale=field[0]/(1+b0*x[1]**2/(6*a0));assert scale>0
    values=[Q(1)]+[v/scale for v in field];slopes=[Q(0)]+[v/scale for v in derivative]
    flux=Q(0);error=Q(0);rows=[]
    for i,(left,right) in enumerate(zip(x,x[1:])):
        h=right-left;change=values[i+1]-values[i]
        p=[values[i],slopes[i],3*change/h**2-(2*slopes[i]+slopes[i+1])/h,
            -2*change/h**3+(slopes[i]+slopes[i+1])/h**2]
        assert at(p,h)==values[i+1] and at([Q(j)*v for j,v in enumerate(p) if j],h)==slopes[i+1]
        a=list(map(rat,data['p_coefficients'][::-1,i]));b=list(map(rat,data['b_coefficients'][::-1,i]))
        flux,err=piece(left,h,a,b,p,flux);error+=ceil_bound(err)
        rows.append(dict(piece=i,exact_residual_bound=str(err),approximate_residual_bound=float(err)))
    prior=json.loads((regular.OUT/'result.json').read_text());r={k:Q(v) for k,v in prior['exact_rationals'].items()}
    L=-mp.iv.ln(1-2*regular.iv(r['mu']))/(2*regular.iv(r['mu']))
    enclosure,delta,normal=response(values[-1],flux,error,r['kappa'],r['B2'],L,r['R_over_M'])
    old=json.loads((regular.original.OUT/'result.json').read_text())['rows'][-1]['alpha_over_phi_infinity']
    value=regular.iv(rat(old));inside=bool(enclosure.a<=value.a and value.b<=enclosure.b)
    save('piece-residuals.json',dict(classification='Proven',rows=rows))
    save('result.json',dict(classification='Proven',completed=True,pieces=len(rows),
        exact_rationals=dict(interpolant_scale=str(scale),p_surface=str(values[-1]),J_surface=str(flux),
            residual_bound=str(error),field_error_bound=str(delta)),
        approximate_residual_bound=float(error),approximate_field_error_bound=float(delta),
        normalization_interval=regular.interval_text(normal),response_interval=regular.interval_text(enclosure),
        approximate_response_interval=[float(enclosure.a),float(enclosure.b)],
        maximum_piece_residual=max(rows,key=lambda r:r['approximate_residual_bound'])['piece'],
        prior_numerical_coefficient=old,prior_numerical_coefficient_inside=inside,
        analytic_control_passed=True,physical_EOS_certified=False,GR_coefficient_error_certified=False,
        finite_amplitude_backreaction=False,dynamic_orbital_observable=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()

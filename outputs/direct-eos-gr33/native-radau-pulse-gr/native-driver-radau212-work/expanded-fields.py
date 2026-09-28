def fields():
    assert read(OUT/'sources.json')['representation_controls_passed']
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    m=Response();fn=FunctionType(base.gr.base.Response.run.__code__,dict(base.gr.base.Response.run.__globals__,OUT=OUT/'gr'))
    rows=[read(OUT/'gr'/f'fields-{n}-g{q}.json') for n,q in [(128,8),(64,8),(128,4)]]
    fine,coarse,low=[dict(np.load(OUT/'gr'/f'fields-{n}-g{q}.npz')) for n,q in [(128,8),(64,8),(128,4)]]
    d=dict(np.load(OUT/'gr/source-128.npz'))
    knots,co=coefficients(d);trace=co['baryon_g']*LD(d['cx'])*LD(C)**2+co['nonrest_trace_erg']
    direct_source=inspect.getsource(base.charge.independent.direct)
    old="H=flow.green.polynomial(d['t'],trace).antiderivative()";assert direct_source.count(old)==1
    direct_source=direct_source.replace(old,"H=trace_poly.antiderivative()")
    ns=dict(base.charge.independent.direct.__globals__,trace_poly=PPoly(np.asarray(trace[::-1],float),knots))
    exec(compile(direct_source,__file__,'exec'),ns);direct,coordinate=ns['direct'](m,dict(d,t=np.asarray(knots,float)),8)
    controls=dict(quadrature=prior.aligned(low,fine,'U'),independent_GR=abs(direct-rows[0]['endpoint_direct'])/max(abs(direct),1e-290))
    times={key:prior.aligned(coarse,fine,key) for key in ['U','U_t','U_x']}
    old=dict(np.load(Path('native-early-gr207-work')/'gr/fields-128-g8.npz'))
    change={key:float(np.max(abs(fine[key]-old[key]))/max(np.max(abs(fine[key])),1e-290)) for key in ['U','U_t','U_x']}
    row=dict(classification='Counterexample candidate',numerical_controls_passed=controls['quadrature']<.002 and controls['independent_GR']<1e-9,controls=controls,time=times,
        repaired_source_time_representation_applied=True,change_from_straight_line=change,potential_relative=max(float(np.max(abs(fine['potential_U']))/max(np.max(abs(fine['U'])),1e-290)),0),
        source_time_passed=read(OUT/'sources.json')['source_time_passed'],field_time_passed=times['U']<.02,rows=rows,physical_steps=0,GR_return_evolution_executed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    row['GR_return_admitted']=row['numerical_controls_passed'] and row['source_time_passed'] and row['field_time_passed']
    write(OUT/'result.json',row);print(json.dumps(row),flush=True);assert row['numerical_controls_passed'],row

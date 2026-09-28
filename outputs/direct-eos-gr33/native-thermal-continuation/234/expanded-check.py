def check():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    m=prior.prior.owner.Model(64);m.precise_tangent=build();z=np.load(OLD/'rejected-joint-stage.npz');t,h=z['time'][()],z['step'][()]
    sol=z['solution'];v=z['initial'];guides=z['guides'];dim=len(v)
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)];Js=[p[0] for p in maps]
    cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    def L(j,value):
        x,g=m.unpack(value);ph,q,*_=m.collision(cs[j],x,g)
        return m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph,q+(Js[j]@g.ravel()).reshape(m.n,4))
    def mat(value):
        vv=value.reshape(2,dim);return (vv-h*(joint.A@np.array([L(j,row) for j,row in enumerate(vv)]))).ravel()
    raw=joint.LinearOperator((2*dim,)*2,mat,dtype=float);op=dynamics.stable_operator(m,raw,[],True)
    affine=[base-(J@g.ravel()).reshape(m.n,4) for (J,base),g in zip(maps,guides)]
    src=np.array([m.pack(s[0]/(m.scale*joint.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
    rhs=(np.tile(v,(2,1))+h*(joint.A@src)).ravel();rhs=flux.precise_rhs(m,t,h,v,maps,guides,rhs)
    linear=rhs-op.matvec(sol);norm=np.linalg.norm(rhs);rates=[];m.precise_values={}
    for j,row in enumerate(sol.reshape(2,dim)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True)
        native=flux.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    actual=(sol.reshape(2,dim)-v-h*(joint.A@np.array(rates))).ravel();actual=flux.precise_defect(m,t,h,v,sol,actual)
    old_actual=z['defect'];assert read(previous/'result.json')['exact_saved_actual_defect']
    physical=joint.physical_norm(m,linear)/z['physical_scales']
    values=[];nonlinear=[]
    for digits in [40,80]:
        with precision.mp.workdps(digits):
            errors=[]
            for j,(c,guide,row) in enumerate(zip(joint.C,guides,sol.reshape(2,dim))):
                _,g=m.unpack(row);now=t+c*h;delta=precision.hp(g)-precision.hp(guide)
                f0=precision.native_B(m,now,guide,m.precise_tangent);f1=precision.native_B(m,now,g,m.precise_tangent)
                errors.append(f1-f0-precision.baryon_product(Js[j],delta))
                if digits==80:
                    for scale in [precision.mp.mpf('.5'),precision.mp.mpf(2)]:
                        fp=precision.native_B(m,now,precision.hp(guide)+scale*delta,m.precise_tangent)
                        curvature=precision.cast(fp-f0-scale*(f1-f0))
                        nonlinear.append(dict(stage=j,scale=float(scale),maximum_stage_effect=float(abs(h)*np.max(abs(curvature))/norm),norm_stage_effect=float(abs(h)*np.linalg.norm(curvature)/norm)))
            values.append(precision.cast(-precision.hp(h)*(precision.hp(joint.A)@np.array(errors))))
    baryon=np.array([m.unpack(row)[1][:,2] for row in (actual+linear).reshape(2,dim)])
    mismatch=values[-1];change=float(np.linalg.norm(values[0]-mismatch)/norm)
    result=dict(classification='Counterexample candidate',original233_reconstruction_preserved=True,same_saved_proposal=True,
        linear_relative=float(np.linalg.norm(linear)/norm),actual_relative=float(np.linalg.norm(actual)/norm),
        affine_native_B_mismatch=float(np.linalg.norm(mismatch)/norm),decomposition_remainder=float(np.linalg.norm(baryon-mismatch)/norm),precision_change=change,
        increment_linearity=nonlinear,new_physical_steps=0,physical_state_accepted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    result.update(old_actual_relative=float(np.linalg.norm(old_actual)/norm),arithmetic_effect=float(np.linalg.norm(actual-old_actual)/norm),constitutive_B=constitutive(m,t,h,sol),linear_passed=bool(np.linalg.norm(linear)/norm<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(actual)/norm<1e-12 and max(joint.physical_norm(m,actual)/z['physical_scales'])<1e-13))
    assert max(result['constitutive_B'])<.002
    np.savez_compressed(OUT/'comparison.npz',linear=linear,actual=actual,affine_native_B_mismatch=mismatch)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert change<1e-25

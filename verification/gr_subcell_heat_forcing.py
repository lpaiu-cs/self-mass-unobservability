"""Propagate the independently measured heat drive through the initial time matrix."""
from types import ModuleType
from fractions import Fraction as F
import json,subprocess,sys
import numpy as np
import sympy as sp
import gr_subcell_heat_drive_coordinates as heat

ROOT=heat.ROOT;OUT=heat.OUT.parent/'gr-subcell-heat-forcing';sha=heat.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def difference(b,j,h,N,drive,kt,tau):
    """drive=(Q_F-Q)/w; Q-rate output is divided by the fixed initial w."""
    den=1-h/tau-j*j*kt/b
    dv=-N*drive/(tau*den)
    return dv,-dv,-2*j*dv/b,den


def symbolic():
    b,h,N,tau=sp.symbols('b h N tau',positive=True)
    j,drive,kt,kr=sp.symbols('j drive kt kr',real=True)
    dv,dq,dT,den=difference(b,j,h,N,drive,kt,tau)
    matrix=sp.Matrix([[1,0,0,0],[0,b,2*j,0],[0,0,1,1],[-j*kr/2,-j*kt/2,h/tau,1]])
    residual=matrix*sp.Matrix([0,dT,dv,dq])-sp.Matrix([0,0,0,N*drive/tau])
    assert all(sp.simplify(x)==0 for x in residual)
    assert sp.simplify(matrix.det()-b*den)==0
    assert all(sp.simplify(x.subs(drive,0))==0 for x in [dv,dq,dT])
    args=[F(1,2),F(1,100),F(1,4),F(9,10),F(1,1000),F(2),F(1)]
    v,q,T,z=difference(*args)
    B,J,H,NN,D,K,TAU=args
    assert B*T+2*J*v==0 and v+q==0 and q+H/TAU*v-J*K*T/2==NN*D/TAU
    assert v<0<q and T>0 and z>0
    assert difference(*args[:4],F(0),*args[5:])[:3]==(0,0,0)
    return dict(classification='Proven',passed=True,
        variables='Differences between two initial time derivatives at the same v=0 state, fixed rho,T,Q,geometry,composition,opacity and proper tau. b=rho*cv*T/w, j=Q/w, h=K*T/(w*c^2), k_T=2+dlnK/dlnT, drive=(Q_F-Q)/w. The heat-rate difference is delta(Q_dot)/w, not the time derivative of Q/w.',
        equations='delta rho_log_dot=0; b*delta T_log_dot+2*j*delta v_dot=0; delta v_dot+delta Q_dot/w=0; delta Q_dot/w+(h/tau)*delta v_dot-j*k_T*delta T_log_dot/2=N*drive/tau.',
        solution='z=1-h/tau-j^2*k_T/b; delta v_dot=-N*drive/(tau*z), delta Q_dot/w=-delta v_dot, delta T_log_dot=-2*j*delta v_dot/b. The full time-matrix determinant is b*z.',
        singular_boundary='At z=0 and nonzero N*drive/tau the reduced equations are inconsistent. With zero drive they admit a nonzero homogeneous direction. z>0 is only a positive time-coefficient condition; it does not certify the entire characteristic cone or nonlinear well-posedness.',
        scope='This difference formula cancels the common momentum and energy forcing. It needs no replacement of spatial pressure forces or assumption that the original absolute rate is zero. It is not a finite GR trajectory or a recalibration of the heat law.')


def source():
    before=(heat.OUT/'candidate.py').read_text()
    a='rows.append(dict(nodes=number,stored_heat_moment_over_enthalpy=float(old),'
    b='rates=forcing(index,number,rho,cv,w,N,Q,QF,K,T,op,a,weights,c)\n        rows.append(dict(nodes=number,forcing=rates,stored_heat_moment_over_enthalpy=float(old),'
    assert before.count(a)==1
    after=before.replace(a,b);assert after.replace(b,a)==before
    return after,{a:b}


def prepare():
    assert not OUT.exists();heat.verify();OUT.mkdir()
    candidate,changes=source();(OUT/'candidate.py').write_text(candidate)
    times=heat.original.closure.OUT/'constant-times.json'
    plan=json.loads((heat.OUT/'plan.json').read_text())
    files=[ROOT/'verification/gr_subcell_heat_forcing.py',OUT/'candidate.py',heat.OUT/'manifest.json',
        heat.OUT/'repair-manifest.json',times,heat.original.closure.OUT/'symbolic.json']
    plan['bindings'].update({p.relative_to(ROOT).as_posix():sha(p) for p in files})
    plan.update(classification='Counterexample candidate',substitutions=changes,
        checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        constant_proper_seconds=json.loads(times.read_text())['exact_proper_seconds'],
        intervention='Change only the heat-law right-hand side from the former collocated Q_F=Q to the independently differentiated temperature/lapse drive. Reevaluate the identical nodal inputs and opacity branch and require every prior heat-audit field to replay exactly.',
        arithmetic='Freeze long-double b,j,h,N,(Q_F-Q)/w,k_T and enthalpy quadrature weights before exact Fraction algebra. Save all these inputs and floating rate displays. Require exactly zero energy, momentum and heat difference residuals whenever the denominator is nonzero.',
        denominator_gate='Count positive, negative and zero time coefficients separately for both unchanged proper times on every node. No negative or singular node is skipped or rescued by retuning tau. Positivity is not a characteristic-cone certificate.',
        boundary='Initial forcing difference on the stored mixed-precision nonuniform profile. No full absolute GR tangent, finite evolution, continuous native/space error, calibrated transport or observational closure.')
    save('plan.json',plan);save('forcing-symbolic.json',symbolic())
    for name in ['symbolic.json','controls.json']:(OUT/name).write_bytes((heat.OUT/name).read_bytes())


def setup(folder):
    global RAW,TIMES
    RAW=OUT/folder;assert not RAW.exists();RAW.mkdir()
    TIMES=list(map(F,json.loads((OUT/'plan.json').read_text())['constant_proper_seconds']))


def forcing(index,number,rho,cv,w,N,Q,QF,K,T,op,a,weights,c):
    values=np.column_stack([rho*cv/w,Q/w,K*T/(w*c*c),N,(QF-Q)/w,5-op[:,2]]).astype(ld)
    norm=a*w*weights;norm/=np.sum(norm)
    assert np.all(np.isfinite(values)) and np.all(values[:,0]>0)
    inputs=[[F(*x.as_integer_ratio()) for x in row] for row in values]
    rates=np.full((2,number,3),np.nan,dtype=ld);records=[]
    for model,tau in enumerate(TIMES):
        positive=negative=zero=0;minimum=None;max_v=max_T=F(0);maxres=F(0)
        for node,(b,j,h,lapse,drive,kt) in enumerate(inputs):
            z=1-h/tau-j*j*kt/b
            minimum=z if minimum is None else min(minimum,z)
            positive+=z>0;negative+=z<0;zero+=z==0
            if z==0:continue
            dv,dq,dT,den=difference(b,j,h,lapse,drive,kt,tau);assert den==z
            res=max(abs(b*dT+2*j*dv),abs(dv+dq),abs(dq+h/tau*dv-j*kt*dT/2-lapse*drive/tau))
            assert res==0;maxres=max(maxres,res);max_v=max(max_v,abs(dv));max_T=max(max_T,abs(dT))
            rates[model,node]=[ld(x.numerator)/ld(x.denominator) for x in [dv,dq,dT]]
        means=None if zero else (rates[model].T@norm).astype(float).tolist()
        records.append(dict(model=model,positive_nodes=int(positive),negative_nodes=int(negative),singular_nodes=int(zero),
            minimum_scaled_time_coefficient=str(minimum),exact_difference_residual=str(maxres),
            maximum_abs_delta_vdot=str(max_v),maximum_abs_delta_logTdot=str(max_T),
            enthalpy_weighted_rate_differences=means))
    path=RAW/f'cell-{index:04d}-nodes-{number}.npz';assert not path.exists()
    np.savez_compressed(path,coefficients=values,enthalpy_weights=norm,rate_displays=rates)
    return dict(models=records,input_columns=['b','j','h_seconds','lapse','drive_over_w','k_T'],
        rate_columns=['delta_vdot','delta_Qdot_over_w','delta_lnTdot'],raw=path.relative_to(ROOT).as_posix(),sha256=sha(path))


def engine():
    text,_=source();assert text==(OUT/'candidate.py').read_text()
    name='gr_subcell_heat_forcing_candidate';obj=ModuleType(name);obj.__file__=str(OUT/'candidate.py')
    sys.modules[name]=obj;exec(compile(text,obj.__file__,'exec'),obj.__dict__);obj.OUT=OUT;obj.forcing=forcing
    return obj


def audit(rows,prior):
    assert [r['cell'] for r in rows]==[r['cell'] for r in prior]
    models=[dict(positive_nodes=0,negative_nodes=0,singular_nodes=0,maximum_abs_delta_vdot=F(0),
        maximum_abs_delta_logTdot=F(0),minimum_scaled_time_coefficient=None) for _ in range(2)]
    for cell,original in zip(rows,prior,strict=True):
        assert cell['finite_8_16_moment_difference_over_enthalpy']==original['finite_8_16_moment_difference_over_enthalpy']
        for current,old in zip(cell['rows'],original['rows'],strict=True):
            assert {k:v for k,v in current.items() if k!='forcing'}==old
            item=current['forcing'];assert sha(ROOT/item['raw'])==item['sha256']
            for target,record in zip(models,item['models'],strict=True):
                assert record['exact_difference_residual']=='0'
                for k in ['positive_nodes','negative_nodes','singular_nodes']:target[k]+=record[k]
                for k in ['maximum_abs_delta_vdot','maximum_abs_delta_logTdot']:target[k]=max(target[k],F(record[k]))
                z=F(record['minimum_scaled_time_coefficient']);k='minimum_scaled_time_coefficient'
                target[k]=z if target[k] is None else min(target[k],z)
    return dict(classification='Counterexample candidate',cells=len(rows),all_previous_heat_fields_identical=True,
        models=[{k:str(v) if isinstance(v,F) else v for k,v in m.items()} for m in models],
        physical_EOS_certified=False,full_GR_evolution=False)


def preflight():
    setup('raw-preflight');engine().preflight()
    rows=json.loads((OUT/'preflight.json').read_text())['rows'];prior=json.loads((heat.OUT/'preflight.json').read_text())['rows']
    save('forcing-preflight.json',audit(rows,prior))
    files=[OUT/'preflight-manifest.json',OUT/'forcing-preflight.json',OUT/'forcing-symbolic.json']
    save('forcing-preflight-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify_preflight()


def verify_preflight():
    engine().verify_preflight()
    for rel,digest in json.loads((OUT/'forcing-preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for cell in json.loads((OUT/'preflight.json').read_text())['rows']:
        for row in cell['rows']:
            raw=row['forcing'];assert sha(ROOT/raw['raw'])==raw['sha256']
    assert json.loads((OUT/'forcing-symbolic.json').read_text())==symbolic()
    r=json.loads((OUT/'forcing-preflight.json').read_text());assert r['cells']==6 and r['all_previous_heat_fields_identical']
    print('PASS independently driven heat-rate preflight; inspect time-coefficient signs',flush=True)


def run():
    verify_preflight();setup('raw-full');engine().run()
    def rows(path):return [c for p in sorted(path.glob('block-*/result.json')) for c in json.loads(p.read_text())['rows']]
    current=rows(OUT);assert [r['cell'] for r in current]==list(range(5735))
    save('forcing-result.json',audit(current,rows(heat.OUT)))
    files=[OUT/'manifest.json',OUT/'forcing-result.json',OUT/'forcing-preflight-manifest.json']
    save('forcing-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify()


def verify():
    verify_preflight();engine().verify()
    for rel,digest in json.loads((OUT/'forcing-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'forcing-result.json').read_text());assert r['cells']==5735 and r['all_previous_heat_fields_identical']
    assert all(sum(m[k] for k in ['positive_nodes','negative_nodes','singular_nodes'])==137640 for m in r['models'])
    print('PASS full nonuniform initial forcing-difference inventory; not a finite GR trajectory',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

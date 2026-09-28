"""Full fixed-material TOV branches with separately retained shell mass changes."""
import json,sys
from types import FunctionType,MethodType,SimpleNamespace
import numpy as np
import gr_mass_increment_replay as replay

g=replay.g;OUT=g.OUT/'gr-increment-structure'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_increment_structure.py',g.ROOT/'verification/common_eos.py',
        g.ROOT/'verification/baryon_entropy.py',g.ROOT/'verification/audit_structured_enthalpy.py',
        g.OUT/'initial-state-17-4.npz',g.OUT/'reference-state.npz',g.OUT/'initial-structure-17-4.json',
        g.OUT/'initial-adiabats-17.npz',replay.OUT/'plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f161e59',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},subdivisions=[4,8],
        startup_gate='The completed independent three-cell mass-addition replay and shifted-mass comparison must pass before running. Their result SHA is bound separately at run time.',
        equations='Same piecewise-constant nuclear abundances and entropy TOV, same native EOS/root and 17-point Pchip with direct fallback. RK4 evaluates the RHS in binary64 but retains cell shifts and passed face coordinates in longdouble. Every cell mass increment is saved independently of total-mass subtraction.',
        parameters='Reuse the same saved centre pressure, surface radius and surface mass without refitting. Evaluate the complete centre/surface branch match at each subdivision; failure is retained, not relabelled as a new matched star.',
        interface_tolerance=1e-8,finite_endpoint_refinement_tolerance=1e-8,
        seed='Retain the original one-centimetre central seed and analytic outer seed formulas. Add their explicitly computed mass terms to the corresponding integrated shell, and report endpoint conversion/seed mismatch.',
        limits='A complete static constraint recalculation at fixed parameters, not physical evolution. No new full EOS, spatial continuum, primitive or radial error certificate; no replacement of the historical saved initial state or the running time paths.'))


class Structure(g.c.Structure):
    def __init__(self,data,sub):
        init=FunctionType(g.c.Structure.__init__.__code__,dict(g.c.Structure.__init__.__globals__,OUT=g.OUT,EOS=g.EOS))
        init(self,'initial',data,17,sub);self.shells=np.zeros(len(self.mat.lp),dtype=np.longdouble)
        self.stats=dict(calls=0,evaluations=0,maximum_score=0.,label=f'increment-structure-{sub}')
        inverse=FunctionType(g.audit.strict_invert.__code__,dict(g.audit.strict_invert.__globals__,
            ROOT_STATS=self.stats,e=SimpleNamespace(OUT=OUT,save=save)))
        state=FunctionType(g.c.Structure.state.__code__,dict(g.c.Structure.state.__globals__,be=SimpleNamespace(invert=inverse)))
        self.state=MethodType(state,self)

    def step(self,a,b,y,i,B,outer):
        n=max(self.sub,int(np.ceil(abs(b-a)/(.12/self.sub))));h=(b-a)/n
        base=np.asarray(y,dtype=np.longdouble);delta=np.zeros(3,dtype=np.longdouble)
        def rhs(x,z):return self.rhs(x,np.asarray(z,dtype=float),i,B,outer)
        for k in range(n):
            x=a+k*h;z=base+delta
            k1=rhs(x,z);k2=rhs(x+h/2,z+h*k1/2)
            k3=rhs(x+h/2,z+h*k2/2);k4=rhs(x+h,z+h*k3)
            delta+=h*(k1+2*k2+2*k3+k4)/6
        self.shells[i]=delta[1]*np.longdouble(B)*(-1 if outer else 1)
        if i%1000==0:print('INCREMENT STRUCTURE CELL',self.sub,i,flush=True)
        return base+delta


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    replay.verify();assert json.loads((replay.OUT/'result.json').read_text())['repair_direct_mass_passed']
    save('startup-binding.json',dict(classification='Proven',sha256=g.c.sha(replay.OUT/'manifest.json')))
    data=dict(np.load(g.OUT/'reference-state.npz'));old=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    parameters=json.loads((g.OUT/'initial-structure-17-4.json').read_text())['parameters'];rows=[];previous=None
    for sub in plan['subdivisions']:
        solver=Structure(data,sub);m=solver.mat;B=m.B
        error,inner,outer=solver.branches(parameters,record=True)
        pc,rs,ms=parameters;rs=m.R*np.exp(rs);ms=g.c.gr.TARGET*np.exp(ms)
        p,e,b=solver.state(m.lp[0],0);f=1-2*ms/rs
        massB=e/b*np.sqrt(f);lpB=-(e+p)*(ms+4*np.pi*rs**3*p)/(4*np.pi*rs**4*b*np.sqrt(f)*p)
        w0=min(m.dm[0]*1e-6,1e-8/abs(lpB*B));outer_seed=massB*B*w0
        central_seed=np.longdouble(inner[0,2])*np.longdouble(B)
        solver.shells[0]+=np.longdouble(outer_seed);solver.shells[-1]+=central_seed
        assert np.all(solver.shells>0) and np.all(np.isfinite(solver.shells))
        n=len(m.lp);faces=np.zeros((n+1,3),dtype=np.longdouble)
        faces[0]=[rs/m.R,ms/B,m.lp[0]];faces[1:m.split+1]=outer[1:,1:]
        faces[m.split:n]=inner[1:,1:][::-1];faces[n]=[0,0,pc]
        mass=faces[:,1]*np.longdouble(B);radius=faces[:,0]*np.longdouble(m.R)
        integral_total=solver.shells.sum(dtype=np.longdouble)
        record=dict(classification='Counterexample candidate',subdivision=sub,cells=n,
            interface_max=float(abs(error).max()),interface_passed=bool(abs(error).max()<plan['interface_tolerance']),
            interface_components=[float(x) for x in error],shell_mass_integral_geom_m=float(integral_total),
            shell_integral_minus_boundary_mass_geom_m=float(integral_total-np.longdouble(ms)),
            shell_integral_relative_boundary_difference=float(integral_total/np.longdouble(ms)-1),
            first_three_shell_mass_geom_cm=[float(x*100) for x in solver.shells[:3]],
            maximum_radius_relative_initial_difference=float(abs(radius-np.asarray(old['radius_faces_m'],dtype=np.longdouble)).max()/rs),
            maximum_mass_relative_initial_difference=float(abs(mass-np.asarray(old['mass_faces_geom'],dtype=np.longdouble)).max()/ms),
            EOS_lookups=solver.calls,direct_inverse_calls=solver.direct_calls,entropy_root_statistics=solver.stats,
            full_GR_evolution=False,physical_EOS_certified=False,continuous_error_certified=False)
        if previous is not None:
            differences=abs(faces-previous);scale=np.array([1,1,1],dtype=np.longdouble)
            gap=float((differences/scale).max())
            record['finite_endpoint_refinement_max']=gap
            record['finite_endpoint_refinement_passed']=gap<plan['finite_endpoint_refinement_tolerance']
        previous=faces.copy();np.savez_compressed(OUT/f'path-{sub}.npz',faces=faces,radius_m=radius,
            mass_geom_m=mass,shell_mass_geom_m=solver.shells,inner=inner,outer=outer,
            original_baryon_g=data['dm'],original_X=data['X'],original_reference_entropy=data['s_B'])
        rows.append(record);save(f'path-{sub}.json',record);print('INCREMENT STRUCTURE RESULT',record,flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_interfaces_passed=all(r['interface_passed'] for r in rows),
        finite_refinement_passed=rows[-1]['finite_endpoint_refinement_passed'],
        old_initial_state_replaced=False,full_GR_evolution=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(replay.OUT/'manifest.json')==json.loads((OUT/'startup-binding.json').read_text())['sha256']
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['completed']
    for row in result['rows']:
        path=np.load(OUT/f"path-{row['subdivision']}.npz")
        assert len(path['shell_mass_geom_m'])==row['cells'] and np.all(path['shell_mass_geom_m']>0)
    print('PASS complete increment-structure bindings; finite static constraint audit only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

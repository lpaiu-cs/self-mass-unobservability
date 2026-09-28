"""The same atmospheric adiabat in density coordinates; no EOS substitution."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
import def_native_atmosphere as parent

h=parent.h
OUT=parent.OUT/'density-coordinate'


def table():
    assert not OUT.exists();OUT.mkdir()
    data,tab=h.inputs();s0=tab['reference'][0,3];X=data['X'][0]
    previous=json.loads((parent.OUT/'pilot-progress.json').read_text())
    base=previous['rows'][0];rho0=base['raw'][0];lt=base['lnT']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),parent.OUT/'plan.json',parent.OUT/'pilot-progress.json',parent.BACKGROUND/'background-0.001.npz',parent.BACKGROUND/'lapse.npz']},
        correction='The pressure-coordinate predictor entered a state with radiation pressure above total prescribed pressure. Native info 125 is eos_calc offset 100 plus density-underflow branch 25. Keep that failure; solve the same entropy in rho,T coordinates, where pressure is an output.',
        entropy='Same native s and same max(2 erg/g,32 ulp(abs(h_thermal))) energy-unit gate. A recorded pressure-coordinate extended-precision fallback is allowed only at an already physical near-root state and must retain the density within 1e-11.',
        rho_drop_step=.25,maximum_log_density_drop=28,stop_temperature_K=1000,
        budget=dict(hard_timeout_seconds=120,maximum_new_entropy_roots=113,new_evolution_steps=0,automatic_expansion=False),
        gate=dict(entropy_energy_score=1,gamma_consistency=1e-8),
        scope='Native isentropic atmosphere segment. The unresolved cold tail, radiative/thermal closure and dynamic matching remain separate requirements.'))
    eos=h.molecular.model.EOS();stats=dict(calls=0,evaluations=0,maximum_score=0.,label='density-neighbour-fallback')
    fallback=h.molecular.inverse(stats,OUT/'precision-roots');rows=[];evaluations=0;start=time.monotonic()
    def sample(lr,t):
        nonlocal evaluations
        a=eos(2,float(lr),float(t),X);evaluations+=1
        error=a[3]-s0;budget=max(2.,32*np.spacing(abs(a[2]+a[1]/a[0])))
        return abs(np.exp(t)*error)/budget,a,t,error
    slope=.3
    for k in range(113):
        lr=np.log(rho0)-k*.25
        if k:lt-=.25*slope
        best=None
        for _ in range(18):
            trial=sample(lr,lt)
            if best is None or trial[0]<best[0]:best=trial
            if trial[0]<.25:break
            correction=np.clip(np.exp(lt)*trial[3]/trial[1][10],-.15,.15)
            nxt=lt-correction
            if nxt==lt:break
            lt=nxt
        if best[0]>.25:
            lo=hi=best[2]
            for _ in range(3):
                lo=np.nextafter(lo,-np.inf);hi=np.nextafter(hi,np.inf)
                for t in [lo,hi]:
                    trial=sample(lr,t)
                    if trial[0]<best[0]:best=trial
        score,a,lt,_=best
        if score>1:
            af,tf,_=fallback(eos,float(np.log(a[1])),s0,X,lt)
            assert abs(np.log(af[0])-lr)<1e-11
            # Density mode is required for rho-coordinate pressure derivatives.
            score,a,lt,_=sample(float(np.log(af[0])),tf)
        assert score<=1,('Unchanged entropy gate',k,score)
        slope=(a[1]/a[0]-a[9])/a[10];gamma=a[5]+a[6]*slope
        assert abs(gamma/a[4]-1)<1e-8,(k,gamma,a[4])
        rows.append(dict(logrho=float(np.log(a[0])),lnT=float(lt),logP=float(np.log(a[1])),raw=a.tolist(),gamma1=float(gamma),entropy_score=float(score)))
        h.write(OUT/'progress.json',dict(classification='Counterexample candidate',rows=rows,evaluations=evaluations,seconds=time.monotonic()-start))
        if np.exp(lt)<=1000:break
        assert time.monotonic()-start<110,'Table wall budget'
    h.write(OUT/'table.json',dict(classification='Counterexample candidate',passed=True,rows=rows,native_evaluations=evaluations,
        extended_statistics=stats,seconds=time.monotonic()-start,reached_cold_threshold=bool(np.exp(lt)<=1000)))
    print('TABLE',len(rows),'T',np.exp(lt),'P',a[1],'seconds',time.monotonic()-start,flush=True)


def run():
    assert not (OUT/'result.json').exists();d=json.loads((OUT/'table.json').read_text());assert d['passed']
    plan=json.loads((OUT/'plan.json').read_text())
    for f,sha in plan['bindings'].items():assert h.digest(h.ROOT/f)==sha
    h.write(OUT/'integration-plan.json',dict(classification='Counterexample candidate',
        table_sha256=h.digest(OUT/'table.json'),source_sha256=h.digest(Path(__file__)),
        equations='Full static DEF mass, scalar and lapse equations plus Jordan isentropic material pressure. Use rho as independent coordinate and integrate added baryon inventory separately. No fitted cold polytrope.',
        tolerances=[1e-10,1e-12],budget=dict(hard_timeout_seconds=60,native_calls=0,new_evolution_steps=0),
        gates=dict(radius_difference_relative=1e-8,baryon_relative_difference=1e-6,maximum_added_baryon_fraction=1e-12),
        scope='Mechanical atmosphere to the last successful native state. A positive-pressure cold endpoint is not relabelled a vacuum surface.'))
    data,tab=h.inputs();bg=np.load(parent.BACKGROUND/'background-0.001.npz');lap=np.load(parent.BACKGROUND/'lapse.npz')
    body=h.Structure(.001);ys=bg['faces'][0];R=ys[0]*body.R;M=ys[1]*body.B;mu=M/R
    phib=.001*(1+body.mu*ys[3]);vb=.001*body.mu*ys[4]/body.R
    nub=float(lap['nu_faces'][0]);Bstar=body.B
    rows=d['rows'];lr=np.array([r['logrho'] for r in rows]);raw=np.array([r['raw'] for r in rows])
    # Interpolate only the declared native curve, rejecting extrapolation.
    lp=np.array([r['logP'] for r in rows]);lt=np.array([r['lnT'] for r in rows]);u=raw[:,2]*1e-4/h.gr.C**2
    curve=PchipInterpolator(lr[::-1],np.c_[lp,lt,u][::-1],axis=0,extrapolate=False);derivative=curve.derivative()
    cx=data['CX'][0]
    def rhs(x,y):
        logP,logT,u=curve(x);gamma=derivative(x)[0];rho=h.gr.G*np.exp(x)*1000/h.gr.C**2
        p=h.gr.G*np.exp(logP)*.1/h.gr.C**4;en=rho*(cx+u)
        r=R*(1+y[0]);mass=R*(mu+y[1]);phi=phib+.001*y[2];v=vb+.001*y[3]/R
        A=np.exp(body.beta*phi*phi/2);pe,ee=A**4*p,A**4*en;b=1-2*mass/r;N=np.exp(nub+y[4])
        nr=mass/(r*r*b)+4*np.pi*r*pe/b+r*v*v/2
        dr=-gamma*p/(en+p)/(nr+body.beta*phi*v)
        vr=4*np.pi/b*(body.beta*phi*(ee-3*pe)+r*v*(ee-pe))-2*(r-mass)/(r*r*b)*v
        matter=4*np.pi*r*r*ee
        source=4*np.pi*body.beta*(phi/.001)*(ee-3*pe)*N*r*r/(R*np.sqrt(b))
        return dr*np.array([1/R,(matter+r*r*b*v*v/2)/R,v/.001,vr*R/.001,nr,
                            4*np.pi*r*r*A**3*rho/(np.sqrt(b)*Bstar),matter/R,source])
    answers=[];started=time.monotonic();h.mp.mp.dps=70
    for tol in [1e-10,1e-12]:
        answer=solve_ivp(rhs,(lr[0],lr[-1]),np.zeros(8),method='DOP853',rtol=tol,
            atol=[tol*1e-5,1e-27,1e-22,1e-21,1e-22,1e-27,1e-27,1e-26],max_step=.0625,dense_output=True)
        assert answer.success,answer.message
        y=answer.y[:,-1];r=R*(1+y[0]);mass=R*(mu+y[1]);phi=phib+.001*y[2];v=vb+.001*y[3]/R;q=r*v
        ext=h.exterior.exact(h.mp.mpf(mass/r),h.mp.mpf(q));pinf=h.mp.mpf(phi)+h.mp.mpf(q)*ext[2]
        shift=float(h.mp.log1p(-2*h.mp.mpf(mass/r))/2-h.mp.mpf(q)**2*ext[3])-(nub+y[4])
        values=dict(tolerance=tol,radius_m=float(r),thickness_m=float(R*y[0]),added_baryon_fraction=float(y[5]),
            added_material_mass_geom_m=float(R*y[6]),scalar_source_increment_normalized=float(y[7]),
            matched_phi_infinity=str(pinf),interior_lapse_shift=float(shift),RHS_evaluations=answer.nfev,
            endpoint_temperature_K=float(np.exp(lt[-1])),endpoint_pressure_dyn_cm2=float(np.exp(lp[-1])))
        answers.append(values);grid=np.linspace(lr[0],lr[-1],401)
        np.savez_compressed(OUT/f'atmosphere-{tol}.npz',logrho=grid,state=answer.sol(grid),thermo=curve(grid))
    radius_error=abs(answers[1]['radius_m']-answers[0]['radius_m'])/R
    baryon_error=abs(answers[1]['added_baryon_fraction']/answers[0]['added_baryon_fraction']-1)
    passed=radius_error<1e-8 and baryon_error<1e-6 and answers[1]['added_baryon_fraction']<1e-12
    result=dict(classification='Counterexample candidate',mechanical_segment_passed=bool(passed),rows=answers,
        radius_comparison_relative=radius_error,baryon_comparison_relative=baryon_error,
        seconds=time.monotonic()-started,physical_free_surface_completed=False,thermal_stationarity=False,
        dynamic_exterior_coupled=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['table','run']);globals()[p.parse_args().action]()

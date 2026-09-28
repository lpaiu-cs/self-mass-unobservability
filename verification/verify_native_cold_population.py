"""Audit stable arithmetic and finish reporting the saved 49-node solution."""
from decimal import Decimal, localcontext
import json
import signal
import time
import numpy as np
import sympy as s
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_simpson
import def_native_cold_population as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.alarm(25)
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',EOS_calls=190,seconds=25,
        controls='Symbolic exact arithmetic; actual failed operands at80decimal digits; native fixed-inventory first-law at three new cold states with two step sizes; independent saved-path reconstruction.',
        gates=dict(first_law=1e-4,step_contrast=1e-4,entropy_root=2e-10,population=1e-12,coarse_fine_speed=.002),
        reuse='All49 states already persisted before the inherited abs(list) reporting error. Reconstruct the report from that NPZ; do not rerun the adiabat.',
        source_sha256=task.sha(__file__)))
    R,D,A,B,Q=s.symbols('R D A B Q',nonzero=True)
    assert s.simplify(R*(D-A*B/Q)/Q-R*(D/Q-(A/Q)*(B/Q)))==0
    assert s.simplify(R*A/Q-R*(A/Q))==0
    operands=['1.6962531274641092e243','-4.5501281107751468e65','1.2244172427407522e88']
    with localcontext() as ctx:
        ctx.prec=80;mu,dx,q=map(Decimal,operands);exact=mu*dx/q
        vals=list(map(float,operands));stable=vals[0]*(vals[1]/vals[2]);original=vals[0]*vals[1]/vals[2]
        error=float(abs((Decimal(stable)-exact)/exact))
    assert not np.isfinite(original) and np.isfinite(stable) and error<1e-15
    task.write(out/'symbolic.json',dict(classification='Proven',passed=True,
        identity='For nonzero Q, R*(D-A*B/Q)/Q=R*(D/Q-(A/Q)*(B/Q)); mu*x/Q=mu*(x/Q). These are algebraic identities, not uniform floating-point error bounds.'))
    task.write(out/'overflow-regression.json',dict(classification='Counterexample candidate',operands=operands,
        original_finite=False,stable_value=stable,decimal_reference=str(exact),relative_error=error,passed=True))
    d=np.load(out/'adiabat.npz');fail=json.loads((out/'adiabat-failure.json').read_text())
    assert fail['nodes']==49 and fail['error'].startswith('TypeError') and len(d['T'])==49
    assert np.all(np.isfinite(d['raw'])) and min(d['T'])>=100
    assert np.max(abs(d['errors'][:,0]))<2e-10 and np.max(abs(d['errors'][:,1]))<1e-12
    ion=task.StableIons(cap=190);fan=ion.fan;lr0=np.log(fan.rho);target=d['number_fractions'][0]
    base=ion.snapshot(lr0,np.log(fan.T),np.zeros(318));molecules=base['molecular_H_fractions'];checks=[]
    for x in [-4.5,-5.,-6.]:
        i=int(np.flatnonzero(d['log_density_ratio']==x)[0]);lr=lr0+x;lt=np.log(d['T'][i]);fields=d['fields'][i]
        a,_,_=ion.constrain(lr,lt,target,fields,target_molecules=molecules,tolerance=1e-12);a=a['eos']
        derivatives=[];laws=[]
        for h in [2e-4,1e-4]:
            neighbors=[]
            for dr,dt in [(h,0),(-h,0),(0,h),(0,-h)]:
                b,_,_=ion.constrain(lr+dr,lt+dt,target,fields,target_molecules=molecules,tolerance=1e-12)
                neighbors.append(b['eos'])
            ar,br,at,bt=neighbors;radial=(ar-br)/(2*h);thermal=(at-bt)/(2*h)
            laws.append([float((np.exp(lt)*thermal[3]-thermal[2])/thermal[2]),
                         float((np.exp(lt)*radial[3]-radial[2]+a[1]/a[0])/(a[1]/a[0]))])
            gamma=radial[1]/a[1]+thermal[1]/a[1]*(a[1]/a[0]-radial[2])/thermal[2]
            derivatives.append([float(thermal[2]),float(gamma)])
        step=float(max(abs(np.array(derivatives[0])/derivatives[1]-1)))
        checks.append(dict(log_density_ratio=x,T=float(np.exp(lt)),first_law=laws,step_relative=step,
            cvT=derivatives[-1][0],fixed_gamma=derivatives[-1][1],
            passed=bool(np.max(abs(np.array(laws)))<1e-4 and step<1e-4 and min(derivatives[-1])>0)))
    curves=[]
    for stride in [2,1]:
        x=-d['log_density_ratio'][::stride];a=d['raw'][::stride]
        gamma=-PchipInterpolator(x,np.log(a[:,1])).derivative()(x)
        enthalpy=fan.cx*fan.c**2+a[:,2]+a[:,1]/a[:,0]
        cs=np.sqrt(gamma*a[:,1]/a[:,0]/enthalpy)
        rap=cumulative_simpson(cs,x=x,initial=0);curves.append(rap)
        assert np.all(gamma>1)
    assert np.array_equal(curves[1],d['rapidity'])
    contrast=float(max(abs(curves[0]-curves[1][::2]))/max(curves[1]))
    assert contrast<.002 and all(c['passed'] for c in checks)
    eq=np.load(task.old.native.OUT/'fine.npz');comparisons=[]
    for x in [-1.,-2.,-4.,-6.]:
        i=int(np.flatnonzero(d['log_density_ratio']==x)[0]);j=int(np.flatnonzero(eq['log_density_ratio']==x)[0])
        comparisons.append(dict(log_density_ratio=x,fixed_T=float(d['T'][i]),LTE_T=float(eq['T'][j]),
            pressure_ratio=float(d['raw'][i,1]/eq['raw'][j,1]),fixed_H_ion=float(d['raw'][i,14]),LTE_H_ion=float(eq['raw'][j,14])))
    result=dict(classification='Counterexample candidate',passed=True,nodes=49,reused_nodes=34,
        new_EOS_calls=fail['EOS_calls'],new_adiabat_seconds=fail['seconds'],minimum_T=float(min(d['T'])),
        coarse_fine_rapidity_relative=contrast,maximum_entropy_scaled=float(max(abs(d['errors'][:,0]))),
        maximum_population_error=float(max(abs(d['errors'][:,1]))),comparisons=comparisons,
        report_recovered_from_saved_solution=True,original_reporting_failure_preserved=True,
        full_fluid_evolution=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    task.write(out/'adiabat.json',result)
    audit=dict(classification='Counterexample candidate',passed=True,checks=checks,EOS_calls=ion.calls,
        seconds=time.monotonic()-start,source_sha256=task.sha(__file__),
        boundary='The fixed-inventory branch is evaluated, not established as the actual finite-reaction gas. Sampled derivative checks are not a whole-domain interval certificate.')
    task.write(out/'audit.json',audit);ion.save('audit-states.npz');signal.alarm(0)
    print(json.dumps(dict(adiabat=result,audit=audit)),flush=True)


if __name__=='__main__':main()

"""Independent saved-ledger and actual-native controls; same direct readout."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json
import signal
import sys
import time
import numpy as np
import def_native_reactive_flow as flow
import def_native_release_charge as read

ex=flow.exchange
OUT=ex.OUT


def audit():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.alarm(30)
    native=ex.Native(cap=280);rows=[];checks=[];saved={}
    for n in [448,896]:
        d=np.load(flow.OUT/f'cells-{n}.npz');f=flow.Flow(n);m=f.base;U=d['U'];h=d['history'];vol=d['volume'];discard=d['conserved_discard'];saved[n]=d
        denom=np.sum(d['initial'][0]*vol);b=np.sum((U[0]-d['initial'][0])*vol)+discard[0]-h[-1,3]
        en=np.sum((U[2]-d['initial'][2])*vol)+discard[2]-h[-1,4]
        species=np.sum((U[3]-d['initial'][3])*vol)+discard[3]-h[-1,6]-h[-1,7]
        V=f.primitive(U);rho,v,lt,y=V;p,u,g,T,kap,cv,s=f.eos.evaluate(rho,lt)
        recovery=np.max(np.abs(f.conserved(*V,m.a)[0]-U));active=rho>=f.eos.floor
        rows.append(dict(cells=n,baryon=float(abs(b)/denom),energy_absolute=float(en),species=float(abs(species)/(denom*f.eos.y0)),conserved_recovery=float(recovery),neutral_inventory=float(np.sum(U[3]*vol)),net_photon_erg=float(h[-1,5]*4*np.pi*m.RJ**2*f.eos.rho0*ex.C**2)))
        assert rows[-1]['baryon']<1e-10 and rows[-1]['species']<1e-9
        if n==896:
            ids=np.flatnonzero(active)
            for i in ids[np.linspace(0,len(ids)-1,6).astype(int)]:
                x=float(np.log(rho[i]));t=float(lt[i]);yy=float(y[i]);a=native.state(x,t,yy);rr=native.rates(a,float(f.eos.d['Trad']),32)
                f.eos.y=np.array([yy]);r=f.eos.reactions(np.array([rho[i]]),np.array([t]))[0]
                hh=1e-4;ar=(native.state(x+hh,t,yy)['raw']-native.state(x-hh,t,yy)['raw'])/(2*hh);at=(native.state(x,t+hh,yy)['raw']-native.state(x,t-hh,yy)['raw'])/(2*hh)
                gamma=ar[1]/a['raw'][1]+at[1]/a['raw'][1]*(a['raw'][1]/a['raw'][0]-ar[2])/at[2]
                errors=[abs(p[i]*f.eos.rho0*ex.C**2/a['raw'][1]-1),abs(u[i]*ex.C**2/a['raw'][2]-1),abs(g[i]/gamma-1),abs(cv[i]*ex.C**2/at[2]-1)]
                checks.append(dict(cell=int(i),rho=float(rho[i]*f.eos.rho0),T=float(T[i]),y=yy,constitutive=list(map(float,errors)),rates=float(np.max(abs(r/rr-1)))))
    # Abundance convergence is independent of the much larger total mass.
    species_grid=abs(rows[0]['neutral_inventory']/rows[1]['neutral_inventory']-1)
    passed=bool(max(max(q['constitutive']) for q in checks)<.002 and max(q['rates'] for q in checks)<.002 and species_grid<.02)
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,actual_native_controls=checks,neutral_inventory_grid_relative=float(species_grid),EOS_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=ex.old.cold.sha(__file__),physical_chemistry_closed=False,final_charge_solved=False)
    np.savez_compressed(OUT/'audit-native-states.npz',**{k:np.array([z[k] for z in native.ion.states]) for k in native.ion.states[0]})
    ex.write(OUT/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


def adapter(n):
    f=flow.Flow(n);m=f.base;m.eos=f.eos;m.primitive=lambda U:f.primitive(U)[:3]
    return m


def charge():
    out=OUT/'charge';assert not out.exists();out.mkdir();start=time.monotonic();signal.alarm(45)
    flow_passed=json.loads((flow.OUT/'result.json').read_text())['passed']
    assert not flow_passed and json.loads((OUT/'audit.json').read_text())['passed']
    ex.write(out/'plan.json',dict(classification='Counterexample candidate',seconds=45,new_fluid_steps=0,new_native_calls=0,
        claim='Diagnostic readout from saved reactive flows whose integrated-trace refinement failed7.77percent. Compare the direct component at448/896 and against the accepted frozen896 history. A diagnostic wave comparison cannot override that failure.',
        boundary='Retained finite H reaction in a prescribed thin photon bath; local frozen acoustic bulk. Direct component only, no full GR feedback or final physical charge.',
        bindings={str(p):ex.old.cold.sha(p) for p in [Path(__file__),Path(flow.__file__),Path(ex.__file__),OUT/'repaired-bank.npz',flow.OUT/'cells-448.npz',flow.OUT/'cells-896.npz']}))
    raw=np.load(ex.old.OUT/'charge/native-bulk.npz')['raw']
    old=np.load(ex.old.OUT/'charge/cells-896-linear-g12.npz');times=old['u_seconds']
    ns=dict(vars(read),OUT=out,prior=SimpleNamespace(Flow=adapter,OUT=flow.OUT));reader=FunctionType(read.readout.__code__,ns,argdefs=read.readout.__defaults__)
    fine,fr=reader(896,'linear',times,raw);coarse,cr=reader(448,'linear',times,raw)
    error=float(max(abs(fine-coarse))/max(abs(fine)));change=float(max(abs(fine-old['normalized_charge']))/max(abs(old['normalized_charge'])))
    result=dict(classification='Counterexample candidate',passed=False,diagnostic_wave_grid_passed=error<.02,registered_flow_passed=flow_passed,wave_grid_relative=error,relative_change_from_frozen=change,
        reactive_endpoint=float(fine[-1]),frozen_endpoint=float(old['normalized_charge'][-1]),endpoint_components_cm=fr['endpoint_components_cm'],same_sign=bool(fine[-1]*old['normalized_charge'][-1]>0),seconds=time.monotonic()-start,
        finite_H_reactions=True,photon_energy_paired=True,full_photon_transport=False,full_GR_scalar_feedback=False,physical_chemistry_closed=False,final_charge_solved=False,full_goal_complete=False)
    ex.write(out/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def chemical_control():
    assert not (OUT/'chemical-control.json').exists();start=time.monotonic();signal.alarm(20);native=ex.Native(cap=100);rows=[]
    for x,T,y in [(0.,native.fan.T,3.5e-6),(-2.,6000.,2e-6),(-6.,1500.,1e-7)]:
        t=np.log(T);a=native.state(x,t,y);derivatives=[]
        for fraction in [.01,.005]:
            h=y*fraction;plus=native.state(x,t,y+h)['raw'];minus=native.state(x,t,y-h)['raw']
            derivative=((plus[2]-minus[2])-T*(plus[3]-minus[3]))/(2*h)
            derivatives.append(derivative)
        expected=-native.nH*ex.K*T*a['affinity']
        rows.append(dict(x=x,T=T,y=y,affinity=a['affinity'],free_energy_derivative_per_y=derivatives,thermodynamic_expected=expected,relative=abs(derivatives[-1]/expected-1),step_relative=abs((derivatives[-1]-derivatives[0])/expected)))
    passed=bool(max(z['relative'] for z in rows)<.001 and max(z['step_relative'] for z in rows)<.001)
    result=dict(classification='Counterexample candidate',passed=passed,rows=rows,EOS_calls=native.ion.calls,seconds=time.monotonic()-start,
        scope='Fixed rho,T,other ionic inventories and molecule fractions; native (du-Tds)/dy equals -nH*k*T*lambda_H. Checks the conjugate affinity used in the inverse rate, not a microscopic Milne derivation.')
    np.savez_compressed(OUT/'chemical-control-states.npz',**{k:np.array([z[k] for z in native.ion.states]) for k in native.ion.states[0]})
    ex.write(OUT/'chemical-control.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':globals()[sys.argv[1]]()

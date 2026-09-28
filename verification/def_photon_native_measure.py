"""Use the published OPLIB grid, then couple absorption to the same material EOS.

This resolves a reader/quadrature error. It does not certify unresolved atomic
lines, in-medium dispersion, scattering redistribution or a stellar atmosphere.
"""
from pathlib import Path
from fractions import Fraction as F
from types import FunctionType
import argparse
import json
import signal
import time
import urllib.request
import numpy as np
from scipy.integrate import trapezoid
import sympy as sp
import def_photon_loss_bounds as loss
import def_photon_matter_split as matter

native=loss.native;ex=loss.ex;h=loss.h
OUT=native.OUT.parent/'def-photon-native-measure'
GRID_URL='https://arxiv.org/html/1203.5832'
SEGMENTS=[(0,1,800,9600),(12,1,200,1600),(20,1,100,1000),
          (30,1,10,700),(100,1,1,900),(1000,10,1,900),(10000,100,1,200)]


def grid():
    # Frey et al., Table 1. Constants are source data, not fitted grid positions.
    u=np.concatenate([start+np.arange(1,n+1)*num/den for start,num,den,n in SEGMENTS])
    assert len(u)==14900 and np.all(np.diff(u)>0)
    return u


def backward(rhsT,rhsE,Cm,Ci,rates,dt):
    """Arrowhead resolvent (I-dt*A)^-1 without a dense spectral matrix."""
    alpha=1/(1+dt*rates);exchange=1-alpha
    T=(Cm*rhsT+exchange@rhsE)/(Cm+exchange@Ci)
    E=alpha*rhsE+exchange*Ci*T
    return T,E


def evolve(Cm,Ci,rates,duration,steps):
    T=1.;E=np.zeros_like(Ci);gamma=1-1/np.sqrt(2);dt=duration/steps
    energy0=Cm;entropy0=Cm;previous=entropy0;defect=0.;growth=0.
    for _ in range(steps):
        U,V=backward(T,E,Cm,Ci,rates,gamma*dt)
        T,E=backward(T+(1-gamma)/gamma*(U-T),E+(1-gamma)/gamma*(V-E),Cm,Ci,rates,gamma*dt)
        entropy=Cm*T*T+np.sum(np.divide(E*E,Ci,out=np.zeros_like(E),where=Ci>0))
        defect=max(defect,abs(Cm*T+E.sum()-energy0)/energy0)
        growth=max(growth,(entropy-previous)/entropy0);previous=entropy
    return T,E,defect,growth


def norm(T,E,Cm,Ci):
    return np.sqrt(Cm*T*T+np.sum(np.divide(E*E,Ci,out=np.zeros_like(E),where=Ci>0)))


def symbolic():
    C,c1,c2,k1,k2,T,E1,E2=sp.symbols('C c1 c2 k1 k2 T E1 E2',positive=True)
    R1=k1*(c1*T-E1);R2=k2*(c2*T-E2);Td=-(R1+R2)/C
    assert sp.expand(C*Td+R1+R2)==0
    diss=sp.expand(C*T*Td+E1*R1/c1+E2*R2/c2)
    assert sp.factor(diss+k1/c1*(c1*T-E1)**2+k2/c2*(c2*T-E2)**2)==0
    return dict(classification='Proven',passed=True,
        equations='dE_i/dt=gamma_i*(C_i*dT-E_i); C_m*d(dT)/dt=-sum_i dE_i/dt. C_i is the derivative of the fixed-physical-frequency-bin Planck energy.',
        invariants='C_m*dT+sum E_i is conserved. The derivative of (C_m*dT^2+sum E_i^2/C_i)/2 is -sum gamma_i*(C_i*dT-E_i)^2/C_i. Bins with C_i=0 are unexcited in this model.',
        opacity_derivative='At LTE, variation of gamma multiplies B_i-E_i=0, so opacity derivatives are not needed for this absorption/emission Jacobian. This does not justify freezing opacity in a nonlinear evolution.',
        rounded_rates='For relative rate error <=eps<1, the loss-resolvent error divided by its DC value is <=eps/(1-eps), uniformly on the imaginary axis. The energy-response error is <=eps/(2-eps). These bounds concern printed-rate uncertainty in the fixed finite measure only.',
        omitted='Spatial transport, Compton/angular redistribution, atomic model/continuum error and the physical surface are separate.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    original=json.loads((native.OUT/'requests.json').read_text())[2]
    rows=[]
    for label,T in [('cold','.001'),('warm','.003')]:
        name='phGrid'+label
        rows.append(dict(original,name=name,fields=dict(original['fields'],mixname=name,temps=T,datype='cont')))
    ex.write(OUT/'requests.json',rows)
    (OUT/'lanl-tops-form.html').write_bytes((native.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(loss.__file__),Path(native.__file__),Path(matter.__file__),matter.old.BRIDGE,
        matter.old.model.LIB,native.OUT/'plan.json',native.OUT/'result.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    for name in ['surf1000T0015','surf1000T002']:
        paths += [native.OUT/(name+'-table.txt'),native.OUT/(name+'.npz')]
    paths += [matter.model.base.thermal.OUT/'coefficients.npz']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='48d04e7',
        claim='Recover source-grid photon absorption/scattering without fitting opacities to means and connect its absorptive LTE Jacobian to the same matter EOS.',
        reassessment='The two existing temperatures are development data: restoring the known reduced-energy grid removed their Planck integral mismatch. Use the published seven grid segments, validate every printed energy interval, and test two newly requested temperatures before accepting the reader.',
        source_grid_URL=GRID_URL,source_grid_table=1,source_grid_segments=SEGMENTS,
        gates=dict(native_gray_relative=.001,printed_energy_interval_score=1,weight_sum_absolute=1e-11,
            repartition_relative=1e-12,coupled_energy_relative=1e-10,entropy_growth=1e-10,
            coupled_time_order_min=1.8,finest_difference_relative_initial=.001),
        discrete_model='One positive absorption/scattering value per native spectral node. Frequency cells are native-node Voronoi bins, with endpoint-node extension at the low/high ends. Integrate blackbody bin weights directly; no opacity or weight calibration to a native mean. This defines a finite spectral candidate, not a continuum certificate.',
        coupled_model='Homogeneous absorptive LTE tangent at each queried native temperature, exact submitted surface composition/density and native gas-only heat capacity. No scattering gain is invented. Evolve the material plus all spectral photon energies for eight harmonic absorption times with SDIRK2 at 16/32/64 steps.',
        budget=dict(fresh_physical_queries=2,forecast_network_seconds=7,hard_network_seconds=60,per_query_seconds=20,
            offline_hard_seconds=120,native_EOS_calls=4,CPU_workers=1,GPU=False,new_stellar_steps=0,automatic_expansion=False),
        bindings={str(p.relative_to(h.ROOT)) if p.is_relative_to(h.ROOT) else str(p):h.digest(p) for p in paths},
        limits='No actual-temperature interpolation, atomic/EOS population match, line-resolution/continuum certificate, general scattering redistribution, physical atmosphere or GR evolution. The original fine-group additivity failure remains preserved.'))
    ex.write(OUT/'symbolic.json',symbolic())


def fetch():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();signal.alarm(20)
    try:
        with urllib.request.urlopen(GRID_URL,timeout=15) as response:(OUT/'frey2013-grid.html').write_bytes(response.read())
    finally:signal.alarm(0)
    env=dict(vars(native.retrieval),OUT=OUT);env['save']=lambda name,value:ex.write(OUT/name,value)
    fn=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__);results=[]
    for row in json.loads((OUT/'requests.json').read_text()):
        assert not list(OUT.glob(row['name']+'-*'));began=time.monotonic();signal.alarm(20)
        try:result=fn(row)
        except Exception as exc:
            ex.write(OUT/'fetch-failure.json',dict(name=row['name'],error=repr(exc)));raise
        finally:signal.alarm(0)
        result['seconds']=time.monotonic()-began;results.append(result);ex.write(OUT/'progress.json',results)
        print('QUERY',result,flush=True);assert time.monotonic()-start<60
    ex.write(OUT/'retrieval.json',dict(seconds=time.monotonic()-start,queries=2,passed=all(not r['warnings_present'] for r in results)))


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    # Check the downloaded independent grid table before using it.
    from bs4 import BeautifulSoup
    soup=BeautifulSoup((OUT/'frey2013-grid.html').read_bytes(),'html.parser')
    tables=[table for table in soup.find_all('table') if '12.005' in table.get_text() and '10100' in table.get_text()]
    assert len(tables)==1
    table=tables[0].get_text(' ',strip=True)
    for token in ['0.00125','9600','1600','1000','700','900','200','30000']:assert token in table
    u=grid();edges=np.r_[0,(u[:-1]+u[1:])/2,u[-1]+50]
    weights=loss.weights(edges,1,24);fine_weights=loss.weights(edges,1,48)
    assert np.max(abs(weights-fine_weights))<1e-12
    assert max(abs(weights.sum(0)-1))<plan['gates']['weight_sum_absolute']
    header=native.reader.namespace()['header'];checks=[]
    original=json.loads((native.OUT/'requests.json').read_text())[2:]
    requests=[(native.OUT,r,False) for r in original]+[(OUT,r,True) for r in json.loads((OUT/'requests.json').read_text())]
    data,physical=matter.model.base.inputs();gas=matter.old.GasEOS()
    rho=float(json.loads((native.OUT/'requests.json').read_text())[2]['fields']['dens'])
    for folder,row,heldout in requests:
        name=row['name'];path=folder/(name+'-table.txt')
        common=header(path,dict(row,fields=dict(row['fields'],datype='gray')));assert common['input_passed']
        parsed=native.audit.parse(path);token=parsed['temperatures'][0];T=float(token);s=parsed['spectra'][token]
        tokens=parsed['tokens'][token][1];assert len(s)==len(u)==14900
        E_exact=[F(int(round(x*800)),800)*F(row['fields']['temps']) for x in u]
        grid_score=max(native.audit.score(e,t[0]) for e,t in zip(E_exact,tokens));assert grid_score<=1
        assert np.isfinite(s).all() and np.all(s[:,1:]>0) and np.ptp(s[:,2])>0
        for t in tokens:
            lt,ut=native.audit.interval(t[1]);la,ua=native.audit.interval(t[2]);ls,us=native.audit.interval(t[3])
            assert lt<=ua+us and ut>=la+ls
        ka=s[:,2];ks=s[:,3];kt=ka+ks # positive source components; no mean fitting
        means=np.array([weights[:,0]@ka,1/(weights[:,1]@(1/kt))])
        gray=np.array(common['means'][token])[[2,1]];relative=abs(means/gray-1)
        printed_u=s[:,0]/T;p=15/np.pi**4*printed_u**3*np.exp(-printed_u)/(-np.expm1(-printed_u))
        printed_planck=float(trapezoid(p*ka,x=printed_u))
        # Disjoint regrouping of one measure cannot change either moment.
        partition=np.array_split(np.arange(len(u)),37)
        regrouped=np.array([sum(weights[ix,0]@ka[ix] for ix in partition),1/sum(weights[ix,1]@(1/kt[ix]) for ix in partition)])
        regroup_error=float(max(abs(regrouped/means-1)))
        Tk=T*1.602176634e-16/matter.model.base.k
        native_gas=gas(2,float(data['lnd'][0]),float(np.log(Tk)),data['X'][0])
        Cm=float(np.exp(data['lnd'][0])*native_gas[10]/Tk);Ci=4*gas.a_rad*Tk**3*weights[:,1]
        assert Cm>0
        rates=gas.c_light*rho*ka;tau=float(weights[:,1]@(1/rates));duration=8*tau
        runs=[evolve(Cm,Ci,rates,duration,n) for n in [16,32,64]]
        differences=[norm(a[0]-b[0],a[1]-b[1],Cm,Ci) for a,b in zip(runs[:-1],runs[1:])]
        order=float(np.log2(differences[0]/differences[1]));error=float(differences[-1]/np.sqrt(Cm))
        energy=max(r[2] for r in runs);growth=max(r[3] for r in runs)
        # Independent exact two-temperature solution for a constant opacity.
        rate=1/tau;Cr=float(Ci.sum());exact=Cm/(Cm+Cr)+Cr/(Cm+Cr)*np.exp(-rate*(1+Cr/Cm)*duration)
        const=evolve(Cm,Ci,np.full(len(Ci),rate),duration,64)
        const_error=abs(const[0]-exact)
        frequency_residual=0.
        for omega in np.array([.01,1,100])/tau:
            z=1j*omega;fraction=rates/(rates+z)
            temp=1/(z*(Cm+Ci@fraction));photon=Ci*fraction*temp
            frequency_residual=max(frequency_residual,float(abs(z*(Cm*temp+photon.sum())-1)))
        passed=bool(relative.max()<.001 and regroup_error<1e-12 and energy<1e-10 and growth<1e-10 and order>=1.8 and error<.001 and const_error<.001 and frequency_residual<1e-12)
        np.savez_compressed(OUT/(name+'-bank.npz'),u=u,energy_keV=u*T,edges_keV=edges*T,weights=weights,
            absorption=ka,scattering=ks,total=kt,printed_total=s[:,1],rate_absorption=rates,
            rho_neutral_g_cm3=rho,T_K=Tk,C_matter_per_volume=Cm,C_photon_bins=Ci,
            initial_radiation_bins=gas.a_rad*Tk**4*weights[:,0],passed=passed,
            final_dT=runs[-1][0],final_dE=runs[-1][1])
        checks.append(dict(name=name,heldout=heldout,passed=passed,T_keV=T,T_K=Tk,grid_printing_score=grid_score,
            native_gray=gray.tolist(),reconstructed=means.tolist(),gray_relative=relative.tolist(),
            rounded_grid_Planck_relative=abs(printed_planck/gray[0]-1),repartition_relative=regroup_error,
            C_matter=Cm,C_radiation=Cr,harmonic_absorption_seconds=tau,duration_seconds=duration,
            time_steps=[16,32,64],time_order=order,finest_difference_relative_initial=error,
            energy_defect=energy,entropy_growth=growth,constant_opacity_exact_temperature_error=const_error,
            frequency_energy_residual=frequency_residual))
        print('STATE',checks[-1],flush=True)
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in checks),checks=checks,
        seconds=time.monotonic()-start,source_grid_restored=True,finite_spectral_measure_constructed=True,
        absorptive_material_photon_tangent_evolved=True,scattering_redistribution_closed=False,
        continuum_opacity_error_certified=False,actual_temperature_interpolation_certified=False,
        physical_surface_flux_closed=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result['passed'],'SECONDS',result['seconds'],flush=True)
    assert result['seconds']<120


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','fetch','run']);globals()[parser.parse_args().action]()

"""Same-solution frozen-exterior charge including the background photons' launch energy.

Counterexample candidate. The unchanged frozen-exterior reading of the phase-265
coupled solution uses only the Radau reference emission of its high/low
histories. The background photons launched during the drive leave with physical
Killing energy E0*(1+delta_nu+u)(r0,t) instead of E0; phase 261 (incident
metric, high) and phase 264 (returned metric, low) computed this full geometric
response, which is dominated by that launch term. Here the launch increment
L0(t,bin)*(delta_nu+u)(r0,t) is read as an additional emission through the same
frozen exterior operator: its exited energy enters kappa, its arrived energy
epsilon and its scalar mass-debit/stress forcing the exterior charge. The
normalization is the original rational readout. Deflection, measure and
propagation-work terms of the full response and the variation of the exterior
scalar operator are not included; their size is reported, not absorbed.
"""
from pathlib import Path
from decimal import Decimal, localcontext
import json,resource,sys,time
import numpy as np
import propagate_exterior_vacuum as photons
import read_returned_exterior as frozen

p=photons.base
ROOT=Path('native-photon-boundary-ld2-exterior265-work');OUT=ROOT/'complete';FROZEN=ROOT/'full'
ACTUAL=Path('native-photon-boundary-ld2-265-work');GEOM=Path('native-geometric-clock265-work')
HIGH261=photons.OUT;LOW264=Path('native-primitive-photon263-work')
p.OUT=photons.OUT=OUT
read,write,sha,LD=p.read,p.write,p.sha,p.LD
C,G=frozen.C,frozen.G
AW=np.arange(1,8,2,dtype=LD)/32
SETTINGS=[(128,8,8),(128,4,8),(128,8,4),(64,8,8)]
GRID=dict(production=16384,control=8192);SUB=dict(production=16,control=8)
GATES=dict(time=.02,angular=.002,radial=.002,representation=.002,launch_quadrature=1e-6)


def emission(edges,f):
    """Piecewise-linear luminosity through the frozen reader's Radau-node Emission.

    Emission.check() builds a query-by-interval matrix, too large for these
    grids. Check instead that its primitives reproduce the independent exact
    cumulative integrals of the same linear pieces at every edge.
    """
    f=np.asarray(f,LD);lum=np.stack([f[:-1]+(f[1:]-f[:-1])/3,f[1:]],axis=1)
    e=frozen.Emission(edges,lum.reshape(-1,4));h=np.diff(np.asarray(edges)).astype(LD)[:,None]
    first=np.vstack([np.zeros((1,4),LD),np.cumsum(h*(f[:-1]+f[1:])/2,axis=0)])
    second=np.vstack([np.zeros((1,4),LD),np.cumsum(h*first[:-1]+h*h*(2*f[:-1]+f[1:])/6,axis=0)])
    a,b=e.primitives(np.asarray(edges))
    errors=[float(np.max(abs(x-y))/max(np.max(abs(y)),LD('1e-290'))) for x,y in [(a,first),(b,second)]]
    assert max(errors)<1e-12,errors
    return e


def launch_incident(m,t):
    d=m.d;values=[]
    for lo in range(0,len(t),512):
        tt=t[lo:lo+512];n=len(tt)
        nu=m.metric.at(tt,np.full(n,d.r0))['nu'];u0=d.z0['alpha'][0]*m.metric.wave(tt,np.zeros(n))[0]/d.r0
        values.append(nu+u0)
    return np.concatenate(values)


def launch_returned(m,n):
    """(delta_nu+u)(r0) of the applied metric actually used by the coupled solution, on its own clock."""
    met=np.load(ACTUAL/f'metric/metric-{n}-g8.npz');fld=np.load(ACTUAL/f'gr/fields-{n}-g8.npz')
    assert np.array_equal(met['t'],fld['t']) and 'outer_photon_geometric_lapse' in met.files
    return met['t'],met['delta_nu_faces'][:,-1]+m.d.z0['alpha'][0]*fld['U'][:,-1]/m.d.r0


def background(m,t):
    return np.column_stack([np.interp(t,m.clock,np.asarray(v,float)) for v in m.lum.T]).astype(LD)


def build(m,T):
    E={}
    for label,N in GRID.items():
        edges=np.linspace(0,T,N+1);assert all(np.min(abs(edges-k))<=1e-18 for k in m.clock)
        E['high',label]=emission(edges,background(m,edges)*launch_incident(m,edges)[:,None])
    for n in [64,128]:
        clock,value=launch_returned(m,n);knots=np.unique(np.r_[clock,m.clock]);knots=knots[np.r_[True,np.diff(knots)>1e-18]]
        knots[-1]=T;assert abs(clock[-1]-T)<=1e-18
        for label,S in SUB.items():
            x=(knots[:-1,None]+np.diff(knots)[:,None]*np.arange(S)[None,:]/S).ravel();edges=np.r_[x,T]
            E['low',n,label]=emission(edges,background(m,edges)*np.interp(edges,clock,value)[:,None])
    return E


def read_through(ext,kernel,e,times):
    rows=[]
    for now in times:
        h,hh=e.primitives(np.maximum(now-kernel['delay'],0));idx=np.arange(len(h));bins=kernel['node_bins']
        mass=-G/C**3*np.sum(kernel['weights']*kernel['mass']*hh[idx,bins],dtype=LD)
        stress=-G/(2*C**4)*np.sum(kernel['weights']*kernel['stress']*h[idx,bins],dtype=LD)
        h,_=e.primitives(np.maximum(now-kernel['infinity'],0));idx=np.arange(len(h))
        arrived=np.sum(kernel['mw']*kernel['mu']*h[idx,kernel['bins']],dtype=LD)
        exited=e.primitives(np.array([now]))[0][0]@AW
        rows.append([-(mass+stress)/LD(ext.M),arrived,exited])
    return np.array(rows,LD)


def D(v):
    num,den=LD(v).as_integer_ratio();return Decimal(num)/Decimal(den)


def run():
    assert not OUT.exists();OUT.mkdir(parents=True);start=time.monotonic()
    audit=read(FROZEN/'audit.json');fro=read(FROZEN/'result.json');actual=read(ACTUAL/'result.json')
    assert audit['passed'] and fro['passed'] and actual['passed'] and actual['photon_geometric_boundary_applied']
    assert read(GEOM/'result.json')['passed']
    files=[Path(__file__),Path(photons.__file__),Path(p.__file__),Path(frozen.__file__),FROZEN/'result.json',FROZEN/'audit.json',
        ACTUAL/'result.json',GEOM/'result.json',GEOM/'boundary-575.npz',HIGH261/'fine/result.json',HIGH261/'photon-boundary/source.npz',LOW264/'fine/result.json']
    files += [ACTUAL/f'{k}/{k2}-{n}-g8.npz' for n in [64,128] for k,k2 in [('metric','metric'),('gr','fields')]]
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Read the same phase-265 coupled solution through the unchanged frozen exterior and rational mass normalization, adding the background photons\' physical launch-energy increment as an emission of that same operator for its high and low components.',
        method='High: L0(t,bin)*(delta_nu+u)(r0,t) of the incident metric, tabulated on16384/8192uniform intervals (a scratch run measured2.9e-7/1.2e-6against the phase261launch energies with4096/2048). Low: the same product with the applied returned metric of this solution (its outer face delta_nu including the new geometric lapse, and alpha*U/r0 of its fields), piecewise-linear on the metric clock plus background knots with16/8sub-intervals; clocks64/128as the time control. Emission primitives, kernel, arrival and kappa/epsilon/scalar readout are the frozen reader\'s own formulas.',
        gates=GATES,cross_checks='Launch energy against the phase261full propagation at16knots and the575-time output; model-versus-full photon mass (J) difference and phase264low launch energy are reported, not gated.',
        budget=dict(seconds=1800,virtual_GiB=16,CPU_threads=1,new_physical_steps=0),
        limits='Deflection, measure and propagation-work parts of the geometric response and the variation of the exterior scalar operator by the metric perturbation are not included. One GR return; no self-GR, EOS/spatial/uniform or observational closure.',
        bindings={str(q):sha(q) for q in files},full_goal_complete=False))
    p.initialize();m=photons.Photons(8,8)
    T=float(np.load(p.saved.full.source.saved(128))['actual_step_edges'][-1]);assert abs(T-m.d.T)<=1e-18
    E=build(m,T)
    ext=frozen.exterior.Exterior();ext.T=T;ext.t=np.linspace(0,T,17)
    high_source=dict(np.load(frozen.full.OUT/'gr/source-128.npz'))
    assert ext.M==float(high_source['M_cm']) and abs(ext.K/float(high_source['K_cm'])-1)<1e-12
    rows={(r['component'],r['clock'],r['angular'],r['radial']):r for r in fro['rows']}
    with localcontext() as ctx:
        ctx.prec=100;alpha=-D(ext.K)/D(ext.M);fac=D(G)/D(C)**4
    readings={};out_rows=[]
    for n,a,r in SETTINGS:
        kernel=ext.kernel(a,r)
        for component in ['high','low']:
            for label in ['production','control']:
                e=E['high',label] if component=='high' else E['low',n,label]
                values=read_through(ext,kernel,e,ext.t);readings[component,n,a,r,label]=values
            v=readings[component,n,a,r,'production'][-1];row=rows[component,n,a,r]
            with localcontext() as ctx:
                ctx.prec=100
                scalar=Decimal(row['scalar_numerator'])+D(v[0])
                kappa=Decimal(row['kappa'])+fac*D(v[2])/D(ext.M);epsilon=Decimal(row['epsilon'])+fac*D(v[1])/D(ext.M)
                normalized=(scalar+alpha*(epsilon-kappa))/(1+kappa-epsilon)
                frozen_only=Decimal(row['normalized_standalone'])
            item=dict(component=component,clock=n,angular=a,radial=r,compact=row['compact'],
                frozen_exterior=row['exterior'],geometric_exterior=float(v[0]),geometric_arrived_energy_erg=float(v[1]),
                geometric_exited_energy_erg=float(v[2]),scalar_numerator=str(scalar),kappa=str(kappa),epsilon=str(epsilon),
                normalized=str(normalized),frozen_only_normalized=str(frozen_only),
                relative_change_from_frozen_only=float((normalized-frozen_only)/frozen_only))
            out_rows.append(item);write(OUT/f'{component}-{n}-a{a}-r{r}.json',item)
            np.savez_compressed(OUT/f'{component}-{n}-a{a}-r{r}.npz',t=ext.t,geometric_exterior=values[:,0],
                geometric_arrived_energy_erg=values[:,1],geometric_exited_energy_erg=values[:,2])
    table={(v['component'],v['clock'],v['angular'],v['radial']):v for v in out_rows}
    controls={}
    for component in ['high','low']:
        ref=table[component,128,8,8];controls[component]={}
        for name,setting in [('time',(64,8,8)),('angular',(128,4,8)),('radial',(128,8,4))]:
            other=table[(component,)+setting]
            controls[component][name]={k:float(abs(Decimal(other[k])-Decimal(ref[k]))/max(abs(Decimal(ref[k])),Decimal('1e-290')))
                for k in ['normalized','kappa','scalar_numerator']}
    representation={}
    for component in ['high','low']:
        for n,a,r in SETTINGS:
            fine=readings[component,n,a,r,'production'][-1];coarse=readings[component,n,a,r,'control'][-1]
            representation[f'{component}-{n}-a{a}-r{r}']=[float(abs(x-y)/max(abs(x),LD('1e-290'))) for x,y in zip(fine,coarse)]
    components=[]
    with localcontext() as ctx:
        ctx.prec=100
        for n in [64,128]:
            h,l=table['high',n,8,8],table['low',n,8,8]
            qh=Decimal(h['normalized']);den=1+Decimal(h['kappa'])-Decimal(h['epsilon']);dm=Decimal(l['kappa'])-Decimal(l['epsilon'])
            increment=(Decimal(l['scalar_numerator'])-(alpha+qh)*dm)/(den+dm);total=qh+increment
            old=next(v for v in fro['components'] if v['clock']==n)
            components.append(dict(clock=n,high=str(qh),same_solution_low_increment=str(increment),total=str(total),
                frozen_only_total=old['total'],change_from_frozen_only=str(total-Decimal(old['total'])),
                relative_change_from_frozen_only=float((total-Decimal(old['total']))/Decimal(old['total'])),
                sign_negative=bool(total<0),compact_sign_unchanged=bool(total*Decimal.from_float(h['compact'])>0)))
    # Cross-checks against the full propagations (reported).
    e=E['high','production'];rows261=read(HIGH261/'fine/result.json')['rows'];knots=m.clock[1:]
    ours=np.array([e.primitives(np.array([t]))[0][0]@AW for t in knots],LD)
    theirs=np.array([r['physical_launch_energy_increment_erg'] for r in rows261],LD)
    launch_knots=float(np.max(abs(ours-theirs))/np.max(abs(theirs)))
    g=np.load(GEOM/'boundary-575.npz');tt=g['evaluation_t']
    ours575=np.array([e.primitives(np.array([t]))[0][0]@AW for t in tt[1:]],LD);theirs575=np.asarray(g['physical_launch_energy_increment_erg'][1:],LD)
    launch_575=float(np.max(abs(ours575-theirs575))/np.max(abs(theirs575)))
    J=np.asarray(g['photon_geometric_mass_cm'],LD)[1:]/(LD(G)/LD(C)**4)
    measure=float(np.max(abs(J-theirs575))/np.max(abs(J)))
    el=E['low',128,'production'];rows264=read(LOW264/'fine/result.json')['rows']
    ours264=np.array([el.primitives(np.array([t]))[0][0]@AW for t in knots],LD)
    theirs264=np.array([r['physical_launch_energy_increment_erg'] for r in rows264],LD)
    low_launch=float(np.max(abs(ours264-theirs264))/np.max(abs(theirs264)))
    passed=(all(max(v.values())<GATES[name] for group in controls.values() for name,v in group.items())
        and max(max(v) for v in representation.values())<GATES['representation'] and launch_knots<GATES['launch_quadrature'] and launch_575<GATES['launch_quadrature'])
    final=next(v for v in components if v['clock']==128)
    result=dict(classification='Counterexample candidate',passed=bool(passed),controls=controls,representation_controls=representation,
        components=components,rows=out_rows,final_total_clock128=final['total'],final_sign_negative=final['sign_negative'],
        cross_checks=dict(high_launch_energy_vs_phase261_knots=launch_knots,high_launch_energy_vs_575_time_propagation=launch_575,
            full_photon_mass_minus_launch_energy_relative=measure,low_launch_energy_vs_phase264_old_metric=low_launch),
        same_solution=str(ACTUAL),frozen_exterior_operator=True,background_photon_launch_energy_applied=True,
        geometric_photon_mass_in_normalization=True,reciprocal_scalar_forcing_on_frozen_operator=True,
        deflection_measure_work_terms_included=False,exterior_scalar_operator_variation_complete=False,self_GR_return_closed=False,
        physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False,seconds=time.monotonic()-start)
    write(OUT/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);assert passed,result


if __name__=='__main__':
    assert sys.argv[1]=='full';start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(1800)
    try:run()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(ROOT/'complete-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))

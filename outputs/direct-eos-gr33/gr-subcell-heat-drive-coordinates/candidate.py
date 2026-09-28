"""Test the Q_F=Q assumption against the unchanged nonuniform reference profile."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import json,sys,traceback
import numpy as np
import sympy as sp
import gr_metric_coupled_subcell as metric
import gr_heat_entropy_closure as closure
import opacity_tables as opacity
import gr_eos_derivative_coordinates as coordinates

ROOT=metric.ROOT;OUT=metric.OUT.parent/'gr-subcell-heat-drive';sha=metric.sha;ld=np.longdouble


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def derivative_rule(number):
    s,w,_,_=metric.integration_matrix(number)
    gap=s[:,None]-s[None,:];np.fill_diagonal(gap,1)
    barycentric=1/np.prod(gap,axis=1)
    D=barycentric[None,:]/barycentric[:,None]/gap;np.fill_diagonal(D,0)
    D[np.diag_indices(number)]=-np.sum(D,axis=1)
    def differentiate(value):return D@(value-value[0])
    errors=[float(np.max(abs(differentiate(s**k)-(k*s**(k-1) if k else np.zeros_like(s))))) for k in range(number)]
    assert max(errors)<1e-12
    return s,w,D,errors


def symbolic():
    P,rho,cvT,ur,cr,ct,w=sp.symbols('P rho cvT ur cr ct w',positive=True)
    slope=(P/rho-ur)/cvT;gamma=cr+ct*slope;ad=slope/gamma
    assert sp.simplify(ur+cvT*slope-P/rho)==0
    assert sp.simplify((ur+P/rho*(cr-1)+(cvT+P/rho*ct)*slope)/gamma-P/rho)==0
    T,K,a,c,pgrad=sp.symbols('T K a c pgrad',positive=True)
    qf=-K*T/(a*c)*(ad-P/w)*pgrad
    assert sp.simplify(qf/(-K*T/(a*c))-(ad-P/w)*pgrad)==0
    tau,h,Q,N,R,X,Y,V,kr,kt,F=sp.symbols('tau h Q N R X Y V kr kt F',real=True)
    den=w*(tau-h)+tau*Q*kt*V/2
    vdot=(tau*(R-Q*(kr*X+kt*Y)/2)-N*(F-Q))/den
    qdot=R-w*vdot;tlog=Y+V*vdot
    assert sp.simplify(tau*qdot+w*h*vdot-N*(F-Q)-tau*Q*(kr*X+kt*tlog)/2)==0
    return dict(classification='Proven',passed=True,
        Fourier_target='At initial v=0, Q_F=-K*T/(a*c)*partial_r ln(T*N), with Q=proper heat flux/c and conductivity K. The thermal force includes the lapse gradient.',
        isentropic_jet='For a differentiable fixed-composition EOS obeying du=T ds+(P/rho)dlnrho, k=(P/rho-u_lnrho)/(cv*T), Gamma=chi_rho+chi_T*k and nabla_ad=k/Gamma. Along ds=0, dH/dlnP=P/rho, so the enthalpy lapse N*(C+H)=constant gives partial_r ln(T*N)=(nabla_ad-P/(epsilon+P))*partial_r lnP.',
        arbitrary_drive='With the previous frozen initial variables, T_log_dot=Y+V*v_dot and Q_dot+w*v_dot=R. The full quadratic heat law gives v_dot={tau*[R-Q*(kr*X+kt*Y)/2]-N*(Q_F-Q)}/{w*(tau-h)+tau*Q*kt*V/2}. Q_F=Q is a specialization, not a consequence of initially zero material velocity.',
        boundary='The jet identity requires the stated Gibbs/differentiability premises; the numerical EOS outputs do not certify those premises on a continuum. The general rate formula is not a computed GR trajectory.')


def prepare():
    assert not OUT.exists();prior=metric.bindings();closure.verify();model=opacity.Opacity();OUT.mkdir()
    files=[Path(__file__),metric.OUT/'plan.json',closure.OUT/'manifest.json',
        ROOT/'verification/opacity_tables.py',opacity.OUT/'enrichment-control.json',
        opacity.o.OUT/'tables-enriched-tables.npz',opacity.o.OUT/'tables-baseline-trace.json',
        metric.g.OUT/'gr-opacity/new-GR-captured.npz',metric.g.OUT/'gr-opacity/evaluation.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='0ab35ac8',
        bindings=dict(prior['bindings'],**{p.relative_to(ROOT).as_posix():sha(p) for p in files}),
        cells=5735,nodes=[8,16],pilot_cells=[0,1,2,2688,2972,5734],workers=2,block_size=32,
        representation='Keep all stored mixed-precision native EOS nodes, reference radii/metric/pressure, compositions and coordinate-volume weights. Reconstruct the same enthalpy lapse and baryon-affine shared luminosity as the metric inverse. Differentiate the interpolation polynomial on the supplied Gauss nodes, subtracting constants before applying its matrix.',
        conductivity='Reevaluate the existing captured radiative-plus-conductive opacity table at every nodal rho,T; retain the same composition inputs and supported table branches. K=16*sigma*T^3/(3*rho*kappa_eff), with the same stored sigma=5.670400e-5. This effective conductivity is not a photon mean free path or a calibrated causal heat time.',
        controls='Exact polynomial derivative identities through the represented degree, a nontrivial Tolman-temperature positive control, and a lapse-omission negative control. Compare the spectral ln(T*N) derivative with the conditional isentropic EOS-jet derivative and compare integrated forces at 8/16 nodes. Report all discrepancies; do not adjust Q or rescue a failed cell.',
        scope='Finite constitutive-drive audit of the specified subcell reference, including every central/surface cell. No numerical pressure force is inserted into the momentum equation. Piecewise entropy/composition interfaces, native derivative error, continuous chart/Gibbs/gradient error, physical opacity, atmosphere and time evolution remain separate.'))
    save('symbolic.json',symbolic())
    controls=[]
    for n in [8,16]:
        s,_,D,error=derivative_rule(n);nu=s*s/10;lnT=2-nu
        force=D@((lnT-lnT[0])+(nu-nu[0]));without=D@(lnT-lnT[0])
        assert np.max(abs(force))<1e-13 and np.max(abs(without))>.1
        controls.append(dict(nodes=n,polynomial_errors=error,Tolman_force_maximum=float(np.max(abs(force))),
            omitted_lapse_force_maximum=float(np.max(abs(without))),passed=True))
    save('controls.json',dict(classification='Counterexample candidate',passed=True,rows=controls))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def initialize():
    global PLAN,REFERENCE,STATE,GRID,MID,PAR,OPACITY,LUMINOSITY
    PLAN=bindings();REFERENCE=json.loads((metric.OUT/'plan.json').read_text())
    STATE=dict(np.load(metric.g.OUT/'initial-state-17-4.npz'));GRID=dict(np.load(metric.g.OUT/'gr-increment-structure/path-4.npz'))
    MID=np.load(metric.g.OUT/'gr-microphysics/auxiliaries.npz')['eos']
    PAR=np.load(metric.g.OUT/'gr-opacity/new-GR-captured.npz')['parameters'];OPACITY=opacity.Opacity()
    LUMINOSITY=np.r_[0.,np.load(metric.g.OUT/'gr-transport/diagnostics.npz')['interior_Linf'],0.].astype(ld)


def one_cell(index):
    source=ROOT/REFERENCE['cell_sources'][str(index)];rows=[]
    for number in PLAN['nodes']:
        data=metric.node_file(str(source.with_name(source.stem+f'-nodes-{number}.npz')))
        j=list(data['cells']).index(index);order=np.argsort(data['radius_cm'][j])
        r=data['radius_cm'][j,order].astype(ld);eos=data['eos'][j,order].astype(ld)
        Tlog=data['lnT'][j,order].astype(ld);Plog=data['logP'][j,order].astype(ld)
        rho=eos[:,0];P=eos[:,1];T=np.exp(Tlog);a=data['metric_a'][j,order].astype(ld)
        weights=data['coordinate_weights_cm3'][j,order].astype(ld);c=ld(metric.g.c.gr.C)*100
        C=ld(data['C_X'][j])*c*c;H=eos[:,2]+P/rho;w=rho*(C+H)
        Hmid=ld(MID[index,2])+ld(MID[index,1])/np.exp(ld(STATE['lnd'][index]))
        nu_delta=np.log1p((Hmid-H)/(C+H));N=np.exp(ld(STATE['nu'][index])+nu_delta)
        outside=index<len(GRID['outer'])-1;branch=GRID['outer'] if outside else GRID['inner']
        k=index if outside else 5734-index;low=ld(0) if k==0 else ld(branch[k,0]);high=ld(branch[k+1,0])
        knots,_=np.polynomial.legendre.leggauss(number)
        if outside:fraction=(1-knots.astype(ld))/2
        else:
            left=np.cbrt(low);right=np.cbrt(high)
            q=((left+right)/2+(right-left)*knots.astype(ld)/2)**3;fraction=(q-low)/(high-low)
        fraction=fraction[order];assert np.all((fraction>0)&(fraction<1))
        L=LUMINOSITY[index+1]+fraction*(LUMINOSITY[index]-LUMINOSITY[index+1])
        Q=L/(4*np.pi*r*r*N*N*c)
        s,_,D,_=derivative_rule(number)
        def derivative(value):return D@(value-value[0])
        rs=derivative(r);assert np.all(rs>0),(index,number,'nonpositive radius derivative')
        force=(derivative(Tlog)+derivative(nu_delta))/rs
        cr,ct,ur,cv=coordinates.pressure_to_density(eos[:,7],eos[:,8],eos[:,9],eos[:,10])
        slope=(P/rho-ur)/cv;gamma=cr+ct*slope
        assert np.all(gamma>0)
        jet_force=(slope/gamma-P/w)*derivative(Plog)/rs
        par=np.repeat(PAR[index:index+1],number,axis=0).copy()
        par[:,3]=np.log(rho).astype(float)/np.log(10);par[:,4]=Tlog.astype(float)/np.log(10)
        op=np.array([OPACITY(p) for p in par],dtype=ld);K=16*ld('5.670400e-5')*T**3/(3*rho*op[:,0])
        QF=-K*T/(a*c)*force;QJ=-K*T/(a*c)*jet_force
        norm=np.sum(a*w*weights);old=np.sum(a*Q*weights)/norm
        mean=np.sum(a*QF*weights)/norm;jet=np.sum(a*QJ*weights)/norm
        rows.append(dict(nodes=number,stored_heat_moment_over_enthalpy=float(old),
            spectral_Fourier_moment_over_enthalpy=float(mean),isentropic_jet_Fourier_moment_over_enthalpy=float(jet),
            spectral_minus_stored_moment_over_enthalpy=float(mean-old),
            maximum_abs_stored_Q_over_enthalpy=float(np.max(abs(Q/w))),
            maximum_abs_spectral_QF_over_enthalpy=float(np.max(abs(QF/w))),
            maximum_abs_jet_QF_over_enthalpy=float(np.max(abs(QJ/w))),
            maximum_spectral_jet_difference_over_enthalpy=float(np.max(abs(QF-QJ)/w)),
            maximum_constitutive_defect_over_enthalpy=float(np.max(abs(QF-Q)/w)),
            spectral_QF_over_enthalpy=(QF/w).astype(float).tolist(),
            stored_Q_over_enthalpy=(Q/w).astype(float).tolist(),
            opacity_minimum=float(np.min(op[:,0])),opacity_maximum=float(np.max(op[:,0]))))
    diff=abs(rows[0]['spectral_Fourier_moment_over_enthalpy']-rows[1]['spectral_Fourier_moment_over_enthalpy'])
    return dict(cell=index,rows=rows,finite_8_16_moment_difference_over_enthalpy=diff)


def preflight():
    initialize();assert not (OUT/'preflight.json').exists();rows=[];failures=[]
    for i in PLAN['pilot_cells']:
        try:rows.append(one_cell(i));print('HEAT DRIVE PREFLIGHT',i,rows[-1]['finite_8_16_moment_difference_over_enthalpy'],flush=True)
        except Exception:failures.append(dict(cell=i,traceback=traceback.format_exc()))
    save('preflight.json',dict(classification='Counterexample candidate',evaluations_completed=not failures,
        rows=rows,failures=failures,constitutive_equality_assumed=False,full_GR_evolution=False))
    save('preflight-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}))
    assert not failures,failures;verify_preflight()


def verify_preflight():
    plan=bindings()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'preflight.json').read_text());assert r['evaluations_completed'] and [c['cell'] for c in r['rows']]==plan['pilot_cells']
    assert json.loads((OUT/'symbolic.json').read_text())==symbolic()
    print('PASS finite heat-drive preflight evaluation and bindings; inspect discrepancies, not Q_F=Q acceptance',flush=True)


def block(cells):
    folder=OUT/f'block-{cells[0]:04d}';folder.mkdir();rows=[];failures=[]
    for i in cells:
        try:rows.append(one_cell(i))
        except Exception:failures.append(dict(cell=i,traceback=traceback.format_exc()))
    value=dict(cells=cells,rows=rows,failures=failures);path=folder/'result.json';path.write_text(json.dumps(value,indent=2)+'\n')
    return dict(start=cells[0],cells=len(cells),evaluated=len(rows),failures=len(failures),sha256=sha(path))


def run():
    verify_preflight();plan=bindings();assert not any(OUT.glob('block-*'));records=[]
    with ProcessPoolExecutor(max_workers=plan['workers'],initializer=initialize) as pool:
        tasks=[pool.submit(block,list(range(i,min(i+plan['block_size'],plan['cells'])))) for i in range(0,plan['cells'],plan['block_size'])]
        for task in as_completed(tasks):
            records.append(task.result());save('progress.json',dict(blocks=records));print('HEAT DRIVE BLOCK',records[-1],flush=True)
    rows=[r for b in records for r in json.loads((OUT/f"block-{b['start']:04d}"/'result.json').read_text())['rows']]
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=sum(b['cells'] for b in records),
        evaluated=len(rows),failure_count=sum(b['failures'] for b in records),blocks=records,
        maximum_constitutive_defect_over_enthalpy=max(r['maximum_constitutive_defect_over_enthalpy'] for c in rows for r in c['rows']),
        maximum_finite_8_16_moment_difference=max(c['finite_8_16_moment_difference_over_enthalpy'] for c in rows),
        constitutive_equality_assumed=False,physical_transport_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file()}));verify()


def verify():
    verify_preflight()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['cells']==5735
    assert r['evaluated']+r['failure_count']==5735
    print('PASS complete heat-drive inventory bindings; retain every failure and finite discrepancy',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

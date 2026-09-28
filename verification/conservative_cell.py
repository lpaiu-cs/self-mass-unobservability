"""Request28: energy-coordinate nuclear continuation and its mass-readout boundary.

Counterexample candidate: a fixed-volume native cell with a declared energy zero.
This does not certify the physical EOS or evolve a self-gravitating star.
"""
from pathlib import Path
import json, math, shutil, sys
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
import native_closure as native
import fresh_microphysics as fresh
import reactive_energy as reaction
import gr_mass as gr
from thermal_restart import sha
from thermal_wd import mesa

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'outputs/native-closure27'
OUT=ROOT/'outputs/conservative-cell28'
CACHE=Path('/home/lpaiu/work/conservative-cell28')
ZONE=2591
DOCS=['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi']


def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')


def species():
    return json.loads((OLD/'explicit-PP-plan.json').read_text())['species']


def rest():
    iso=reaction.isotopes();names=species()
    q=np.array([float(iso[k]['qex'])*reaction.QCONV/iso[k]['a'] for k in names])
    w=np.array([(float(iso[k]['w'])/iso[k]['a']-1)*(gr.C*100)**2 for k in names])
    return q,w


def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='bf2f600',
        previous_manifest_sha256=sha(OLD/'manifest.json'),zone=ZONE+1,
        duration_proper_days=1.6294,step_counts=[1,2,4],
        energy='Define E=u_native + sum(qex_i Qconv X_i/A_i), up to the conserved baryon constant. '
        'This explicitly chooses the native internal energy as u_Q. In neutral-W notation u_W=u_Q-g(X). '
        'The physical composition-dependent EOS energy zero is not certified by this definition.',
        reaction='All 26 species evolve at fixed density; no intermediate/bulk reservoir. '
        'Exponential Rosenbrock Euler freezes T per reaction step. Temperature is inverted from total energy '
        'with trapezoidal native reaction plus thermal neutrino loss. No heat source is added twice.',
        controls=dict(strict_energy_relative_to_released=1e-6,temperature_log_absolute=2e-6,
            energy_resolution_relative_to_total_released=.05,composition_absolute=1e-16,
            composition_relative_to_total_change=1e-3,maximum_temperature_iterations=5),
        policy='Both strict and resolution-limited diagnostic gates are registered before continuation. '
        'Keep every failed gate. A finite EOS residual is an energy defect, never exact conservation. '
        'If resolution-limited gate fails, stop that continuation and analyze the actual attainable gap.',
        scalar='Derive fixed-metric static susceptibility first variation and weak-field Q/M cancellation; '
        'validate the former on the frozen GR interpolant with finite perturbations and a uniform sphere.',
        exclusions='No runtime rebuild/install, full stellar dynamics, scalar drive assignment, or new observations.'))
    history={}
    dest=OUT/'request27-notes';dest.mkdir()
    for rel in ['docs/'+s+'.md' for s in DOCS]+['paper/revision-manifest.json']:
        expected=json.loads((OLD/'manifest.json').read_text())['sha256'][rel]
        assert sha(ROOT/rel)==expected
        p=dest/Path(rel).name;shutil.copy2(ROOT/rel,p)
        history[rel]=dict(snapshot=p.relative_to(ROOT).as_posix(),sha256=expected,historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',history)


def context():
    native.OUT=OUT;native.CACHE=CACHE;native.context()


def evaluate(label,X,lnT):
    context();base=dict(np.load(OLD/'PP-steady-2-input.npz'))
    base['X']=base['X'].copy();base['lnT']=base['lnT'].copy()
    base['X'][ZONE]=X;base['lnT'][ZONE]=lnT
    path=OUT/(label+'-native.npz')
    if not path.exists():
        extra=' use_eosDT_HELMEOS=.true.' if label.startswith('HELM-') else ''
        native.setup(label,base,extra=extra,species=species(),network=(OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        native.trace(label)
    else:
        inp=np.load(OUT/(label+'-input.npz'))
        assert np.array_equal(inp['X'],base['X']) and np.array_equal(inp['lnT'],base['lnT'])
    d={k:v[ZONE] for k,v in np.load(path).items()}
    _,profile=mesa(OUT/(label+'-profile.data.gz'))
    p={k:float(v[ZONE]) for k,v in profile.items()}
    assert max(abs(d['X']-X))<1e-12 and abs(np.log(d['T'])-lnT)<1e-12
    assert d['heat']==p['eps_nuc'] and d['neutrino']==p['eps_nuc_neu_total']
    return d,p


def initial():
    d={k:v[ZONE] for k,v in np.load(OLD/'PP-steady-2-native.npz').items()}
    _,profile=mesa(OLD/'PP-steady-2-profile.data.gz')
    return d,{k:float(v[ZONE]) for k,v in profile.items()}


def ledger():
    q,w=rest();d,p=initial();f=d['dxdt'];extra=[species().index(k) for k in ['h2','li7','be7','b8']]
    he4=species().index('he4');closed_four=np.zeros(26);closed_four[extra]=f[extra];closed_four[he4]=-sum(f[extra])
    qmass=-q@f-d['neutrino'];correction=d['heat']-qmass
    rows=[]
    for h in [1e-4,5e-5]:
        _,a=mesa(OLD/f'PP-T-{h}--1-profile.data.gz');_,b=mesa(OLD/f'PP-T-{h}-1-profile.data.gz')
        fd=(b['energy'][ZONE]-a['energy'][ZONE])/(2*h)
        rows.append(dict(h=h,energy_dlnT=fd,relative_difference_from_native_cvT=fd/(p['cv']*d['T'])-1))
    save('energy-ledger.json',dict(classification='Counterexample candidate',zone=ZONE+1,
        heat_native_erg_g_s=float(d['heat']),reaction_neutrino_erg_g_s=float(d['neutrino']),thermal_neutrino_erg_g_s=p['non_nuc_neu'],
        Q_rest_release_erg_g_s=float(-q@f),neutral_W_rest_release_erg_g_s=float(-w@f),
        energy_reference_shift_rate_erg_g_s=float((w-q)@f),
        Q_conservative_deposition_erg_g_s=float(qmass),native_minus_Q_conservative_erg_g_s=float(correction),
        correction_is_independent_physical_weak_work_measurement=False,
        four_state_closed_rest_release_erg_g_s=float(-q@closed_four),
        bulk_plus_reservoir_rest_release_erg_g_s=float(-q@(f-closed_four)),
        cvT_erg_g=float(p['cv']*d['T']),temperature_derivative_checks=rows,
        conclusion='The four-state intermediate steady subsystem is not a closed burning cell: bulk fuel supplies its nonzero heat. '
        'The native-minus-Q difference is measured, not assumed to be an independently certified electron-work correction. '
        'Native EOS energy derivatives are not promoted to exact derivatives.'))
    print('ENERGY LEDGER',float(-q@f),float(correction),float(-q@closed_four),rows,flush=True)


def evolve(helm=False,coupled=False):
    plan=json.loads((OUT/'plan.json').read_text());gate=plan['controls'];q,w=rest();d0,p0=initial()
    prefix='HELM-coupled-' if coupled else 'HELM-' if helm else ''
    if helm:
        assert json.loads((OUT/'HELM-plan.json').read_text())['default_failure_sha256']==sha(OUT/'continuation.json')
        d0,p0=evaluate('HELM-initial',d0['X'],float(np.log(d0['T'])))
    X0=d0['X'].copy();lnT0=float(np.log(d0['T']));duration=plan['duration_proper_days']*86400
    total_release=duration*abs(q@d0['dxdt']);summaries=[];previous=None
    if coupled:
        assert json.loads((OUT/'coupled-plan.json').read_text())['split_failure_sha256']==sha(OUT/'HELM-continuation.json')
        thermo=json.loads((OUT/'HELM-composition-energy.json').read_text());assert thermo['passed']
        uX=np.array(thermo['u_X'])
    for steps in plan['step_counts']:
        dt=duration/steps;change=np.zeros(26,dtype=np.longdouble);X=X0.copy();lnT=lnT0;d,p=d0,p0
        lost=0.;records=[];complete=True
        for step in range(steps):
            loss_increment=None
            if coupled:
                b=d['dxdt_T']*d['T'];capacity=p['cv']*d['T'];nuX=-q@d['jacobian']-d['heat_X'];nuT=-q@b-d['heat_T']*d['T']
                block=np.zeros((29,29));block[:26,:26]=d['jacobian'];block[:26,26]=b;block[:26,-1]=d['dxdt']
                block[26,:26]=(d['heat_X']-uX@d['jacobian'])/capacity
                block[26,26]=(d['heat_T']*d['T']-uX@b)/capacity
                block[26,-1]=(-(q+uX)@d['dxdt']-d['neutrino']-p['non_nuc_neu'])/capacity
                block[27,:26]=nuX;block[27,26]=nuT;block[27,-1]=d['neutrino']+p['non_nuc_neu']
                increment=expm(dt*block)[:28,-1];delta=increment[:26];loss_increment=float(increment[27])
                assert loss_increment>=0
            else:
                block=np.zeros((27,27));block[:26,:26]=d['jacobian'];block[:26,-1]=d['dxdt']
                delta=expm(dt*block)[:26,-1]
            change+=delta.astype(np.longdouble)
            proposed=np.asarray(X0.astype(np.longdouble)+change,dtype=float)
            negative=float(min(proposed));assert negative>=-1e-24,negative
            proposed=np.maximum(0,proposed);ix=int(proposed.argmax());proposed[ix]+=1-sum(proposed)
            projection=proposed-(X0.astype(np.longdouble)+change)
            change=proposed.astype(np.longdouble)-X0.astype(np.longdouble)
            rest_change=float(q@(proposed-X0));nu_before=float(d['neutrino']+p['non_nuc_neu'])
            # ponytail: use native cv only as an inversion slope; actual energy residual decides acceptance.
            guess=lnT+float(increment[26]) if coupled else lnT+(-float(q@(proposed-X))-dt*nu_before)/(p['cv']*d['T'])
            samples=[];best=None
            for iteration in range(gate['maximum_temperature_iterations']):
                label=f'{prefix}cell-{steps}-{step}-{iteration}';nd,npf=evaluate(label,proposed,guess)
                nloss=lost+(loss_increment if coupled else dt*.5*(nu_before+float(nd['neutrino']+npf['non_nuc_neu'])))
                defect=(npf['energy']-p0['energy'])+rest_change+nloss
                row=dict(label=label,lnT=guess,defect_erg_g=float(defect),total_loss_erg_g=float(nloss))
                samples.append(row)
                if best is None or abs(defect)<abs(best[0]): best=(defect,guess,nd,npf,nloss,label)
                if abs(defect)<=gate['strict_energy_relative_to_released']*total_release: break
                guess-=defect/(npf['cv']*nd['T'])
            defect,lnT,d,p,lost,label=best;X=proposed
            score=abs(defect)/total_release;ok=score<gate['energy_resolution_relative_to_total_released']
            records.append(dict(step=step+1,selected=label,samples=samples,energy_score=score,
                strict_energy_passed=bool(score<gate['strict_energy_relative_to_released']),resolution_energy_passed=bool(ok),
                max_projection=float(np.max(abs(projection))),minimum_before_projection=negative))
            print('CELL',steps,step+1,'dlnT',lnT-lnT0,'defect/release',score,flush=True)
            if not ok: complete=False;break
        row=dict(steps=steps,completed=complete,records=records,delta_lnT=lnT-lnT0,
            max_delta_X=float(max(abs(X-X0))),energy_defect_erg_g=float(defect),
            energy_defect_over_total_initial_release=score,neutrino_loss_erg_g=float(lost),
            Q_rest_change_erg_g=float(q@(X-X0)),internal_energy_change_erg_g=p['energy']-p0['energy'],
            strict_energy_passed=all(r['strict_energy_passed'] for r in records))
        np.savez_compressed(OUT/f'{prefix}cell-{steps}-end.npz',X=X,lnT=lnT,energy=p['energy'],neutrino_loss=lost)
        if complete and previous is not None:
            oldX,oldT,oldLoss=previous;scoreX=float(max(abs(X-oldX)/(gate['composition_absolute']+gate['composition_relative_to_total_change']*abs(X-X0))))
            row['refinement']=dict(composition_score=scoreX,lnT_difference=abs(lnT-oldT),
                passed=bool(scoreX<1 and abs(lnT-oldT)<gate['temperature_log_absolute'] and (not coupled or abs(lost-oldLoss)/max(lost,1e-30)<1e-3)))
            if coupled: row['refinement']['neutrino_loss_relative_difference']=abs(lost-oldLoss)/max(lost,1e-30)
        previous=(X,lnT,lost) if complete else None;summaries.append(row)
        save(prefix+'continuation.json',dict(classification='Counterexample candidate',rows=summaries,
            total_initial_rest_release_scale_erg_g=float(total_release),exact_energy_conservation=False,
            physical_EOS_energy_zero_certified=False,full_GR_evolution=False))


def helm_plan():
    save('HELM-plan.json',dict(classification='Counterexample candidate',
        default_failure_sha256=sha(OUT/'continuation.json'),
        model='Existing use_eosDT_HELMEOS=.true. with the same 26-species network. '
        'This changes the EOS to the fully ionized approximation; it is not a repair or certification of the default blended evaluator. '
        'The original peak rho, T and all abundances are retained at initialization. Reevaluate all EOS-dependent rates.',
        gates='Keep the original strict energy and time-refinement gates unchanged. All original failures remain.',
        domain='Only the peak hot fixed-volume cell is interpreted. Other profile zones are a native-call carrier. '
        'No cold atmosphere, hydrostatic equilibrium, global Helmholtz derivative consistency, or real-star EOS calibration is asserted.'))


def helm(): evolve(True)


def coupled_plan():
    save('coupled-plan.json',dict(classification='Counterexample candidate',
        split_failure_sha256=sha(OUT/'HELM-continuation.json'),
        changes='Couple all 26 reaction variables, delta lnT and an integrated neutrino-energy variable in one exponential step. '
        'Project endpoint T using the actual HELM energy. This resolves the fast linear nuclear transients in the loss quadrature.',
        Jacobian='Use native f_X,f_T,heat_X,heat_T; use the source identity nu=-q*f-heat for an approximate reaction-neutrino Jacobian. '
        'Use a checked initial HELM energy composition gradient. Heat-capacity derivatives, thermal-neutrino derivatives, '
        'weak-Q state corrections, and subsequent composition-gradient changes are omitted from this approximate step matrix; '
        'the base vector and endpoint EOS/nu are evaluated natively. This is not a full physical Jacobian certificate.',
        controls='Same strict energy and composition/temperature refinement gates; additionally require integrated neutrino step-refinement <1e-3. '
        'Check an analytic irreversible conversion with finite heat capacity and known neutrino fraction.',
        HELM_gradient='At fixed rho,T, HELM internal energy depends on abar,zbar. Use independent H1/He4 and C12/He4 directions '
        'at absolute steps 1e-6 and 5e-7, check refinement <1e-4; reconstruct its abundance gradient up to a baryon constant.'))


def helm_gradient():
    d,p=initial();d,p=evaluate('HELM-initial',d['X'],float(np.log(d['T'])))
    X=d['X'];names=species();iso=reaction.isotopes();A=np.array([iso[k]['a'] for k in names]);Z=np.array([iso[k]['z'] for k in names])
    number=sum(X/A);electron=sum(X*Z/A);abar=1/number;zbar=electron/number
    da=-abar**2/A;dz=(Z-zbar)*abar/A;he=names.index('he4');directions=[names.index(k) for k in ['h1','c12']];rows=[]
    for h in [1e-6,5e-7]:
        rhs=[]
        for j in directions:
            values=[]
            for sign in [-1,1]:
                xp=X.copy();xp[j]+=sign*h;xp[he]-=sign*h
                _,pp=evaluate(f'HELM-gradient-{names[j]}-{h}-{sign}',xp,float(np.log(d['T'])))
                values.append(pp['energy'])
            rhs.append((values[1]-values[0])/(2*h))
        matrix=np.array([[da[j]-da[he],dz[j]-dz[he]] for j in directions]);slopes=np.linalg.solve(matrix,rhs)
        grad=slopes[0]*da+slopes[1]*dz;grad-=grad[he]
        rows.append(dict(h=h,slopes=slopes.tolist(),u_X=grad.tolist()))
    a,b=[np.array(row['u_X']) for row in rows];err=float(max(abs(a-b))/max(abs(b)));assert err<1e-4
    save('HELM-composition-energy.json',dict(classification='Counterexample candidate',rows=rows,u_X=b.tolist(),
        relative_refinement=err,passed=True,scope='Initial fixed-rho,T fully ionized HELM composition gradient, baryon-constant gauge u_He4=0.'))
    # Energy-conserving analytic reaction X0 -> X1 with fixed heat capacity.
    lam=.3;Q=7.;fraction=.2;cap=5.;t=12.
    matrix=np.zeros((5,5));matrix[0,0]=-lam;matrix[1,0]=lam;matrix[2,0]=(1-fraction)*Q*lam/cap;matrix[3,0]=fraction*Q*lam
    got=expm(t*matrix)@np.array([1.,0,0,0,1]);converted=1-np.exp(-lam*t)
    expected=np.array([1-converted,converted,(1-fraction)*Q*converted/cap,fraction*Q*converted,1])
    assert max(abs(got-expected))<1e-13 and abs(Q*got[0]+cap*got[2]+got[3]-Q)<1e-13
    save('coupled-positive-control.json',dict(classification='Proven',passed=True,max_absolute_error=float(max(abs(got-expected))),
        test='Irreversible fuel conversion, temperature rise, integrated neutrino loss and exact total-energy invariant.'))
    print('HELM GRADIENT',err,'PASS coupled analytic control',flush=True)


def coupled(): evolve(True,True)


def symbolic():
    import sympy as s
    q0,q1,w0,w1,x0,x1,u,nloss=s.symbols('q0 q1 w0 w1 x0 x1 u nloss')
    Rq=q0*x0+q1*x1;Rw=w0*x0+w1*x1;g=Rw-Rq
    assert s.expand(Rw+(u-g)-(Rq+u))==0
    f0,f1,nu=s.symbols('f0 f1 nu');udot=-q0*f0-q1*f1-nu
    assert s.expand(q0*f0+q1*f1+udot+nu)==0
    # Normalized monopole: universal zero-binding charge is alpha*M.
    a,M,E,eps=s.symbols('a M E eps',nonzero=True)
    assert s.diff(a*(M+eps*E)/(M+eps*E),eps)==0
    # Wronskian variation: (p phi')'=v phi and delta phi(infinity)=0.
    p,dp,phi,phip,phipp,z,zp,zpp,v,dv=s.symbols('p dp phi phip phipp z zp zpp v dv')
    divergence=s.expand(dp*(phi*zp-z*phip)+p*(phi*zpp-z*phipp))
    assert s.simplify(divergence.subs({phipp:(v*phi-dp*phip)/p,zpp:(v*z+dv*phi-dp*zp)/p})-dv*phi**2)==0
    save('symbolic.json',dict(classification='Proven',energy_reference_invariant=True,
        closed_cell_rest_plus_internal_source_cancels=True,universal_charge_mass_ratio_cancellation=True,
        fixed_metric_susceptibility_variation='delta chi = - integral phi^2 delta v dr, phi(infinity)=1, exterior phi=1+chi/r+O(r^-2)',
        assumptions='Smooth static fixed metric, regular center, compactly supported delta v, fixed outer mass and asymptotic normalization. '
        'General stellar perturbations also change the metric, radius, ADM mass and matter; this is not their full derivative.',
        no_go='In a closed system internal nuclear conversion is not an independent source of total mass. '
        'Escaping radiation, external work and self-gravitating structure must be included. '
        'A universal leading charge Q=alpha*M gives zero variation of Q/M under an internal redistribution.'))
    print('PASS symbolic energy, normalization, and susceptibility variation',flush=True)


def scalar():
    raw=np.load(ROOT/'outputs/remaining-closure26/scalar-background.npz')
    r=raw['r'];R=r[-1];M=raw['mass'][-1];x=r/R
    N=PchipInterpolator(x,np.exp(raw['nu']));f=PchipInterpolator(x,1-2*np.divide(raw['mass'],r,out=np.zeros_like(r),where=r>0))
    tr=PchipInterpolator(x,raw['trace']);mu=M/R;H=-np.log1p(-2*mu)/(2*mu)
    beta=-4.;v=lambda t:4*np.pi*gr.G*beta*float(tr(t))*R*R/gr.C**4*float(N(t))/np.sqrt(float(f(t)))*t*t
    # Fixed compact-support shape on the actual GR interpolant; not a stellar perturbation.
    shape=lambda t:16*t*t*(1-t)**2
    def solve(epsilon,tol):
        start=1e-6;aa=4*np.pi*gr.G*beta*float(tr(0))*R*R/gr.C**4
        def rhs(t,y):
            p=t*t*float(N(t))*np.sqrt(float(f(t)))
            return [y[1]/p,(1+epsilon*shape(t))*v(t)*y[0],v(t)*shape(t)*y[0]**2]
        sol=solve_ivp(rhs,(start,1),[1+aa*start**2/6,float(N(0))*aa*start**3/3,0],
            method='DOP853',rtol=tol,atol=[tol*.01,tol*1e-6,tol*1e-8],max_step=.001)
        assert sol.success;a,j,I=sol.y[:,-1];normal=a+H*j
        return -R*j/normal,-R*I/normal**2
    rows=[]
    for tol in [1e-10,2e-12]:
        chi,derivative=solve(0,tol)
        for h in [1e-3,5e-4]:
            plus=solve(h,tol)[0];minus=solve(-h,tol)[0];fd=(plus-minus)/(2*h)
            rows.append(dict(tol=tol,h=h,chi_m=chi,variational_derivative_m=derivative,
                finite_difference_m=fd,relative_error=abs(fd/derivative-1)))
    assert max(row['relative_error'] for row in rows)<1e-4
    # Independent flat uniform sphere: chi=tan(k)/k-1, v=-k^2*r^2.
    k=.04;eps=1e-5
    chi=lambda z:np.tan(k*np.sqrt(1+z))/(k*np.sqrt(1+z))-1
    analytic=(1/np.cos(k)**2-np.tan(k)/k)/2
    from scipy.integrate import quad
    integral=quad(lambda t:k*k*t*t*(np.sinc(k*t/np.pi)/np.cos(k))**2,0,1,epsabs=1e-15)[0]
    assert abs(integral/analytic-1)<1e-10
    assert abs((chi(eps)-chi(-eps))/(2*eps)/analytic-1)<1e-6
    save('scalar-variation.json',dict(classification='Counterexample candidate',rows=rows,
        uniform_sphere_relative_error=abs(integral/analytic-1),passed=True,
        full_mass_normalized_stellar_derivative=False,
        meaning='Actual fixed-GR-interpolant potential derivative verified independently. '
        'It cannot be multiplied by the PP heat pole without deriving scalar drive, energy/pressure perturbations and metric response.'))
    print('SCALAR VARIATION',rows,flush=True)


def bounds():
    from fractions import Fraction as F
    pi=F(22,7);G=F('6.67428e-11');rho=F(200000000);R=F(70000000);c=F(299792458)
    b=F(4);floor=F('0.99');radius=F('0.01')
    kappa=2*pi*G*b*rho*R**2/(c*c*floor**2)
    B=4*pi*G*b*rho*R**3/(3*c*c*floor)
    q=(1+radius)*kappa;assert q<1
    raw=np.load(ROOT/'outputs/remaining-closure26/scalar-background.npz');r=raw['r']
    assert r[-1]<=float(R) and np.max(raw['trace'])<=float(rho*c*c) and np.min(raw['trace'])>=0
    assert min(np.exp(raw['nu']))>=float(floor) and max(np.exp(raw['nu']))<=1
    assert min(1-2*np.divide(raw['mass'],r,out=np.zeros_like(r),where=r>0))>=float(floor)
    derivative_bounds={str(n):F(math.factorial(n))*B*kappa**(n-1)/(1-q)**(n+1) for n in [1,2,3]}
    save('derivative-bounds.json',dict(classification='Proven',kappa=str(kappa),kappa_decimal=float(kappa),
        B_m=str(B),parameter_radius=str(radius),
        bounds_exact_m={k:str(v) for k,v in derivative_bounds.items()},
        bounds_decimal_m={k:float(v) for k,v in derivative_bounds.items()},
        central_difference_truncation_bound_m={str(h):float(derivative_bounds['3']*F(str(h))**2/6) for h in [1e-3,5e-4]},
        assumptions='Fixed static spherical metric globally N>=0.99, f>=0.99, N<=1; matter radius<=7e7m; '
        '0<=trace<=2e8*c^2 SI, beta=-4, v(epsilon)=v0*(1+epsilon*s), |s|<=1, |epsilon|<=0.01. '
        'Schwarzschild exterior and regular center; derivative with respect to this potential parameter only.',
        proof='The positive static Green operator K has sup norm <=kappa and the asymptotic-charge functional L norm <=B. '
        'phi=(I-K0-epsilon*K1)^(-1)1; ||K1||<=kappa, ||L1||<=B. Differentiate its convergent Neumann resolvent. '
        '|chi^(n)| <= n! B kappa^(n-1)/(1-(1+a)kappa)^(n+1) for n>=1 and |epsilon|<=a.',
        exclusions='Bounds do not include ODE solver error or derivatives of EOS/TOV reconstruction, variable metric, or observations.'))
    print('PASS exact rational shape-parameter derivative bounds',float(kappa),flush=True)


def finite_variation():
    raw=np.load(ROOT/'outputs/remaining-closure26/scalar-background.npz');r=raw['r'];R=r[-1];M=raw['mass'][-1];x=r/R
    N=PchipInterpolator(x,np.exp(raw['nu']));f=PchipInterpolator(x,1-2*np.divide(raw['mass'],r,out=np.zeros_like(r),where=r>0))
    tr=PchipInterpolator(x,raw['trace']);mu=M/R;H=-np.log1p(-2*mu)/(2*mu);rows=[]
    bound=json.loads((OUT/'derivative-bounds.json').read_text())['central_difference_truncation_bound_m']
    for tol in [1e-10,2e-12]:
        for h in [1e-3,5e-4]:
            def rhs(t,y):
                p=t*t*float(N(t))*np.sqrt(float(f(t)));v=-16*np.pi*gr.G*float(tr(t))*R*R/gr.C**4*float(N(t))/np.sqrt(float(f(t)))*t*t
                s=16*t*t*(1-t)**2;minus,jm,zero,j0,plus,jp,_,_=y
                return [jm/p,(1-h*s)*v*minus,j0/p,v*zero,jp/p,(1+h*s)*v*plus,
                    v*s*minus*plus,v*s*zero*zero]
            start=1e-6;aa=-16*np.pi*gr.G*float(tr(0))*R*R/gr.C**4
            y=[1+aa*start**2/6,float(N(0))*aa*start**3/3]*3+[0,0]
            sol=solve_ivp(rhs,(start,1),y,method='DOP853',rtol=tol,
                atol=[tol*.01,tol*1e-6]*3+[tol*1e-8]*2,max_step=.001)
            assert sol.success;a,jm,b,j0,c,jp,I,J=sol.y[:,-1]
            fd=-R*I/((a+H*jm)*(c+H*jp));derivative=-R*J/(b+H*j0)**2
            error=abs(fd-derivative);assert error<bound[str(h)]
            rows.append(dict(tol=tol,h=h,central_difference_via_exact_identity_m=fd,derivative_m=derivative,
                absolute_difference_m=error,analytic_truncation_bound_m=bound[str(h)]))
    save('finite-variation.json',dict(classification='Counterexample candidate',rows=rows,passed=True,
        identity='(chi(+h)-chi(-h))/(2h) = -integral phi(+h)*phi(-h)*v1 dr',
        advantage='Computes the exact nonlinear finite difference from a Wronskian integral, without subtracting independently integrated nearly equal monopoles.',
        independent_original_subtracted_results='scalar-variation.json',
        full_numerical_error_certificate=False,
        limit='Analytic truncation bound is rigorous under the envelope; floating ODE errors still lack an interval certificate.'))
    print('FINITE VARIATION',rows,flush=True)


def interval_plan():
    save('interval-plan.json',dict(classification='Conjectural',
        target='Outward interval enclosure of the fixed-GR-PCHIP potential-shape derivative, without treating floating ODE agreement as a proof.',
        coefficients='Treat saved binary64 PCHIP coefficients and radial knots as exact numbers defining this interpolant. '
        'Interval Horner evaluation encloses each polynomial on each subcell; no unproved continuum EOS/TOV error is omitted from the scope statement.',
        subdivisions=[1,4],interval_decimal_precision=40,
        method='Bound the static Green operator sup norm using positive envelopes of its absolute kernel. '
        'The convergent Neumann series gives a global interval for phi. Integrate phi^2 times the exact shape source with outward intervals.',
        gates='Operator norm <1/2; both enclosures contain independently integrated derivative; '
        '4-subcell enclosure narrower than 1-subcell enclosure; analytic constant-density sphere is enclosed.',
        excludes='Physical profile, EOS, atomic data, moving-star and observation derivatives. Exact interval enclosure is conditional on the frozen interpolant and specified constants.'))


def interval_derivative():
    from mpmath import iv
    from fractions import Fraction as F
    iv.dps=40
    def rational(endpoint):
        sign,man,exponent,_=endpoint._mpi_[0]
        return F((-1 if sign else 1)*man)*F(2)**exponent
    def enclosure(cells,H,R):
        accumulated=iv.mpf(0);central=iv.mpf(0)
        for left,right,a,b,weight in cells:
            source=abs(b);volume=(right**3-left**3)/3
            if left==0: central+=a*source*right**2/6
            else:
                inv=1/left-1/right
                central+=a*(accumulated*inv+source/3*((right**2-left**2)/2-left**3*inv))
            accumulated+=source*volume
        norm=central+H*accumulated;k=norm.b
        assert rational(k)<F(1,2)
        phi=1+iv.mpf([-1,1])*k/(1-k);derivative=iv.mpf(0)
        for _,_,_,b,weight in cells: derivative+=b*weight*phi**2
        derivative*=R
        return norm,derivative
    # Analytic uniform sphere with a fractional density perturbation.
    kk=iv.mpf('0.04')**2
    norm,control=enclosure([(iv.mpf(0),iv.mpf(1),iv.mpf(1),kk,iv.mpf(1)/3)],iv.mpf(1),iv.mpf(1))
    k=iv.mpf('0.04');analytic=(1/iv.cos(k)**2-iv.tan(k)/k)/2
    assert rational(control.a)<=rational(analytic.a)<=rational(analytic.b)<=rational(control.b)
    raw=np.load(ROOT/'outputs/remaining-closure26/scalar-background.npz');r=raw['r'];R=float(r[-1]);M=float(raw['mass'][-1]);x=r/R
    splines=[PchipInterpolator(x,np.exp(raw['nu'])),PchipInterpolator(x,1-2*np.divide(raw['mass'],r,out=np.zeros_like(r),where=r>0)),PchipInterpolator(x,raw['trace'])]
    coefficients=OUT/'scalar-polynomials.npz'
    if not coefficients.exists(): np.savez_compressed(coefficients,x=x,R=R,M=M,N=splines[0].c,f=splines[1].c,trace=splines[2].c)
    frozen=dict(np.load(coefficients))
    assert np.array_equal(x,frozen['x']) and R==frozen['R'] and M==frozen['M']
    assert all(np.array_equal(spline.c,frozen[key]) for spline,key in zip(splines,['N','f','trace']))
    mu=iv.mpf(M)/iv.mpf(R);H=-iv.ln(1-2*mu)/(2*mu)
    factor=16*iv.pi*iv.mpf('6.67428e-11')*iv.mpf(R)**2/iv.mpf(299792458)**4
    rows=[];reference=json.loads((OUT/'scalar-variation.json').read_text())['rows'][-1]['variational_derivative_m']
    for subdivisions in [1,4]:
        cells=[]
        for i in range(len(x)-1):
            edges=np.linspace(x[i],x[i+1],subdivisions+1)
            for lo,hi in zip(edges[:-1],edges[1:]):
                left=iv.mpf(float(lo));right=iv.mpf(float(hi));u=iv.mpf([float(lo),float(hi)])-iv.mpf(float(x[i]));values=[]
                for key in ['N','f','trace']:
                    value=iv.mpf(float(frozen[key][0,i]))
                    for coefficient in frozen[key][1:,i]: value=value*u+iv.mpf(float(coefficient))
                    values.append(value)
                n,f,tr=values;assert rational(n.a)>0 and rational(f.a)>0
                a=1/(n*iv.sqrt(f));b=factor*tr*n/iv.sqrt(f)
                primitive=lambda t:16*(t**5/5-t**6/3+t**7/7)
                weight=primitive(right)-primitive(left)
                assert rational(weight.a)>0
                cells.append((left,right,a,b,weight))
        norm,derivative=enclosure(cells,H,iv.mpf(R));lo,hi=rational(derivative.a),rational(derivative.b)
        assert lo<F(reference)<hi
        row=dict(subdivisions=subdivisions,cells=len(cells),operator_norm_upper_exact=str(rational(norm.b)),
            operator_norm_upper=float(rational(norm.b)),derivative_lower_exact_m=str(lo),derivative_upper_exact_m=str(hi),
            derivative_lower_m=float(lo),derivative_upper_m=float(hi),width_m=float(hi-lo),
            independent_ODE_value_inside=True)
        rows.append(row);print('INTERVAL',subdivisions,len(cells),row['operator_norm_upper'],float(lo),float(hi),flush=True)
    assert rows[1]['width_m']<rows[0]['width_m']
    save('interval-derivative.json',dict(classification='Proven',rows=rows,passed=True,
        coefficient_sha256=sha(coefficients),
        interval_library='mpmath.iv with 40 decimal digits; exact dyadic interval endpoints stored as rational numbers.',
        control_uniform_sphere_enclosed=True,
        control_lower_exact=str(rational(control.a)),control_upper_exact=str(rational(control.b)),
        target='Derivative of the exact binary64-coefficient PCHIP mathematical model under v(epsilon)=v0*(1+epsilon*16x^2(1-x)^2). '
        'G=6.67428e-11 and c=299792458 are treated as the specified exact constants.',
        proof='K has kernel G(r,s)*(-v(s)); its absolute integral is largest at r=0. '
        '||K||<=k<1/2 gives |phi-1|<=k/(1-k). Interval polynomial evaluation, exact monomial integrals, '
        'and interval arithmetic enclose integral (-v1)*phi^2, including the analytic Schwarzschild exterior Green factor.',
        physical_EOS_TOV_uncertainty_included=False,full_stellar_or_observational_derivative_certified=False))


def metric_variation():
    raw=np.load(ROOT/'outputs/remaining-closure26/scalar-background.npz');r=raw['r'];R=r[-1];M=raw['mass'][-1];x=r/R
    N=PchipInterpolator(x,np.exp(raw['nu']));f=PchipInterpolator(x,1-2*np.divide(raw['mass'],r,out=np.zeros_like(r),where=r>0))
    tr=PchipInterpolator(x,raw['trace']);mu=M/R;H=-np.log1p(-2*mu)/(2*mu)
    rows=[]
    for ps,vs in [(1,0),(1,1),(0,1)]:
        def rhs(t,y):
            p=t*t*float(N(t))*np.sqrt(float(f(t)));v=-16*np.pi*gr.G*float(tr(t))*R*R/gr.C**4*float(N(t))/np.sqrt(float(f(t)))*t*t
            shape=16*t*t*(1-t)**2;phi,j,z,jz,I=y
            return [j/p,v*phi,jz/p-ps*shape*j/p,v*z+vs*shape*v*phi,shape*(vs*v*phi**2+ps*j*j/p)]
        start=1e-6;aa=-16*np.pi*gr.G*float(tr(0))*R*R/gr.C**4
        sol=solve_ivp(rhs,(start,1),[1+aa*start**2/6,float(N(0))*aa*start**3/3,0,0,0],
            method='DOP853',rtol=2e-12,atol=[2e-14,2e-18,2e-18,2e-20,2e-22],max_step=.001)
        assert sol.success;a,j,z,jz,I=sol.y[:,-1];norm=a+H*j
        direct=-R*(jz*norm-j*(z+H*jz))/norm**2;adjoint=-R*I/norm**2
        err=abs(direct-adjoint)/max(abs(adjoint),1e-30);assert err<1e-8
        rows.append(dict(delta_p_shape=ps,delta_v_shape=vs,direct_variational_m=direct,integral_variation_m=adjoint,relative_error=err))
    import sympy as s
    p,pp,phi,phip,phipp,z,zp,zpp,v,dv,h,hp=s.symbols('p pp phi phip phipp z zp zpp v dv h hp')
    Wprime=pp*(phi*zp-z*phip)+p*(phi*zpp-z*phipp)
    fluxprime=phip*h*phip+phi*hp*phip+phi*h*phipp
    expr=Wprime+fluxprime-dv*phi**2-h*phip**2
    expr=expr.subs(zpp,(v*z+dv*phi-hp*phip-h*phipp-pp*zp)/p).subs(phipp,(v*phi-pp*phip)/p)
    assert s.simplify(expr)==0
    save('metric-variation.json',dict(classification='Counterexample candidate',rows=rows,passed=True,
        formula='delta chi = -integral [phi^2 delta v + (phi prime)^2 delta p] dr',
        definitions='p=r^2 N sqrt(f), v=4pi G beta r^2 N (e-3P)/(c^4 sqrt(f)); phi(infinity)=1.',
        charge_convention='phi_physical=phi_infinity + phi_infinity*chi/r+... = phi_infinity - alpha_A*M_geom/r+...; alpha_A=-phi_infinity*chi/M_geom.',
        normalized_variation='At fixed phi_infinity: delta alpha_A=-phi_infinity*(delta chi/M-chi*delta M/M^2).',
        boundary='Fixed areal radial coordinate with exterior included. Smooth perturbations, or distributional moving interfaces with their boundary terms retained. '
        'Asymptotic delta p=O(r), delta phi=O(1/r) makes the boundary correction vanish. '
        'The numerical shapes have compact support; no Einstein constraint satisfaction is asserted.',
        theorem_symbolically_checked=True,physical_metric_variation_derived=False))
    print('PASS metric variation and mass-normalized convention',rows,flush=True)


def untraced():
    data=np.load(OUT/'HELM-coupled-cell-4-end.npz');X=data['X'];lnT=float(data['lnT']);label='HELM-untraced-final'
    base=dict(np.load(OLD/'PP-steady-2-input.npz'));base['X']=base['X'].copy();base['lnT']=base['lnT'].copy()
    base['X'][ZONE]=X;base['lnT'][ZONE]=lnT;context()
    native.setup(label,base,extra=' use_eosDT_HELMEOS=.true.',species=species(),network=(OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    fresh.run(label);fresh.collect(label)
    _,a=mesa(OUT/(label+'-profile.data.gz'));row=json.loads((OUT/'HELM-coupled-continuation.json').read_text())['rows'][-1]
    _,b=mesa(OUT/(row['records'][-1]['selected']+'-profile.data.gz'))
    errors={k:float(max(abs(a[k]-b[k])/np.maximum(1e-100,abs(b[k])))) for k in ['rho','logT','energy','pressure','eps_nuc','eps_nuc_neu_total','non_nuc_neu']}
    assert max(errors.values())<1e-12
    save('untraced-control.json',dict(classification='Counterexample candidate',passed=True,errors=errors))


def summary():
    q,w=rest();rows=json.loads((OUT/'HELM-coupled-continuation.json').read_text())['rows'];assert all(r['completed'] and r['strict_energy_passed'] for r in rows)
    assert all(r['refinement']['passed'] for r in rows[1:])
    original=json.loads((OUT/'continuation.json').read_text());split=json.loads((OUT/'HELM-continuation.json').read_text())
    assert not original['rows'][-1]['completed'] and not any(r['strict_energy_passed'] for r in original['rows'])
    assert not any(r['refinement']['passed'] for r in split['rows'][1:])
    largest=0.;thermal=0.
    for path in sorted(OUT.glob('HELM-coupled-*-native.npz')):
        d={k:v[ZONE] for k,v in np.load(path).items()};_,p=mesa(OUT/path.name.replace('-native.npz','-profile.data.gz'))
        largest=max(largest,abs(float(q@d['dxdt']+d['heat']+d['neutrino']))/max(1,abs(float(d['heat']))))
        thermal=max(thermal,abs(float(p['non_nuc_neu'][ZONE]))/max(1,abs(float(d['heat']))))
    last=rows[-1];c2=(gr.C*100)**2
    save('closure-summary.json',dict(classification='Counterexample candidate',
        all_coupled_energy_and_refinement_gates_passed=True,
        worst_coupled_step_energy_defect_over_initial_release=max(t['energy_score'] for r in rows for t in r['records']),
        final_temperature_log_change=last['delta_lnT'],final_max_abundance_change=last['max_delta_X'],
        final_rest_mass_fraction_change=last['Q_rest_change_erg_g']/c2,
        final_thermal_mass_fraction_change=last['internal_energy_change_erg_g']/c2,
        final_cell_total_mass_fraction_change=(last['Q_rest_change_erg_g']+last['internal_energy_change_erg_g'])/c2,
        final_escaping_neutrino_energy_over_rest_mass=last['neutrino_loss_erg_g']/c2,
        max_native_Q_source_identity_relative=largest,max_thermal_neutrino_fraction=thermal,
        original_default_failure_retained=True,split_HELM_time_failure_retained=True,
        physical_common_EOS_certified=False,full_GR_evolution=False,scalar_or_force_observable_derived=False,
        status='Conditional hot isochoric reaction/thermal continuation passed. Strict interval derivative of frozen scalar interpolant obtained. '
        'Physical stellar EOS, metric/fluid evolution and scalar drive/readout remain separate unresolved requirements.'))
    print('SUMMARY native Q identity residual',largest,'thermal nu fraction',thermal,flush=True)


def maintain():
    bodies=[
        '분류: Counterexample candidate. Request28은 26종 전체 조성·온도·누적 중성미자 에너지를 결합한 고정 부피 대조를 수행했다. 기본 EOS의 에너지 역산 실패와 온도를 분리한 HELM의 시간 정밀화 실패를 보존했다. 별도 고온 HELM 결합 대조는 같은 엄격 에너지·조성·온도 기준과 추가 중성미자 적분 기준을 통과했다. 물리적 공통 EOS나 GR 항성 진화의 인증은 아니다.',
        '분류: Proven. 정적 감수율의 변분은 δχ=−∫[φ²δv+(φ′)²δp]dr이며, 고정 φ∞에서 α_A=−φ∞χ/M의 변분에는 질량 정규화 항도 들어간다. 보편 결합의 영 결합에너지 한계에서 Q=αM이면 내부 에너지 재분배로 δ(Q/M)이 생기지 않는다. 고정 GR 보간 모형의 지정 퍼텐셜 모양 미분은 [251.36919,251.70105] m로 구간 보증했다.',
        '분류: Proven. 내부 핵 전환은 닫힌 계의 총에너지에 독립적인 원천으로 더해지지 않는다. 반응열·정지에너지·조성별 내부에너지와 경계 복사 손실을 함께 계수해야 한다. 이전 PP 준정상 부분계의 가열은 벌크 연료를 필요로 하므로 네 중간 핵종의 상태를 닫힌 질량 저장고로 읽을 수 없다. 영 scalar 가지와 유한 carrier 정적 보간 경계는 유지한다.',
        '분류: Counterexample candidate. 외부 열욕을 제거한 별도 HELM 고정 부피 계산에서 26종 반응과 온도를 함께 진행했다. 2→4분할 조성 오차/허용량은 0.02855, 온도 ln 차이는 2.35e-11, 중성미자 손실 상대 차이는 1.25e-5다. 이는 지정 국소 근사의 시간 대조이며 실제 항성의 열·유체 모드나 scalar 힘의 비선형 전달함수는 아니다.',
        '분류: Counterexample candidate. 기본 EOS의 온도 역산은 엄격 에너지 기준에 실패했고 4분할 첫 단계에서 별도 진단 한계도 넘어 중단했다. HELM 분리 적분도 조성 시간 기준에 실패했다. 이 실패를 보존하고 조성·온도·손실의 결합 적분을 별도로 통과했다.\n\n분류: Proven. 고정 GR 보간 모형의 지정 scalar 미분은 구간 연산으로 보증했으므로 이 좁은 항목은 해소되었다.\n\n분류: Conjectural. 물리적 공통 EOS와 에너지 영점 교정, 보존된 26종 GR 초기화·동적 진화, 실제 scalar 구동 및 질량 정규화 전하·관측 추론은 남는다.'
    ]
    rev=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in zip(DOCS,bodies):
        rel='docs/'+stem+'.md';binding=json.loads((OUT/'historical-note-bindings.json').read_text())[rel]
        assert sha(ROOT/rel)==binding['sha256'],rel
        with (ROOT/rel).open('ab') as f: f.write(('\n\n## Request 28 반응 에너지 보존과 scalar 미분 구간\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST28_CONSERVATIVE_CELL_KO.md).\n').encode())
        rev['sha256'][rel]=sha(ROOT/rel)
    rev['request28_supporting_note_update']=dict(evidence_manifest='outputs/conservative-cell28/manifest.json',
        historical_notes='outputs/conservative-cell28/historical-note-bindings.json',
        status='고온 HELM 26종 반응·열 결합 대조 및 고정 GR scalar 미분 구간 보증; 전체 물리·관측 폐쇄 미완료',
        artifact_status='Request12 원고 PDF/ZIP은 역사 산출물로 보존한다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(rev,ensure_ascii=False,indent=2)+'\n')


def seal():
    import mpmath, scipy, sympy
    context();native.audit();summary()
    deps=Path(mpmath.__file__).parent
    save('provenance.json',dict(classification='Proven',checkpoint='bf2f600',previous_manifest_sha256=sha(OLD/'manifest.json'),
        executable_sha256=sha(fresh.BINARY),runtime_data_binding='outputs/fresh-microphysics25/runtime-data-bindings.json',
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,sympy=sympy.__version__,mpmath=mpmath.__version__,
        interval_source_root=str(deps),interval_sources={p.relative_to(deps).as_posix():sha(p) for p in sorted(deps.rglob('*.py'))},
        historical_policies='Preserve all previous failures; no runtime rebuild/install and no new empirical observations.',
        literature=[dict(title='MESA energy equations',url='https://docs.mesastar.org/en/latest/reference/controls.html'),
            dict(title='Damour and Esposito-Farese (1996)',url='https://arxiv.org/abs/gr-qc/9602056'),
            dict(title='mpmath interval arithmetic',url='https://mpmath.org/doc/current/contexts.html#arbitrary-precision-interval-arithmetic-iv')]))
    save('gates.json',dict(classification='Proven',conditional_HELM_cell_energy_and_time_controls=True,
        fixed_GR_interpolant_scalar_derivative_interval=True,metric_variation_and_mass_normalization_derived=True,
        default_EOS_strict_energy_gate=False,split_HELM_composition_time_gate=False,
        physical_common_EOS_certified=False,full_GR_thermal_fluid_metric_evolution=False,
        scalar_drive_charge_map=False,complete_nonlinear_observation=False,final_submission_package_updated=False,
        classification_detail='Theorem progress: interval derivative, metric variation and normalization no-go. '
        'Loophole progress: native 26-species isochoric reaction/temperature/energy continuation in the specified HELM approximation.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/conservative_cell.py',ROOT/'notes/REQUEST28_CONSERVATIVE_CELL_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27','conservative-cell28']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/conservative-cell28/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,v in before.items():
                        if k!='sha256': assert after[k]==v,k
                    for k,v in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==v,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name);count+=1
    prov=json.loads((OUT/'provenance.json').read_text());assert sha(fresh.BINARY)==prov['executable_sha256']
    for path,digest in prov['interval_sources'].items(): assert sha(Path(prov['interval_source_root'])/path)==digest,path
    data=json.loads((ROOT/'outputs/fresh-microphysics25/runtime-data-bindings.json').read_text())
    for path,digest in data['sha256'].items(): assert sha(Path(data['root'])/path)==digest,path
    assert sha(OUT/'continuation.json')==json.loads((OUT/'HELM-plan.json').read_text())['default_failure_sha256']
    assert sha(OUT/'HELM-continuation.json')==json.loads((OUT/'coupled-plan.json').read_text())['split_failure_sha256']
    assert json.loads((OUT/'interval-derivative.json').read_text())['coefficient_sha256']==sha(OUT/'scalar-polynomials.npz')
    gates=json.loads((OUT/'gates.json').read_text())
    for key in ['default_EOS_strict_energy_gate','split_HELM_composition_time_gate','physical_common_EOS_certified',
        'full_GR_thermal_fluid_metric_evolution','scalar_drive_charge_map','complete_nonlinear_observation','final_submission_package_updated']: assert not gates[key]
    print('PASS',count,'artifact/history SHA;',len(data['sha256']),'native data SHA;',len(prov['interval_sources']),'interval source SHA',flush=True)


def recheck():
    import tempfile,contextlib,io
    global OUT,CACHE
    original,cache=OUT,CACHE
    labels=['cell-4-0-2','HELM-cell-4-3-1','HELM-coupled-cell-4-3-1','HELM-gradient-c12-5e-07-1']
    with tempfile.TemporaryDirectory(prefix='conservative-cell28-') as temp:
        OUT=Path(temp)/'out';CACHE=Path(temp)/'runs';OUT.mkdir();CACHE.mkdir();context()
        try:
            for p in original.glob('*.json'): shutil.copy2(p,OUT/p.name)
            shutil.copy2(original/'scalar-polynomials.npz',OUT/'scalar-polynomials.npz')
            # Cached replay checks all coupled steps; four selected states are re-executed below.
            for p in original.glob('HELM-coupled-*'):
                if p.is_file(): shutil.copy2(p,OUT/p.name)
            for suffix in ['-input.npz','-native.npz','-profile.data.gz']:
                shutil.copy2(original/('HELM-initial'+suffix),OUT/('HELM-initial'+suffix))
            worst=0.
            for label in labels:
                folder=CACHE/label;folder.mkdir();shutil.copy2(fresh.BINARY,folder/'binary')
                shutil.copytree(original/'inputs'/label,OUT/'inputs'/label)
                shutil.copy2(original/(label+'-input.npz'),OUT/(label+'-input.npz'))
                for p in (original/'inputs'/label).iterdir():
                    if p.name=='input.mod.gz':
                        with fresh.gzip.open(p,'rb') as src,(folder/'input.mod').open('wb') as dst: shutil.copyfileobj(src,dst)
                    else: shutil.copy2(p,folder/p.name)
                with contextlib.redirect_stdout(io.StringIO()): native.trace(label)
                a=np.load(OUT/(label+'-native.npz'));b=np.load(original/(label+'-native.npz'))
                for key in a.files:
                    err=float(np.max(abs(a[key]-b[key])/np.maximum(1e-100,abs(b[key]))));worst=max(worst,err);assert err<1e-10,(label,key)
                _,a=mesa(OUT/(label+'-profile.data.gz'));_,b=mesa(original/(label+'-profile.data.gz'))
                for key in ['energy','pressure','eps_nuc','eps_nuc_neu_total','non_nuc_neu']:
                    assert max(abs(a[key]-b[key])/np.maximum(1e-100,abs(b[key])))<1e-12,(label,key)
                print('RECHECK actual native state',label,flush=True)
            with contextlib.redirect_stdout(io.StringIO()):
                coupled();symbolic();scalar();bounds();finite_variation();interval_derivative();metric_variation()
            a=json.loads((OUT/'HELM-coupled-continuation.json').read_text());b=json.loads((original/'HELM-coupled-continuation.json').read_text())
            assert a==b
            assert json.loads((OUT/'interval-derivative.json').read_text())==json.loads((original/'interval-derivative.json').read_text())
            print('PASS 4 actual native reruns, all coupled-step replay, symbolic/variational/interval recomputation; max native difference',worst,flush=True)
        finally: OUT,CACHE=original,cache;context()


def git_blobs():
    import hashlib,subprocess
    expected=dict(json.loads((OUT/'manifest.json').read_text())['sha256'])
    expected['outputs/conservative-cell28/manifest.json']=sha(OUT/'manifest.json')
    with subprocess.Popen(['git','cat-file','--batch'],cwd=ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for rel,digest in expected.items():
            proc.stdin.write(('HEAD:'+rel+'\n').encode());proc.stdin.flush()
            header=proc.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(rel,header)
            remaining=int(header[2]);actual=hashlib.sha256()
            while remaining:
                chunk=proc.stdout.read(min(remaining,1048576));assert chunk
                actual.update(chunk);remaining-=len(chunk)
            assert proc.stdout.read(1)==b'\n' and actual.hexdigest()==digest,rel
        proc.stdin.close();assert proc.wait(timeout=10)==0
    result=dict(classification='Proven',passed=True,raw_blobs=len(expected),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
    (CACHE/'git-object-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS Git raw blobs',result,flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

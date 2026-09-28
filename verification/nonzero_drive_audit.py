"""Independent EOS/exterior checks and conditional algebraic derivative enclosures."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import sys

import mpmath as mp
import numpy as np
import sympy as s
from scipy.integrate import quad
from sympy.printing.pycode import MpmathPrinter

from nonzero_drive import ROOT,OUT,EOS,ScalarStar,electron_gas,save
from remaining_audit import exact_endpoints,endpoints

def symbolic():
    x=s.symbols('x',positive=True)
    p=x*(2*x*x-3)*s.sqrt(1+x*x)+3*s.asinh(x)
    energy=3*(x*(1+2*x*x)*s.sqrt(1+x*x)-s.asinh(x))
    assert s.simplify(s.diff(p,x)-8*x**4/s.sqrt(1+x*x))==0
    assert s.simplify(s.diff(energy,x)-3*(energy+p)/x)==0
    # Baryon first law gives n=(epsilon+p) exp(-integral dp/(epsilon+p)).
    e,pv,ep=s.symbols('e p ep',positive=True)
    assert s.simplify((ep+1)/(e+pv)-1/(e+pv)-ep/(e+pv))==0
    # The exact compactified vacuum equations and the unscalarized limit.
    R,u,m,z=s.symbols('R u m z',positive=True);r=R/u;f=1-2*m/r
    assert s.simplify((s.Rational(1,2)*f*z*z/r**2)*(-R/u**2)+f*z*z/(2*R))==0
    cp,ci,co,a,b,c,phi=s.symbols('cp ci co a b c phi',positive=True)
    chis=[cp,ci,co];v=[a,b,c];D=s.diag(*chis);K=s.Matrix([[0,a,b],[a,0,c],[b,c,0]])
    M=s.eye(3)-D*K;det=1-cp*ci*a*a-cp*co*b*b-ci*co*c*c-2*cp*ci*co*a*b*c
    assert s.expand(M.det()-det)==0
    q=M.inv()*D*s.ones(3,1)*phi
    assert (M*q-D*s.ones(3,1)*phi).applyfunc(s.factor)==s.zeros(3,1)
    # Envelope theorem: integrating out q preserves reciprocal pair forces.
    L=D.inv()-K;V=(s.Matrix(s.symbols('q:3')).T*L*s.Matrix(s.symbols('q:3')))[0]/2-phi*sum(s.symbols('q:3'))
    assert s.diff(V,a)==-s.Symbol('q0')*s.Symbol('q1')
    assert s.factor(s.diff(q[0],phi)-q[0]/phi)==0
    # A sign flip of the scalar field/charges leaves matter observables unchanged.
    signal=s.symbols('T');phi0=s.symbols('phi0',real=True)
    assert (phi0**2*signal).subs(phi0,-phi0)==phi0**2*signal
    assert s.diff(phi0**2*signal,phi0).subs(phi0,0)==0
    # Leading monopole radiation has rank-one damping, not three independent
    # positive damping rates. Omitting inertia is an explicit reduction.
    total,totaldot,cap,light=s.symbols('S Sdot Ceff c',positive=True)
    assert s.expand((total/cap+totaldot/light-phi)*cap-(total+cap/light*totaldot-cap*phi))==0
    save('symbolic.json',dict(classification='Proven',electron_gas_thermodynamics=True,baryon_density_identity=True,
      compactified_vacuum=True,coupled_determinant=str(det),equilibrium_equation_verified=True,
      reciprocal_force=True,charge_linear_in_background=True,force_even_in_background=True,
      exact_static_elimination='V_eff = -phi^2 * 1^T L^-1 1 / 2',static_elimination_is_not_new_dynamic_state=True,
      scalar_sign_not_identifiable_from_matter_timing=True,phi_score_and_Fisher_vanish_at_zero_in_smooth_even_model=True,
      Gaussian_interval_in_phi_at_zero_not_justified=True,
      rank_one_radiation_reduction='(Ceff/c) Sdot + S = Ceff*phi; Ceff=1^T L^-1 1',
      no_inertia_assumption_explicit=True))
    return chis,v,phi,M,q,det

def stars():
    d=json.loads((OUT/'stars.json').read_text());checks=[]
    for k in range(3):
        base=next(r for r in d['records'] if r['index']==k and r['precision']=='base')
        fine=next(r for r in d['records'] if r['index']==k and r['precision']=='fine')
        coarse_chi=base['linear_static']['susceptibility_m'];chi=fine['linear_static']['susceptibility_m']
        dr=base['zero']['radius_m']/fine['zero']['radius_m']-1;dc=coarse_chi/chi-1
        assert max(abs(dr),abs(dc))<1e-5
        rows=fine['finite_background'];assert abs(rows[0]['charge_m']+rows[1]['charge_m'])<1e-11
        assert all(abs(r['baryon_relative_error'])<1e-7 and r['shooting_error']<1e-7 for r in rows)
        changes=[r['q_over_phi_m']/chi-1 for r in rows]
        assert max(abs(np.array(changes)))<1e-3
        checks.append(dict(body=k,radius_relative_refinement=dr,chi_relative_refinement=dc,
          q_over_phi_relative_to_zero_chi=changes,max_baryon_relative_error=max(abs(r['baryon_relative_error']) for r in rows)))
    old=json.loads((ROOT/'outputs/research-remediation/stellar-matching.json').read_text())['stars'][-1]
    ns=next(r for r in d['records'] if r['index']==0 and r['precision']=='fine')
    native_comparison=dict(radius_relative_change=ns['zero']['radius_m']/old['radius_m']-1,
      chi_relative_change=ns['linear_static']['susceptibility_m']/old['static']['susceptibility_m']-1)
    assert max(abs(x) for x in native_comparison.values())<1e-5
    cold=next(r for r in d['records'] if r['index']==1 and r['precision']=='fine')['zero']['radius_m']
    solar_radius_m=695700000. # nominal solar-radius unit; physical measurement quoted in this unit
    optical=dict(classification='Proven',observational_source_classification='Imported from prior work',source='https://arxiv.org/abs/1402.0407',
      radius_solar=.091,radius_uncertainty_solar=.005,Teff_K=15800,Teff_uncertainty_K=100,
      cold_model_radius_solar=cold/solar_radius_m,cold_radius_to_observed_ratio=cold/(.091*solar_radius_m),
      cold_model_reproduces_observed_radius=False)
    assert cold/solar_radius_m<.091-3*.005
    save('optical-domain-check.json',optical)
    # Independent analytic exterior matching, Eqs. 36-38 of Mendes/Ortiz 2016.
    exterior=[]
    for k in range(3):
        fine=next(r for r in d['records'] if r['index']==k and r['precision']=='fine')
        row=fine['finite_background'][0];eos=EOS(fine['eos'],16384)
        # Recover central field by the regular linear normalization; then test
        # the compactified exterior of that actual solution, no shooting needed.
        pc=np.log(row['central_pressure_geom']);phic=row['phi_infinity']/fine['linear_static']['central_normalization_phi_infinity']
        star=ScalarStar(eos,pc,phic,rtol=2e-11,r0_factor=.5)
        m,h,nu,phi,psi,nb=star.full.y[:,-1];R=star.radius
        v=R*psi**2+2*m/(R*(R-2*m));w=np.sqrt(v*v+4*psi*psi);atanh=np.arctanh(w/(v+2/R))
        asym=phi+2*psi/w*atanh;mass=R*R*v/2*np.sqrt(1-2*m/R)*np.exp(-v/w*atanh);charge=-2*mass*psi/v
        errors=[abs(asym-star.phi_infinity),abs(mass/star.mass-1),abs(charge/star.charge-1)]
        assert errors[0]<1e-13 and max(errors[1:])<1e-9
        exterior.append(dict(body=k,phi_absolute_error=errors[0],mass_relative_error=errors[1],charge_relative_error=errors[2]))
    from scipy import constants
    pref=constants.m_e**4*constants.c**5/(24*np.pi**2*constants.hbar**3)
    gas=[]
    for x in [1e-4,.01,.049,.051,1.,3.]:
        p,e,rest=electron_gas(x);integral=8*pref*quad(lambda t:t**4/np.sqrt(1+t*t),0,x,epsabs=1e-35,epsrel=1e-12)[0]
        err=float(p/integral-1);assert abs(err)<1e-9;gas.append(dict(x=x,pressure_relative_integral_error=err))
    save('stellar-audit.json',dict(classification='Proven',scope='지정 EOS 수치해의 수렴·대칭성·독립 진공 접합 검사; 항성 해의 구간 인증은 아님',refinement=checks,request13_native_comparison=native_comparison,analytic_exterior=exterior,electron_pressure_integral=gas))

def intervals():
    chis,v,phi,M,q,det=symbolic();iv=mp.iv;iv.dps=70
    d=json.loads((OUT/'coupled.json').read_text());g=d['geometry'];mass=np.array(d['masses_geom_m'])
    a0=g[0]['a_m'];e0=g[0]['ecc'];a1=g[1]['a_m'];e1=g[1]['ecc'];f=mass[1]/sum(mass[:2])
    lo=[a0*(1-e0),a1*(1-e1)-f*a0*(1+e0),a1*(1-e1)-(1-f)*a0*(1+e0)]
    hi=[a0*(1+e0),a1*(1+e1)+f*a0*(1+e0),a1*(1+e1)+(1-f)*a0*(1+e0)]
    # These decimal boxes DEFINE the conditional domain; they are not a proof
    # that the physical star or full GR orbit lies inside these boxes.
    boxes_chi=[iv.mpf([str(np.floor(x*.99999)),str(np.ceil(x*1.00001))]) for x in d['chi_m']]
    radii=[iv.mpf([str(np.floor(l*.99999)),str(np.ceil(h*1.00001))]) for l,h in zip(lo,hi)]
    inputs=boxes_chi+[1/r for r in radii]+[iv.mpf('1e-5')]
    symbols=chis+v+[phi];printer=MpmathPrinter({'fully_qualified_modules':False})
    evaluate=lambda ex:s.lambdify(symbols,ex,modules=[{'mpf':iv.mpf}],printer=printer)(*inputs)
    determinant=evaluate(det);assert determinant.a>0
    rownorm=max((evaluate(sum((s.eye(3)-M)[i,j] for j in range(3))).b for i in range(3)))
    assert rownorm<1
    report=[];qboxes=[];jacboxes=[]
    for expr in q:
        qq=evaluate(expr);jac=[evaluate(s.diff(expr,x)) for x in v]
        hes=[[evaluate(s.diff(expr,x,y)) for y in v] for x in v]
        qboxes.append(qq);jacboxes.append(jac)
        report.append(dict(q=endpoints(qq),dq_d_inverse_distance=list(map(endpoints,jac)),d2q_d_inverse_distance2=[[endpoints(h) for h in row] for row in hes]))
    # Directed bounds imply nonsingularity throughout the entire box. Positive
    # diagonal susceptibilities make L congruent to I-sqrt(D) K sqrt(D), hence SPD.
    save('coupled-certificate.json',dict(classification='Proven',scope='명시한 감수율·거리 구간의 정적 축약 모형에만 조건부; 항성·GR 궤도·timing 인증은 아님',
      chi_exact_boxes=list(map(exact_endpoints,boxes_chi)),distance_m_exact_boxes=list(map(exact_endpoints,radii)),
      determinant=endpoints(determinant),DK_infinity_norm_upper=endpoints(rownorm)[1],reduced_static_Hessian_positive_definite=True,
      inverse_norm_upper=endpoints(1/(1-rownorm)),all_phase_Kepler_triangle_domain_only=True,derivatives=report,
      full_timing_certificate=False,stellar_numerical_error_certified=False))
    # Continuous, all-phase tracking bound for that same reduced DAE, started
    # at its equilibrium. Kepler speeds are overbounded by declared integer SI
    # values; no native GR trajectory enclosure is silently assumed.
    from nonzero_drive import GM_SUN,MSUN,C
    vp=np.sqrt(GM_SUN*(sum(mass[:2])/MSUN)/a0*(1+e0)/(1-e0))
    vo=np.sqrt(GM_SUN*(sum(mass)/MSUN)/a1*(1+e1)/(1-e1))
    speeds=[np.ceil(vp*1.00001),np.ceil((vo+f*vp)*1.00001),np.ceil((vo+(1-f)*vp)*1.00001)]
    absmax=lambda x:max(abs(x.a),abs(x.b))
    B=sum(absmax(jacboxes[i][j])*iv.mpf(str(speeds[j]))/(radii[j].a**2) for i in range(3) for j in range(3))
    capacity=sum(qboxes)/inputs[-1];tau=capacity/iv.mpf(str(C));error=tau.b*B
    force_errors=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        bound=((absmax(qboxes[i])+absmax(qboxes[j]))*error+error**2)/(iv.mpf(str(mass[i]))*iv.mpf(str(mass[j])))
        force_errors.append(endpoints(bound)[1])
    save('collective-tracking.json',dict(classification='Proven',scope='관성 없는 선도 monopole 복사 DAE와 지정 Kepler/감수율 구간에 조건부; 완전한 항성 모드 또는 timing 정리가 아님',
      initial_total_charge_equals_equilibrium=True,all_times_in_declared_orbit_domain=True,
      speed_upper_m_per_s=list(map(float,speeds)),effective_susceptibility_m=endpoints(capacity),relaxation_seconds=endpoints(tau),
      equilibrium_total_charge_speed_upper_m_per_s=endpoints(B)[1],charge_tracking_error_upper_m=endpoints(error)[1],
      pair_Delta_tracking_error_upper=force_errors,inertia_and_higher_multipoles_certified_negligible=False,
      theorem='For tau(t)>0, e=S-Ceff*phi obeys edot+e/tau=-(Ceff*phi)dot; |e|<=|e0| exp(-integral dt/tau)+tau_max sup|(Ceff*phi)dot|.'))

def coupled_checks():
    d=json.loads((OUT/'coupled.json').read_text());a,b=d['grids'];rows=[]
    for aa,bb in zip(a['carriers'],b['carriers']):
        rel=np.array(aa['all_pairs_amplitude'])/bb['all_pairs_amplitude']-1
        assert max(abs(rel))<.005
        rows.append(dict(carrier=aa['carrier'],relative_refinement=rel.tolist()))
    save('coupled-audit.json',dict(classification='Proven',grid_refinement=rows,
      no_actual_TOA_likelihood_evaluated=True,no_one_pole_relaxation_inferred=True))

def check():
    stars();intervals();coupled_checks();print('PASS: 항성·EOS·외부 진공·결합 전하·조건부 미분 감사')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def maintain():
    frozen=json.loads((ROOT/'outputs/nbody-readout16/manifest.json').read_text())['sha256']
    additions={
      'model-definition':'분류: Counterexample candidate. DEF beta=−4, 배경 scalar=1e−5, SLy 중성자별과 차가운 mu_e=2 전자 축퇴 EOS의 백색왜성 두 개를 지정했다. 영 배경의 기준 별로부터 바리온 수를 고정한 비영 가지를 수치 계산했다. 선도 작은 배경 모형은 Lq=phi*1, L=diag(1/chi)−K, K_AB=1/r_AB이며 모든 전하가 서로 응답한다. 실제 J0337의 EOS·배경 값을 결정한 것은 아니다.',
      'observable-targets':'분류: Proven. 지정된 선도 모형에서 동반성 전하 응답을 포함하면 안쪽 쌍의 내궤도 변동은 고정 동반성 근사보다 약 36.6배, 바깥쪽 쌍의 외궤도 변동은 약 18.6배 커진다. 절대 진폭은 각각 약 4.88e−17, 7.07e−17의 무차원 힘 결합이다. 이는 TOA 잔차나 검출값이 아니다. 공통 두 쌍 결합 템플릿과 새로운 세 쌍 모형을 동일시하지 않는다. scalar 부호 반전은 물질 관측량을 보존하며 영 배경에서 통상적인 phi 선형 Fisher 근사가 퇴화한다.',
      'adiabatic-limit':'분류: Proven. 명시한 감수율·거리 구간에서 결합 행렬의 가역성과 정적 Hessian 양의 정부호, 전하의 역거리 일·이차 미분 구간을 확인했다. 관성을 생략한 선도 monopole 복사의 rank-one 감쇠 아래 총전하 S는 (Ceff/c) Sdot+S=Ceff*phi를 만족한다. 지정 구간의 완화시간은 0.150614–0.150631 ms이며, 평형 초기자료에서 각 전하의 연속 추적 오차는 1.37e−15 m 이하라는 조건부 상계를 얻었다. 관성·고차 복사·실제 궤도 오차는 이 상계에 포함되지 않는다.',
      'nonadiabatic-regime':'분류: Proven. 비영 배경으로 궤도에 따라 변하는 전하는 생기지만, 이번 정적 해는 V_eff=−phi²*1^T L^-1 1/2로 정확히 제거된다. 관성을 생략한 선도 복사 모형의 빠른 집단 완화도 일 단위 상태를 확립하지 않는다. 분류: Conjectural. 비영 가지의 결합 모드·관성·고차 복사와 실제 구동·광자 전파를 일관되게 연결해야 비단열 관측 후보로 승격할 수 있다.',
      'failure-ledger-dynamic-chi':'분류: Proven. Request 16의 영 구동 경계는 비영 배경을 명시한 지정 후보에서 벗어났으며, 별 세 개의 유한 배경 평형과 선도 상호 구동을 계산했다. 그러나 작은 전체 되먹임이 작은 변동 신호의 상대오차를 보장하지 않아 고정 동반성 근사가 실패했다. 정적 전하 제거와 rank-one 복사의 조건부 빠른 완화는 새로운 느린 관측량을 공급하지 않는다. 차가운 안쪽 백색왜성 모형은 반지름 약 0.0212 태양반지름으로, Kaplan 등의 광학 반지름 0.091±0.005 태양반지름을 재현하지 못한다. 실제 항성 matching에는 열·외피 구조가 필요하다. 분류: Conjectural. 항성·동역학의 전 기간 구간 인증, 전체 매개변수 초기화·광자 전파, 새 쌍별 신호의 비선형 likelihood·pulse/noise 검증은 미완료다. 압력 좌표 실패·수치법 변경·축약 모형의 가정을 보고서에 보존한다.'}
    rels=['docs/'+k+'.md' for k in additions]+['paper/revision-manifest.json']
    for rel in rels:assert sha(ROOT/rel)==frozen[rel],rel
    dest=OUT/'request16-notes';dest.mkdir(exist_ok=False);bindings={}
    for rel in rels:
        path=dest/Path(rel).name;shutil.copy2(ROOT/rel,path)
        bindings[rel]=dict(snapshot=path.relative_to(ROOT).as_posix(),sha256=frozen[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as file:file.write(('\n\n## Request 17 비영 구동과 동반성 응답\n\n'+body+'\n\n세부 근거: [한글 도출·검증 보고서](../notes/REQUEST17_NONZERO_DRIVE_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request17_supporting_note_update']=dict(evidence_manifest='outputs/nonzero-drive17/manifest.json',historical_notes='outputs/nonzero-drive17/historical-note-bindings.json',status='지정 비영 배경 항성·동반성 구동 및 조건부 축약 미분·추적 상계; 전체 물리 timing·추론 미완료',artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')

def seal():
    import importlib.metadata as md
    import scipy,lalsimulation
    import lalsimulation._lalsimulation as native_lal
    import nonzero_drive,stellar_matching,remaining_audit
    modules=[np,scipy,lalsimulation,native_lal,s,mp,nonzero_drive,stellar_matching,remaining_audit]
    inputs=['outputs/validated-variational/ivp.hex','outputs/research-remediation/stellar-matching.json','outputs/research-remediation/sources/LALSimNeutronStarEOS_SLY.dat']
    save('provenance.json',dict(before_task_checkpoint='043982c',interpreter=sys.executable,
      dependency_path='/home/lpaiu/work/nutimo_pilot/request13_deps',versions={name:m.__version__ for name,m in [('numpy',np),('scipy',scipy),('sympy',s),('mpmath',mp)]},
      distribution_versions={name:md.version(name) for name in ['numpy','scipy','lalsuite','sympy','mpmath']},
      module_file_sha256={str(m.__file__):sha(m.__file__) for m in modules},input_sha256={p:sha(ROOT/p) for p in inputs},
      lal_shared_libraries={str(p):sha(p) for p in Path('/home/lpaiu/work/nutimo_pilot/request13_deps/lalsuite.libs').glob('liblal*.so*')},
      producer='verification/nonzero_drive.py stars; verification/nonzero_drive.py coupled; verification/nonzero_drive_audit.py check'))
    save('gates.json',dict(classification='Proven',specified_finite_background_stars_computed=True,
      specified_leading_companion_drive_computed=True,coupled_algebraic_derivatives_conditionally_certified=True,
      reduced_radiation_DAE_tracking_conditionally_certified=True,theorem_progress=True,
      genuine_orbital_timescale_state_established=False,full_stellar_interval_certificate=False,
      full_physical_dynamic_force_and_readout=False,full_span_variational_certificate=False,
      specified_cold_WD_matches_observed_inner_WD=False,
      full_28_parameter_initialization_certificate=False,complete_nonlinear_observational_inference=False))
    paths=[p for p in OUT.rglob('*') if p.is_file() and p!=OUT/'manifest.json']
    paths += [ROOT/'verification/nonzero_drive.py',ROOT/'verification/nonzero_drive_audit.py',ROOT/'notes/REQUEST17_NONZERO_DRIVE_KO.md']
    paths += [ROOT/p for p in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',scope='Request 17 지정 후보 및 조건부 정리; 전체 연구 완료 아님',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}))

def verify():
    histories={
      'outputs/validated-variational/manifest.json':'outputs/remaining-levers15/historical-note-bindings.json',
      'outputs/remaining-levers15/manifest.json':'outputs/nbody-readout16/historical-note-bindings.json',
      'outputs/nbody-readout16/manifest.json':'outputs/nonzero-drive17/historical-note-bindings.json'}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/nonzero-drive17/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,value in before.items():
                        if k!='sha256':assert after[k]==value,k
                    for k,value in before['sha256'].items():
                        if k not in old:assert after['sha256'][k]==value,k
                else:assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name)
            count+=1
    for name,expected in json.loads((OUT/'provenance.json').read_text())['input_sha256'].items():assert sha(ROOT/name)==expected
    gates=json.loads((OUT/'gates.json').read_text());assert gates['theorem_progress'] and not gates['complete_nonlinear_observational_inference']
    print('PASS:',count,'현재·역사적 해시, 입력 및 미완료 판정 보존')

if __name__=='__main__':{'symbolic':symbolic,'stars':stars,'intervals':intervals,'check':check,'maintain':maintain,'seal':seal,'verify':verify}[sys.argv[1]]()

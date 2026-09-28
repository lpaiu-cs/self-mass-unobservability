"""Request29 independent saved-state checks and fresh native replay."""
import json, sys, shutil
from pathlib import Path
import numpy as np
import common_eos as c
from thermal_wd import mesa


def audit():
    root=c.OUT;native_rows=[];states=[];sources=[];eos=c.EOS()
    initial=dict(np.load(root/'restored-state-17-4.npz'))
    for path in sorted(root.glob('*-native.npz')):
        label=path.name.removesuffix('-native.npz');d=np.load(path)
        inp=np.load(root/(label+'-input.npz'));_,p=mesa(root/(label+'-profile.data.gz'))
        assert all(np.all(np.isfinite(d[k])) for k in d.files),label
        errors=[float(max(abs(np.log(d['T'])-inp['lnT']))),
            float(max(abs(np.log(d['rho'])-inp['lnd']))),float(np.max(abs(d['X']-inp['X'])))]
        assert max(errors)<1e-12,(label,errors)
        assert np.array_equal(d['heat'],p['eps_nuc']) and np.array_equal(d['neutrino'],p['eps_nuc_neu_total']),label
        native_rows.append(dict(label=label,calls=len(d['T']),input_errors=errors))
    for path in sorted(root.glob('*-state-*.npz')):
        if path.name.startswith(('base-','epoch0-')): continue
        d=np.load(path);label=path.name.removesuffix('.npz')
        assert np.array_equal(d['dm'],initial['dm']),label
        assert d['X'].min()>=0 and max(abs(d['X'].sum(1)-1))<1e-12,label
        f=1-2*d['m_mid_geom']/d['r_mid_m'];assert f.min()>0
        states.append(dict(label=label,min_X=float(d['X'].min()),sum_X_error=float(max(abs(d['X'].sum(1)-1))),
            min_metric_f=float(f.min()),max_composition_change=float(np.max(abs(d['X']-initial['X'])))))
    Linf=np.r_[initial['L']*np.exp(2*initial['nu_faces'][:-1]),0.]
    rest=(c.W/c.A-1)*(c.gr.C*100)**2
    for path in sorted(root.glob('*evolve-*-source.npz')):
        label=path.name.removesuffix('-source.npz');out=dict(np.load(path))
        report=json.loads((root/(label+'-source.json')).read_text())
        parts=label.split('-');step=int(parts[-1]);count=int(parts[-2]);iso=bool(report.get('isobaric',False))
        prefix='P-' if iso else ''
        start=initial if step==0 else dict(np.load(root/f'{prefix}evolve-{count}-{step-1}-state-17-4.npz'))
        h=report['dt']*np.exp(start['nu'])
        qflux=(Linf[:-1]-Linf[1:])/(start['dm']*np.exp(2*start['nu']))
        defects=[];scores=[]
        for i in range(len(h)):
            mode=1 if iso else 2;coordinate=start['logP'][i] if iso else start['lnd'][i]
            a=eos(mode,coordinate,start['lnT'][i],start['X'][i])
            b=eos(mode,coordinate,out['lnT'][i],out['X'][i])
            e0=a[2]+(a[1]/a[0] if iso else 0);e1=b[2]+(b[1]/b[0] if iso else 0)
            target=e0-rest@(out['X'][i]-start['X'][i])-out['loss_per_baryon_gram'][i]-h[i]*qflux[i]
            budget=max(2.,32*np.spacing(abs(target)),abs(target-e0)*1e-8)
            defects.append(e1-target);scores.append(abs(e1-target)/budget)
        defects=np.array(defects);scores=np.array(scores)
        np.savez_compressed(root/(label+'-independent-energy.npz'),defects=defects,scores=scores)
        sources.append(dict(label=label,max_local_score=float(max(scores)),
            local_gate_passed=bool(max(scores)<=1),
            redshifted_absolute_residual_over_release=float(start['dm']@(np.exp(start['nu'])*abs(defects))/report['absolute_rest_release_erg'])))
        print('INDEPENDENT ENERGY',sources[-1],flush=True)
    c.save('independent-audit.json',dict(classification='Counterexample candidate',native=native_rows,
        total_native_calls=sum(r['calls'] for r in native_rows),states=states,source_energy=sources,
        scope='Fresh EOS evaluation of stored endpoints and independent arithmetic checks. '
            'EOS root warm-start variation is included; failures are reported without widening the registered local tolerance.'))
    print('AUDIT',len(native_rows),'native states;',len(states),'GR states;',len(sources),'source endpoints',flush=True)


def recheck():
    import native_eos_bridge as shim
    c.native.OUT=c.OUT;c.native.CACHE=c.CACHE;c.native.context()
    data=dict(np.load(c.OUT/'restored-state-17-4.npz'));rows=[]
    for tag,reference in [('pass','EOS-pass-through'),('common','EOS-common')]:
        label='EOS-recheck-'+tag
        assert not (c.CACHE/label).exists(),label
        c.native.setup(label,data,species=c.NAMES,network=(c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
        aux=None if tag=='pass' else np.load(c.OUT/'EOS-common-replacement.npy')
        if aux is not None: np.save(c.OUT/(label+'-replacement.npy'),aux)
        shim.trace(label,aux)
        a=np.load(c.OUT/(label+'-native.npz'));b=np.load(c.OUT/(reference+'-native.npz'))
        errors={k:float(np.max(abs(a[k]-b[k])/np.maximum(1e-100,abs(b[k])))) for k in a.files}
        assert max(errors.values())<1e-10,(label,errors)
        rows.append(dict(label=label,reference=reference,errors=errors))
    c.save('fresh-recheck.json',dict(classification='Proven',rows=rows,
        scope='Two independent disposable-process executions using the identical executable and inputs. Reproducibility, not independent physical validation.'))


def derivative_coordinates():
    base=np.load(c.OUT/'EOS-common-native.npz');h=5e-5
    minus=np.load(c.OUT/'EOS-attribution-frozen--1-native.npz')
    plus=np.load(c.OUT/'EOS-attribution-frozen-1-native.npz')
    dx=(plus['X']-minus['X'])/(2*h);dr=(plus['rho']-minus['rho'])/(2*h)
    dt=(plus['T']-minus['T'])/(2*h)
    fd=(plus['dxdt']-minus['dxdt'])/(2*h)
    chain=np.einsum('nij,nj->ni',base['jacobian'],dx)+base['dxdt_rho']*dr[:,None]+base['dxdt_T']*dt[:,None]
    scale=np.maximum(1e-30,np.maximum(abs(base['dxdt_T']*base['T'][:,None]).max(1),abs(base['dxdt']).max(1)))
    score=abs(fd-chain)/scale[:,None];mask=abs(base['heat'])>1
    where=np.unravel_index(np.argmax(np.where(mask[:,None],score,-1)),score.shape);i,j=where
    c.save('derivative-coordinate-audit.json',dict(classification='Counterexample candidate',
        max_input_X_difference=float(np.max(abs(plus['X']-minus['X']))),
        max_input_density_log_difference=float(max(abs(np.log(plus['rho']/minus['rho'])))),
        native_coordinate_chain_error=float(score[i,j]),cell=int(i),species=c.NAMES[j],
        finite_difference=float(fd[i,j]),reported_coordinate_chain=float(chain[i,j]),
        input_composition_chain=float((base['jacobian'][i]@dx[i])[j]),
        scope='Audit actual native coordinate differences, including roundoff and composition normalization; no declaration of a corrected physical derivative.'))


def composition_diagnostics():
    original=np.load(c.OUT/'restored-state-17-4.npz')['X'];rows=[]
    pairs=[('raw-1-2','P-evolve-1-0','P-evolve-2-1'),('raw-2-4','P-evolve-2-1','P-evolve-4-3'),
        ('extrapolated','extrapolated-1-2','extrapolated-2-4')]
    for label,left,right in pairs:
        a=np.load(c.OUT/(left+'-state-17-4.npz'))['X'];b=np.load(c.OUT/(right+'-state-17-4.npz'))['X']
        score=abs(b-a)/(1e-16+1e-3*abs(b-original));index=np.unravel_index(np.argmax(score),score.shape)
        i,j=index;copy=score.copy();copy[:,c.NAMES.index('he4')]=0
        second=np.unravel_index(np.argmax(copy),copy.shape)
        rows.append(dict(label=label,cell=int(i),species=c.NAMES[j],score=float(score[i,j]),
            difference=float(b[i,j]-a[i,j]),total_change=float(b[i,j]-original[i,j]),
            max_except_He4=float(copy.max()),except_He4_cell=int(second[0]),except_He4_species=c.NAMES[second[1]],
            failures_by_species={n:int(np.sum(score[:,j]>1)) for j,n in enumerate(c.NAMES)}))
    c.save('composition-diagnostics.json',dict(classification='Counterexample candidate',rows=rows))
    print('COMPOSITION DIAGNOSTICS',rows,flush=True)


def provenance():
    for rel in ['star/private/net.f90','net/private/net_eval.f90']:
        dest=c.OUT/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(c.fresh.MESA/rel,dest)
    paths=[c.fresh.BINARY,c.CACHE/'common_eos_aux.so',c.gr.CACHE/'build/src/libfree_eos.so.1.0.0',
        c.gr.CACHE/'gr_eos_bridge.so']
    c.save('provenance.json',dict(classification='Imported from prior work',
        executable_and_libraries={str(p):c.sha(p) for p in paths},
        MESA_source_root=str(c.fresh.MESA),MESA_source_version=(c.fresh.MESA/'data/version_number').read_text().strip(),
        source_caveat='Archived source signatures are evidence, not proof of compilation into the preserved binary. '
            'Both eta derivative slots gave unchanged native outputs in the functional controls.',
        FreeEOS_source_sha256={p.relative_to(c.gr.SOURCE).as_posix():c.sha(p) for p in c.gr.SOURCE.rglob('*') if p.is_file()},
        unchanged_native_data_binding='outputs/fresh-microphysics25/runtime-data-bindings.json'))


def mass_quadrature():
    assert (c.OUT/'mass-quadrature-plan.json').exists()
    descriptions=[('restored','restored-input.npz'),('P-evolve-4-3','P-evolve-4-3-source.npz')]
    models=[]
    for label,filename in descriptions:
        data=dict(np.load(c.OUT/filename));solver=c.Structure(label,data,17,4)
        pars=json.loads((c.OUT/(label+'-structure-17-4.json')).read_text())['parameters']
        residual,inner,outer=solver.branches(pars,record=True)
        assert max(abs(residual))<1e-8
        models.append((solver,inner,outer))
    original=np.load(c.OUT/'restored-state-17-4.npz');final=np.load(c.OUT/'P-evolve-4-3-state-17-4.npz')
    summary=json.loads((c.OUT/'P-evolution.json').read_text())['rows'][-1]
    dcx=((final['X']-original['X'])/c.A)@(c.W-c.A);cx=(original['X']/c.A)@c.W
    rows=[]
    for count in [2,4]:
        nodes,weights=np.polynomial.legendre.leggauss(count);terms=[]
        for i in range(len(original['dm'])):
            values=[]
            for solver,inner,outer in models:
                m=solver.mat;outside=i<m.split
                start=outer[i] if outside else inner[len(m.lp)-1-i]
                low,high=(m.outer[i],m.outer[i+1]) if outside else (m.inner[i+1],m.inner[i])
                fractions=(low+high)/2+(high-low)*nodes/2;sample=[]
                for q in fractions:
                    y=solver.step(np.log(start[0]),np.log(q),start[1:],i,m.B,outside)
                    a,_,_=c.be.invert(solver.eos,y[2],solver.ref[i,3],m.eps[i],m.lt[i])
                    sample.append([y[0]*m.R,y[1]*m.B,a[2]])
                values.append(np.array(sample))
            a,b=values;r0,m0,u0=a.T;r1,m1,u1=b.T
            s0=np.sqrt(1-2*m0/r0);s1=np.sqrt(1-2*m1/r1)
            ds=-2*((m1-m0)-m0/r0*(r1-r0))/r1/(s0+s1)
            integrand=(dcx[i]+(u1-u0)/(c.gr.C*100)**2)*s1+(cx[i]+u0/(c.gr.C*100)**2)*ds
            terms.append(float(weights@integrand/2))
        measured=float(original['dm']@np.array(terms))*(c.gr.C*100)**2
        score=abs(measured-summary['expected_mass_change_energy_erg'])/summary['release_erg']
        row=dict(points=count,mass_change_energy_erg=measured,energy_score=score,energy_passed=bool(score<1e-6))
        rows.append(row);np.save(c.OUT/f'mass-quadrature-{count}.npy',terms)
        print('MASS QUADRATURE',row,flush=True)
    c.save('mass-quadrature.json',dict(classification='Counterexample candidate',rows=rows,
        refinement_over_release=abs(rows[1]['mass_change_energy_erg']-rows[0]['mass_change_energy_erg'])/summary['release_erg'],
        midpoint_score=summary['energy_score'],
        scope='Gauss integration of the same declared piecewise-isentropic TOV models at identical baryon coordinates. '
            'This independently checks the midpoint mass-change estimator, not the full physical dynamics.'))


def finalize():
    assert (c.OUT/'independent-audit.json').exists()
    for filename in ['evolution-resolution.log','input-attribution.log','mass-quadrature.log','isobaric.log','extrapolation.log','recheck.log']:
        shutil.copy2(c.CACHE/filename,c.OUT/filename)
    c.save('failure-register.json',dict(classification='Counterexample candidate',
        fixed_volume_initial_inverse=dict(cell=1875,defect_erg_g=-2,completed=False),
        fixed_volume_after_32ulp_plan=dict(steps=4,step_zero_based=0,cell=3587,defect_erg_g=19,completed=False),
        original_derivative_slot_classifier='Assertion stopped at unchanged outputs for slot 15. Both slots subsequently recorded as unidentified; no order assignment.',
        source_precision='Independent endpoint EOS re-evaluation is checked against the unchanged registered local inverse tolerance, separately from the global energy criterion.',
        original_verdicts_preserved=True))
    independent=json.loads((c.OUT/'independent-audit.json').read_text())
    mass=json.loads((c.OUT/'mass-quadrature.json').read_text());p=json.loads((c.OUT/'P-evolution.json').read_text())
    ext=json.loads((c.OUT/'extrapolation.json').read_text())
    c.save('gates.json',dict(classification='Counterexample candidate',
        finite_initial_EOS_derivative_checks=json.loads((c.OUT/'full-EOS-audit.json').read_text())['passed'],
        physical_common_EOS_certified=False,continuous_EOS_derivative_error_certified=False,
        native_common_EOS_input_intervention_run=True,native_analytic_derivative_gate=False,
        fixed_baryon_GR_initialization=True,fixed_volume_time_sequence_completed=False,
        prescribed_flux_isobaric_time_sequence_completed=len(p['rows'])==3,
        isobaric_four_step_energy_gate=p['rows'][-1]['energy_passed'],
        independent_mass_quadrature_energy_gate=all(r['energy_passed'] for r in mass['rows']),
        isobaric_composition_time_gate=p['rows'][-1]['refinement']['passed'],
        Richardson_all_gates=ext['passed'],
        independent_local_EOS_inverse_gate=all(r['local_gate_passed'] for r in independent['source_energy']),
        whole_star_physical_transport=False,full_GR_thermal_fluid_metric_evolution=False,
        conditional_mass_normalized_scalar_readout=True,scalar_actual_drive_charge_map=False,
        complete_nonlinear_observation=False,final_submission_package_updated=False))
    additions={
        'model-definition':'분류: Counterexample candidate. 불소를 복원한 26종 핵 재고를 보존하고, 이온 수·전하 수 보존 사상과 조성별 원자 결합에너지 기준 이동으로 FreeEOS 입력을 연결했다. 모든 초기 5,735구역의 유한 미분 대조는 통과했다. 미지원 원소의 부분 이온화·혼합 엔트로피와 물리 EOS 오차 인증은 남는다.\n\n분류: Proven. 고정 조성의 합성 자유에너지 F_B=C F_FE(Cρ_B,T,ε)+ΔI는 원래 자유에너지의 열역학 항등식을 보존한다. 조성 변화에는 화학 에너지 항이 필요하다.',
        'observable-targets':'분류: Proven. 등방 에너지 방출의 단극 운동량 변화가 u^μ에 평행하면 정규화 투영 후 횡방향 가속도는 0이다. 순수 GR의 이 단극 경계에서 핵 가열·질량 손실만으로 새 자유낙하 힘을 얻지 못한다.\n\n분류: Counterexample candidate. β=−4의 지정 GR 보간 모형과 고정 광도 열 경로에서 질량 정규화 scalar 읽기를 계산했다. 4분할의 δ(α_A/φ∞)는 약 −5.02184e−11이다. 유한 Wronskian 식은 외부 꼬리·계량·질량 정규화를 포함한다. 실제 구동·비선형 역반응·관측 likelihood 및 이 미소 잔여항의 물리적 오차 인증은 아니다.',
        'adiabatic-limit':'분류: Proven. 같은 핵 반응의 정지에너지 손실과 열량을 별도 총질량 원천으로 중복 계수할 수 없다. 고정 바리온 열 경로에는 화학 에너지, 적색편이된 중성미자·광자 손실과 표면 압력 일을 포함해야 한다.\n\n분류: Conjectural. 표면 광도가 정해졌다고 내부 수송이 닫히거나 단일 orbital relaxation 상태가 식별되는 것은 아니다. 영 scalar 가지와 정적 carrier 보간 경계는 그대로 유지한다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 5,735개 물질 구역에서 1.6294일의 26종 반응·열·동일 바리온 TOV 경로를 계산했다. 고정 압력 4분할은 에너지와 온도 대조를 통과했으나 조성 시간 오차/허용량 3.37977로 실패했다. 별도 외삽은 Be7 오차를 줄였지만 H1·He4의 수 ulp 차이로 엄격 조성 기준을 실패했다.\n\n분류: Conjectural. 고정 광도·준정적 계량 경로는 자체 수송·유체·계량 진화가 아니다. 조건부 scalar 읽기와 실제 주파수 구동의 phase lag 또는 pole 식별은 구별한다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 공통 EOS의 초기 유한 미분 검사는 통과했지만 반응 함수의 반환 온도 미분은 EOS 보조 입력을 고정해도 벡터 척도 약 0.107의 차이를 보였다. 입력 좌표 차이는 0이므로 단순 입력 정규화로 설명하지 않는다. 두 eta 미분 인자 대조는 변화가 없어 순서를 판별하지 못했다.\n\n분류: Counterexample candidate. 고정 부피 4분할의 첫 원천 단계는 에너지 역산 +19 erg/g에서 중단했다. 완료된 1·2분할은 에너지·시간 기준을 실패했다. 별도 고정 압력 4분할과 외삽도 엄격 조성 시간 기준을 실패했다. 저장 단계의 독립 EOS 역산 잔차는 전역 질량 에너지 점수와 별도로 판정한다. 원래 실패를 변경하지 않는다.\n\n분류: Proven. 등방 GR 단극 질량 손실만으로는 횡방향 힘을 만들지 않는다.\n\n분류: Conjectural. 최소 추가 연결은 미지원 원소의 물리 EOS·연속 오차 보증, 정확한 반응 미분과 보존적 시간 오차 제어, 자체 수송·유체·계량·대기, 비영 scalar 구동 및 전체 관측 전방 모형이다.'}
    bindings={}
    for rel,digest in json.loads((c.OUT/'historical-bindings.json').read_text()).items():
        snap=c.OUT/'previous-notes'/Path(rel).name
        assert c.sha(snap)==digest
        bindings[rel]=dict(snapshot=snap.relative_to(c.ROOT).as_posix(),sha256=digest,historical_manifest=rel.startswith('paper/'))
    c.save('historical-note-bindings.json',bindings)
    for name,paragraph in additions.items():
        path=c.ROOT/'docs'/(name+'.md');prior=c.OUT/'previous-notes'/path.name
        assert path.read_bytes()==prior.read_bytes(),name
        text='\n\n## Request 29 공통 EOS와 GR 열 경로의 관측 연결\n\n'+paragraph+'\n\n세부 근거: [한글 보고서](../notes/REQUEST29_COMMON_EOS_GR_KO.md).\n'
        with path.open('ab') as f: f.write(text.encode('utf-8'))
    manifest_path=c.ROOT/'paper/revision-manifest.json';paper=json.loads(manifest_path.read_text())
    paper['request29_supporting_note_update']=dict(evidence_manifest='outputs/common-eos29/manifest.json',
        historical_notes='outputs/common-eos29/historical-note-bindings.json',
        status='26종 EOS 유한 감사·같은 바리온 GR 열 경로와 질량 정규화 scalar 읽기; 미분·시간 실패 보존, 물리 EOS·실제 GR/관측 폐쇄 미완료',
        artifact_status='Request12 원고 PDF/ZIP은 역사 산출물로 보존한다.')
    for name in additions:
        rel='docs/'+name+'.md';paper['sha256'][rel]=c.sha(c.ROOT/rel)
    manifest_path.write_text(json.dumps(paper,ensure_ascii=False,indent=2)+'\n')
    pack()


def pack():
    paths=[p for p in c.OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [c.ROOT/'verification'/n for n in ['common_eos.py','verify_common_eos.py','common_eos_bridge.f90','native_eos_bridge.py']]
    paths += [c.ROOT/'notes/REQUEST29_COMMON_EOS_GR_KO.md',c.ROOT/'.gitattributes']
    paths += [c.ROOT/n for n in json.loads((c.OUT/'historical-note-bindings.json').read_text())]
    c.save('manifest.json',dict(classification='Proven',sha256={p.relative_to(c.ROOT).as_posix():c.sha(p) for p in sorted(paths)}))


def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27','conservative-cell28','common-eos29']
    history={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*history,'paper/revision-manifest.json','outputs/common-eos29/manifest.json']:
        old=json.loads((c.ROOT/history[label]).read_text()) if label in history else {}
        for name,expected in json.loads((c.ROOT/label).read_text())['sha256'].items():
            path=c.ROOT/name
            if name in old:
                bind=old[name];path=c.ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((c.ROOT/name).read_text())
                    for k,v in before.items():
                        if k!='sha256': assert after[k]==v,k
                    for k,v in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==v,k
                else: assert (c.ROOT/name).read_bytes().startswith(path.read_bytes())
            assert c.sha(path)==expected,(label,name);count+=1
    prov=json.loads((c.OUT/'provenance.json').read_text())
    for path,digest in prov['executable_and_libraries'].items(): assert c.sha(Path(path))==digest,path
    for path,digest in prov['FreeEOS_source_sha256'].items(): assert c.sha(c.gr.SOURCE/path)==digest,path
    binding=json.loads((c.ROOT/prov['unchanged_native_data_binding']).read_text())
    for path,digest in binding['sha256'].items(): assert c.sha(Path(binding['root'])/path)==digest,path
    gates=json.loads((c.OUT/'gates.json').read_text())
    for key in ['physical_common_EOS_certified','continuous_EOS_derivative_error_certified','native_analytic_derivative_gate',
        'fixed_volume_time_sequence_completed','isobaric_composition_time_gate','Richardson_all_gates',
        'whole_star_physical_transport','full_GR_thermal_fluid_metric_evolution','scalar_actual_drive_charge_map',
        'complete_nonlinear_observation','final_submission_package_updated']: assert gates[key] is False,key
    print('PASS',count,'artifact/history SHA;',len(binding['sha256']),'native data SHA;',
        len(prov['FreeEOS_source_sha256']),'FreeEOS source SHA',flush=True)


def git_blobs():
    c.cell.OUT=c.OUT;c.cell.CACHE=c.CACHE
    # Reuse the same raw-object checker, with only the manifest filename
    # corrected for this phase; no Git text conversion is accepted.
    import hashlib,subprocess
    expected=dict(json.loads((c.OUT/'manifest.json').read_text())['sha256'])
    expected['outputs/common-eos29/manifest.json']=c.sha(c.OUT/'manifest.json')
    with subprocess.Popen(['git','cat-file','--batch'],cwd=c.ROOT,stdin=subprocess.PIPE,stdout=subprocess.PIPE) as proc:
        for rel,digest in expected.items():
            proc.stdin.write(('HEAD:'+rel+'\n').encode());proc.stdin.flush()
            header=proc.stdout.readline().split();assert len(header)==3 and header[1]==b'blob',(rel,header)
            remaining=int(header[2]);actual=hashlib.sha256()
            while remaining:
                chunk=proc.stdout.read(min(remaining,1048576));assert chunk
                actual.update(chunk);remaining-=len(chunk)
            assert proc.stdout.read(1)==b'\n' and actual.hexdigest()==digest,rel
        proc.stdin.close();assert proc.wait(timeout=10)==0
    print('PASS',len(expected),'raw Git blobs at',subprocess.check_output(['git','rev-parse','HEAD'],cwd=c.ROOT,text=True).strip(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

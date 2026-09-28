"""남은 레버의 격리된 빌드·검증 실행. WSL에서 실행한다."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import shlex
import os
import difflib
import re

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/remaining-levers15'
RUN=Path.home()/'work/nutimo_pilot'
CAPD=RUN/'request14_capd'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def multiprecision_build():
    OUT.mkdir(exist_ok=True,parents=True)
    deps=RUN/'request15_deps';deps.mkdir(exist_ok=True)
    with (OUT/'mp-dependencies.log').open('w') as log:
        subprocess.run(['apt-get','download','libmpfr-dev','libgmp-dev'],cwd=deps,stdout=log,stderr=subprocess.STDOUT,check=True)
    for p in deps.glob('*.deb'):
        subprocess.run(['dpkg-deb','-x',str(p),str(deps)],check=True)
    build=CAPD/'build-request15-mp'
    include=f'-I{deps}/usr/include -I{deps}/usr/include/x86_64-linux-gnu'
    commands=[['cmake','-S',str(CAPD),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release','-DCAPD_ENABLE_MULTIPRECISION=ON',f'-DCMAKE_CXX_FLAGS={include}',f'-DCMAKE_LIBRARY_PATH={deps}/usr/lib/x86_64-linux-gnu'],['cmake','--build',str(build),'-j','4']]
    with (OUT/'mp-build.log').open('w') as log:
        for cmd in commands:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'mp-build.json').write_text(json.dumps(dict(commands=commands,capd_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=CAPD,text=True).strip(),packages={p.name:sha(p) for p in deps.glob('*.deb')},library_sha256=sha(build/'libcapd.a')),indent=2)+'\n')


def mp_client():
    flags=shlex.split(subprocess.check_output([str(CAPD/'build-request15-mp/bin/capd-config'),'--cflags','--libs'],text=True))
    deps=RUN/'request15_deps/usr'
    cmd=['g++',str(ROOT/'verification/remaining_variational.cpp'),f'-I{deps}/include',f'-I{deps}/include/x86_64-linux-gnu',f'-L{deps}/lib/x86_64-linux-gnu',*flags,'-o',str(RUN/'request15_variational_v2')]
    with (OUT/'mp-client-build.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'mp-client-build.json').write_text(json.dumps(dict(command=cmd,binary_sha256=sha(RUN/'request15_variational_v2'),source_sha256=sha(ROOT/'verification/remaining_variational.cpp')),indent=2)+'\n')

def mp_run():
    label,horizon,bits,seconds=sys.argv[2:6]
    cmd=[str(RUN/'request15_variational_v2'),str(ROOT/'outputs/validated-variational/ivp.hex'),horizon,bits,str(OUT/(label+'.jsonl')),seconds,*sys.argv[6:]]
    with (OUT/(label+'.log')).open('w') as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    (OUT/(label+'-run.json')).write_text(json.dumps(dict(command=cmd,exit_code=result.returncode),indent=2)+'\n')

def init_build():
    source=RUN/'nutimo_request13_stable/src'
    target=RUN/'nutimo_request15_initialization/src'
    assert not target.exists(), 'Preserve existing runtime evidence'
    shutil.copytree(source,target)
    p=target/'Parameters.cpp';old=p.read_text();new=old
    replacements={
      'Mp * nrpt2 * rp + Mi * nrit2 * rit - GMsol': 'Mp * nrpt2 * rp + Mi * nrit2 * ri - GMsol',
      'rpt = rpt + rBit - un/(deux* pow( clight, 2) ) * ( deux * dotprod3d<value_type>(rpt, rBit) - GMsol * Gg[0][1] * Mi / rip * dotprod3d<value_type>( rBit, nip ) * nip ) ;':
      'rpt = rpt + rBit; // Request 15: published osculating initialization convention (Voisin 2020, footnote 8).',
      'rit = rit + rBit - un/(deux* pow( clight, 2) ) * ( deux * dotprod3d<value_type>(rit, rBit) - GMsol * Gg[0][1] * Mp / rip * dotprod3d<value_type>( rBit, nip ) * nip ) ;':
      'rit = rit + rBit; // Same velocity-addition convention; subsequent COM correction retained.'}
    for a,b in replacements.items():
        assert new.count(a)==1,a
        new=new.replace(a,b)
    p.write_text(new)
    (OUT/'initialization.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='request13_stable/Parameters.cpp',tofile='request15_initialization/Parameters.cpp')))
    snapshots=OUT/'native-source';snapshots.mkdir(exist_ok=True)
    for name in ['Parameters.cpp','Orbital_elements.cpp','Diagnostics.cpp','Delay_brut.cpp','Fittriple-compute.cpp','Fittriple-init.cpp','Fittriple-IO.cpp','Spline.cpp']:
        shutil.copy2(source/name,snapshots/name)
    shutil.copy2(p,OUT/'Parameters-corrected.cpp')
    command=json.loads((ROOT/'outputs/research-remediation/stable-build.json').read_text())['command']
    with (OUT/'initialization-build.log').open('w') as log:
        subprocess.run(command,cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'initialization-build.json').write_text(json.dumps(dict(command=command,source_sha256=sha(source/'Parameters.cpp'),corrected_source_sha256=sha(p),library_sha256=sha(target/'libFittriplecpp.so')),indent=2)+'\n')

def init_launch():
    label=sys.argv[2]
    source=RUN/('nutimo_request13_stable/src' if label=='original' else 'nutimo_request15_initialization/src')
    run=RUN/('run_request15_'+label)
    assert not run.exists(), 'Do not overwrite an evaluation'
    shutil.copytree(RUN/'run_request13_stable',run)
    env=dict(os.environ,LD_LIBRARY_PATH=str(source),PYTHONPATH=str(source),TEMPO2=str(RUN/'install/third_party/tempo2'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    for k in ['SEPDYN_A','SEPDYN_W','SEPDYN_PH','SEPDYN_TAU']:env.pop(k,None)
    action={'remapped':'state-refit','spin-calibrated':'spin-calibrate','reload':'reload'}.get(label,'init-evaluate')
    with (OUT/(label+'-initialization.log')).open('w') as log:
        subprocess.run([sys.executable,str(Path(__file__).resolve()),action,label],cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)

def state_refit():
    sys.path.insert(0,str(RUN/'request13_deps'))
    import numpy as np
    from scipy.optimize import least_squares
    import python_Fittriple_interface as pfi
    fit=pfi.PyFittriple('parfile-planetGR-max-bestfit','0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
    target=np.load(OUT/'original-initialization.npz')['state']
    fit.Compute_initial_state_vectors();mass=fit.Get_masses().copy()
    base=np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    scales=base['scales'][base['fmap']].astype(float)
    h=np.array(json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())['abs_steps'])
    # Fit the 22 orbital coordinates to 18 relative state components and 4 masses.
    # This is a coordinate remapping, not an observational likelihood fit.
    targetrel=target[1:]-target[0]
    divisor=np.array([1e3]*3+[1e-3]*3)
    calls=0
    def value(z):
        nonlocal calls
        zz=np.zeros(28);zz[6:]=z
        fit.Set_fitted_parameter_relativeshifts(h*zz/scales)
        st=fit.Compute_initial_state_vectors();m=fit.Get_masses()
        calls+=1
        return np.r_[((st[1:]-st[0]-targetrel)/divisor).ravel(),(m/mass-1)/1e-8]
    def jac(z):
        cols=[]
        for k in range(22):
            dz=np.zeros(22);dz[k]=1e-3
            cols.append((value(z+dz)-value(z-dz))/(2e-3))
        return np.column_stack(cols)
    result=least_squares(value,np.zeros(22),jac=jac,x_scale='jac',max_nfev=150,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    final=value(result.x);st=fit.Compute_initial_state_vectors();m=fit.Get_masses()
    fit.Compute_lnposterior(0);res=fit.Get_time_residuals().copy()
    old=np.load(OUT/'original-initialization.npz')['res'];dr=res-old
    np.savez(OUT/'remapped-initialization.npz',state=st,res=res,masses=m,target_masses=mass,scaled_orbital_shift=result.x,steps=h,parameters=fit.Get_parameters())
    report=dict(classification='Proven',scope='동일 물리 초기 상태로의 수치 좌표 재매핑; 관측 재적합이나 엄밀한 미분 인증 아님',optimizer_success=bool(result.success),message=result.message,calls=calls,max_scaled_state_mass_mismatch=float(np.max(abs(final))),relative_position_max_m=float(np.max(abs((st[1:]-st[0]-targetrel)[:,:3]))),relative_velocity_max_m_s=float(np.max(abs((st[1:]-st[0]-targetrel)[:,3:]))),mass_relative_max=float(np.max(abs(m/mass-1))),residual_rms_us=float(np.sqrt(np.mean(res*res))),difference_from_original_rms_us=float(np.sqrt(np.mean(dr*dr))),difference_from_original_max_us=float(np.max(abs(dr))),physical_inference_complete=False)
    (OUT/'remapped-initialization.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')

def readout_export():
    source=RUN/'nutimo_request15_initialization/src'
    target=RUN/'nutimo_request15_readout/src'
    assert not target.exists()
    shutil.copytree(source,target)
    path=target/'Delay_brut.cpp';old=path.read_text();new=old
    # Read-only exports at actual native call boundaries. Each mode writes once.
    for name,env,body in [
      ('Delays_Brut_nogeometric','TIMING15_NONGEOM',r'''
    static bool exported15=false;
    if (!exported15 && !einstein && shapiro && aberration && spinaxis==NULL) if(const char* file15=getenv("TIMING15_NONGEOM")) {
        FILE* f15=fopen(file15,"w");if(!f15)throw runtime_error("readout export");
        fprintf(f15,"%La %La %La %La %La %La %La %La\n",tis[0],Mp,Mi,Mo,freq,clight,deuxpi,daysec);
        for(int k=0;k<3;++k)fprintf(f15,"%La ",SSB_to_PSB[k]);
        for(int k=0;k<6;++k)fprintf(f15,"%La %La %La ",sp[0][k],si[0][k],so[0][k]);
        fprintf(f15,"%La\n",delay[0]);fclose(f15);exported15=true;
    }
'''),
      ('Delays_Brut_geometric_local','TIMING15_GEOM',r'''
    static bool exported15=false;
    if (!exported15 && kopeikin && shklovskii) if(const char* file15=getenv("TIMING15_GEOM")) {
        FILE* f15=fopen(file15,"w");if(!f15)throw runtime_error("readout export");
        fprintf(f15,"%La %La %La %La %La %La %La %La %La\n",tis[0],currentBAT,posepoch,distance,distance_derivative,clightdays,yrsec,clight,radmasdeg);
        for(int k=0;k<3;++k)fprintf(f15,"%La %La %La %La ",SSB_to_PSB[k],proper_motion[k],r_obs[k],sp[0][k]);
        fprintf(f15,"%La\n",delay[0]);fclose(f15);exported15=true;
    }
''')]:
        start=new.index('void '+name+'(')
        end=new.index('\n}',start)
        section=new[start:end]
        assert section.count('        return ;')==1
        section=section.replace('        return ;',body+'        return ;')
        new=new[:start]+section+new[end:]
    path.write_text(new)
    (OUT/'readout-export.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='initialization/Delay_brut.cpp',tofile='readout/Delay_brut.cpp')))
    shutil.copy2(path,OUT/'Delay-brut-export.cpp')
    shutil.copy2(source/'Constants.h',OUT/'native-source/Constants.h')
    command=json.loads((ROOT/'outputs/research-remediation/stable-build.json').read_text())['command']
    with (OUT/'readout-build.log').open('w') as log:subprocess.run(command,cwd=target,stdout=log,stderr=subprocess.STDOUT,check=True)
    run=RUN/'run_request15_readout';assert not run.exists();shutil.copytree(RUN/'run_request13_stable',run)
    env=dict(os.environ,LD_LIBRARY_PATH=str(target),PYTHONPATH=str(target),TEMPO2=str(RUN/'install/third_party/tempo2'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',TIMING15_NONGEOM=str(OUT/'nongem-native.hex'),TIMING15_GEOM=str(OUT/'geom-native.hex'))
    for k in ['SEPDYN_A','SEPDYN_W','SEPDYN_PH','SEPDYN_TAU']:env.pop(k,None)
    code="import python_Fittriple_interface as p; f=p.PyFittriple('parfile-planetGR-max-bestfit','0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'); f.Compute_lnposterior(0)"
    with (OUT/'readout-run.log').open('w') as log:subprocess.run([sys.executable,'-c',code],cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'readout-build.json').write_text(json.dumps(dict(command=command,library_sha256=sha(target/'libFittriplecpp.so'),source_sha256=sha(path)),indent=2)+'\n')

def spin_calibrate():
    import numpy as np
    import python_Fittriple_interface as pfi
    fit=pfi.PyFittriple('parfile-planetGR-max-bestfit','0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
    base=np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    saved=np.load(OUT/'remapped-initialization.npz')
    z=np.zeros(28);z[6:]=saved['scaled_orbital_shift'];h=saved['steps'];scales=base['scales'][base['fmap']].astype(float)
    old=np.load(OUT/'original-initialization.npz')['res']
    def residual(shift):
        fit.Set_fitted_parameter_relativeshifts(h*shift/scales)
        fit.Compute_lnposterior(0)
        return fit.Get_time_residuals().copy()
    before=residual(z);cols=[]
    for k in [4,5]:
        dz=np.zeros(28);dz[k]=.01
        cols.append((residual(z+dz)-residual(z-dz))/.02)
    J=np.column_stack(cols);shift=np.linalg.lstsq(J,old-before,rcond=None)[0]
    z[4:6]+=shift;after=residual(z)
    fit.Save_parfile(str(OUT/'corrected-remapped.par'))
    np.savez(OUT/'spin-calibrated.npz',res=after,scaled_parameters=z,parameters=fit.Get_parameters(),spin_jacobian=J)
    result=dict(classification='Proven',scope='기존 잔차를 재현하는 초기화·스핀 좌표 보정. 관측 자료의 신규 최적 적합 아님.',rms_difference_before_us=float(np.sqrt(np.mean((before-old)**2))),rms_difference_after_us=float(np.sqrt(np.mean((after-old)**2))),max_difference_after_us=float(np.max(abs(after-old))),residual_rms_us=float(np.sqrt(np.mean(after**2))),spin_scaled_shift=shift.tolist(),fresh_native_verification=True,full_physical_inference_complete=False)
    (OUT/'spin-calibrated.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')

def reload_check():
    import numpy as np
    import python_Fittriple_interface as pfi
    fit=pfi.PyFittriple(str(OUT/'corrected-remapped.par'),'0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
    fit.Compute_lnposterior(0);res=fit.Get_time_residuals().copy()
    reference=np.load(OUT/'spin-calibrated.npz')['res'];difference=res-reference
    assert np.max(abs(difference))<1e-3, 'Saved parameter roundtrip must agree within the declared 1 ns check'
    (OUT/'reload-check.json').write_text(json.dumps(dict(classification='Proven',saved_parameter_reload_passed=True,declared_max_difference_us=1e-3,measured_max_difference_us=float(np.max(abs(difference))),measured_rms_difference_us=float(np.sqrt(np.mean(difference*difference))),par_sha256=sha(OUT/'corrected-remapped.par')),indent=2)+'\n')

def maintain():
    previous=OUT/'request14-notes';previous.mkdir(exist_ok=True)
    frozen=json.loads((ROOT/'outputs/validated-variational/manifest.json').read_text())['sha256']
    additions={
      'model-definition':'분류: Proven. 활성 초기화의 질량 중심 위치 항과 회전 공변성을 깨는 속도 결합을 수정했다. 기존 물리 초기 상태를 재현하는 궤도·스핀 좌표를 별도로 계산했으며, 이는 새로운 힘이나 동적 관측량이 아니다. 질량 함수와 Kepler 근의 구간 미분은 인증했으나 초기 상태 전체로의 미분 연쇄는 미완료다.',
      'observable-targets':'분류: Proven. 실제 호출 입력에서 기하 지연 및 Shapiro·수차 지연의 구간 편미분을 검증했다. 방출시각과 Einstein 변환의 연쇄법칙 및 잔차 기반 오차식을 도출했다. 분류: Conjectural. 전체 기간의 운동·Einstein 누적 적분·스플라인·Tempo2 의존성을 연결하고 네 번째 천체의 생략 지연을 포함하거나 상계로 정당화해야 한다.',
      'adiabatic-limit':'분류: Proven. 초기화 규약의 수정과 매개변수 재매핑은 새 완화 pole을 만들지 않는다. 이번 질량·Kepler·대수 관측식의 부분 인증은 기존 단일 pole의 단열 붕괴 조건을 변경하지 않는다. 분류: Conjectural. EOS 응답의 전체 이체 구동·힘·관측량 연결은 여전히 별도의 물리 조건이다.',
      'nonadiabatic-regime':'분류: Proven. 다중 정밀도 전파, 초기화 수정·재매핑, 질량·Kepler·대수 관측식의 미분 검증을 수행했다. 이번 계산은 동적 SEP가 없는 GR 기준선에 조건부이며, 전체 기간 타이밍 인증과 비단열 검출 주장을 완성하지 않는다.',
      'failure-ledger-dynamic-chi':'분류: Proven. 추가 실패 경계를 기록한다. (1) 수정 전 초기화에는 질량 중심 항의 차원 오류와 회전 비공변 속도 항이 있다. 수정·재매핑 후에도 기존 수치 결과를 새 물리 제약으로 재해석하지 않는다. (2) 4체 운동에 3체 Einstein·Shapiro 지연 API가 연결되어 있으므로 전체 4체 관측 모형이라는 승격은 허용되지 않는다. (3) 질량·Kepler·호출 지점 관측식의 부분 미분 인증은 전체 초기화·누적 적분·보간·역변환의 결합 오차를 대신하지 못한다. 다중 정밀도 실행의 실패·종료 원인은 새 보고서에 개별 보존한다.'}
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    mapping={}
    for stem,body in additions.items():
        rel='docs/'+stem+'.md';path=ROOT/rel
        assert sha(path)==frozen[rel], 'Preserve a changed document; do not overwrite'
        snap=previous/(stem+'.md');shutil.copy2(path,snap)
        extra='\n\n## Request 15 후속 검증\n\n'+body+'\n\n세부 근거: [한글 실행·검증 보고서](../notes/REQUEST15_REMAINING_LEVERS_KO.md).\n'
        with path.open('ab') as f:f.write(extra.encode())
        mapping[rel]=dict(snapshot=str(snap.relative_to(ROOT)).replace(os.sep,'/'),sha256=frozen[rel],prefix_bytes=snap.stat().st_size)
        revision['sha256'][rel]=sha(path)
    (OUT/'historical-note-bindings.json').write_text(json.dumps(mapping,indent=2)+'\n')
    revision['request15_supporting_note_update']=dict(evidence_manifest='outputs/remaining-levers15/manifest.json',historical_notes='outputs/remaining-levers15/historical-note-bindings.json',status='초기화 결함 수정·좌표 재매핑 및 질량·Kepler·대수 관측식 부분 인증; 전체 타이밍·물리 추론 미완료',artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')

def seal():
    from variational_audit import radius,sqrt_up
    from fractions import Fraction
    mp_runs={}
    for label in ['mp256-long','mp256-local-v2','mp256-local-small']:
        execution=json.loads((OUT/(label+'-run.json')).read_text())
        rows=[json.loads(line) for line in (OUT/(label+'.jsonl')).read_text().splitlines()] if (OUT/(label+'.jsonl')).exists() else []
        log=(OUT/(label+'.log')).read_text()
        assert '"requested_horizon_completed":false' in log or '"requested_horizon_completed":true' in log
        report=dict(execution=execution,requested_horizon_completed='"requested_horizon_completed":true' in log,recorded_steps=[r['step'] for r in rows],last_recorded_time_internal=rows[-1]['t'] if rows else None)
        if rows:
            last=rows[-1]
            report.update(state_midpoint_error_l2_upper=sqrt_up(sum((radius(p)**2 for p in last['state']),Fraction(0))),jacobian_midpoint_error_operator_upper=sqrt_up(sum((radius(p)**2 for p in last['jacobian']),Fraction(0))))
        mp_runs[label]=report
    (OUT/'mp-audit.json').write_text(json.dumps(dict(classification='Proven',runs=mp_runs,full_timing_certificate=False,note='시각은 실제 저장한 인증 끝점만 보고한다. 내부 구간 폭과 binary64로 바깥 반올림해 저장한 구간 폭은 다를 수 있다.'),ensure_ascii=False,indent=2)+'\n')
    command=json.loads((OUT/'initialization-build.json').read_text())['command']
    def source_closure(src):
        names={name for name in command if name.endswith('.cpp')}
        names.add('python_Fittriple_interface.pyx')
        pending=list(names)
        while pending:
            for name in re.findall(r'#include\s*["<]([^">]+)',(src/pending.pop()).read_text()):
                if (src/name).is_file() and name not in names:
                    names.add(name);pending.append(name)
        return sorted(names)
    provenance={}
    for label,folder in [('original','nutimo_request13_stable'),('corrected','nutimo_request15_initialization'),('readout','nutimo_request15_readout')]:
        src=RUN/folder/'src'
        interface=list(src.glob('python_Fittriple_interface*.so'));assert len(interface)==1
        provenance[label]=dict(library_sha256=sha(src/'libFittriplecpp.so'),interface_sha256=sha(interface[0]),all_compiled_sources={name:sha(src/name) for name in source_closure(src)})
    src=RUN/'nutimo_request13_stable/src'
    names=source_closure(src);snapdir=(OUT/'native-source').resolve()
    for p in snapdir.iterdir():
        assert p.resolve().parent==snapdir
        if p.is_file() and p.name not in names:p.unlink() # only this producer's new snapshots
    for name in names:shutil.copy2(src/name,snapdir/name)
    provenance['mp_v1']=dict(binary_sha256=sha(RUN/'request15_variational'),source_sha256=sha(OUT/'mp-v1-client.cpp'),build_record='mp-v1-build.json')
    provenance['mp_v2']=dict(binary_sha256=sha(RUN/'request15_variational_v2'),source_sha256=sha(ROOT/'verification/remaining_variational.cpp'),build_record='mp-client-build.json')
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    gates=dict(classification='Proven',initialization_defects_repaired=True,physical_state_coordinate_remapping_checked=True,mass_kernel_derivatives_certified=True,kepler_roots_and_tails_certified=True,algebraic_readout_partials_certified=True,implicit_time_formulas_verified=True,full_span_variational_certificate=False,full_28_parameter_initialization_certificate=False,complete_delay_integral_interpolation_certificate=False,full_EOS_binary_force_readout_matching=False,complete_nonlinear_observational_inference=False,theorem_progress=True,genuine_new_observable_established=False)
    (OUT/'gates.json').write_text(json.dumps(gates,indent=2)+'\n')
    # Request 14 also bound the then-current manuscript manifest. Preserve its
    # exact bytes and verify the later annotation separately, without editing it.
    bindings_path=OUT/'historical-note-bindings.json';bindings=json.loads(bindings_path.read_text())
    rel='paper/revision-manifest.json'
    if rel not in bindings:
        expected=json.loads((ROOT/'outputs/validated-variational/manifest.json').read_text())['sha256'][rel]
        blob=subprocess.check_output(['git','show','cb531c1:'+rel],cwd=ROOT)
        candidates=[blob,blob.replace(b'\n',b'\r\n')]
        exact=[b for b in candidates if hashlib.sha256(b).hexdigest()==expected]
        assert exact, 'Historical manuscript manifest must match its frozen hash'
        snap=OUT/'request14-notes/revision-manifest.json';snap.write_bytes(exact[0])
        bindings[rel]=dict(snapshot=str(snap.relative_to(ROOT)),sha256=expected,historical_manifest=True,source_commit='cb531c1')
        bindings_path.write_text(json.dumps(bindings,indent=2)+'\n')
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [ROOT/'verification'/n for n in ['remaining_levers.py','remaining_audit.py','remaining_variational.cpp']]
    paths += [ROOT/'notes/REQUEST15_REMAINING_LEVERS_KO.md',ROOT/'paper/revision-manifest.json']
    paths += [ROOT/n for n in json.loads((OUT/'historical-note-bindings.json').read_text())]
    manifest=dict(classification='Proven',before_task_checkpoint='24556e7',scope='Request 15 실행·근거 묶음; 전체 타이밍과 물리 추론의 완료를 뜻하지 않는다.',sha256={str(p.relative_to(ROOT)).replace(os.sep,'/'):sha(p) for p in sorted(paths)})
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')

def verify():
    historical=json.loads((OUT/'historical-note-bindings.json').read_text());count=0
    for label in ['outputs/research-remediation/manifest.json','outputs/validated-variational/manifest.json','paper/revision-manifest.json','outputs/remaining-levers15/manifest.json']:
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if label=='outputs/validated-variational/manifest.json' and name in historical:
                binding=historical[name];path=ROOT/binding['snapshot']
                assert sha(path)==binding['sha256']==expected
                if binding.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for key,value in before.items():
                        if key!='sha256':assert after[key]==value,key
                    for key,value in before['sha256'].items():
                        if key not in historical:assert after['sha256'][key]==value,key
                else:assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name)
            count+=1
    gates=json.loads((OUT/'gates.json').read_text())
    assert gates['theorem_progress'] and not gates['complete_nonlinear_observational_inference'] and not gates['full_span_variational_certificate']
    provenance=json.loads((OUT/'provenance.json').read_text())
    for label in ['original','corrected','readout']:
        assert provenance[label]['interface_sha256']==provenance['original']['interface_sha256']
        for name,expected in provenance[label]['all_compiled_sources'].items():
            path=OUT/'native-source'/name
            if label!='original' and name=='Parameters.cpp':path=OUT/'Parameters-corrected.cpp'
            if label=='readout' and name=='Delay_brut.cpp':path=OUT/'Delay-brut-export.cpp'
            assert sha(path)==expected,(label,name)
    for version,build_name in [('mp_v1','mp-v1-build.json'),('mp_v2','mp-client-build.json')]:
        assert provenance[version]['binary_sha256']==json.loads((OUT/build_name).read_text())['binary_sha256']
    print('PASS:',count,'현재·동결 해시; 이전 문서 원문 보존; 전체 인증 승격 금지 유지')

def init_evaluate():
    import numpy as np
    import python_Fittriple_interface as pfi
    label=sys.argv[2]
    fit=pfi.PyFittriple('parfile-planetGR-max-bestfit','0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')
    fit.Compute_lnposterior(0)
    states=fit.Compute_initial_state_vectors();res=fit.Get_time_residuals().copy()
    base=np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    # A mapping diagnostic only: these finite differences are never certificates.
    h=np.array(json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())['abs_steps'])
    scale=base['scales'][base['fmap']].astype(float)
    jac=[];half=[]
    for k in range(28):
        deriv=[]
        for factor in [1.,.5]:
            z=np.zeros(28);z[k]=h[k]*factor/scale[k]
            fit.Set_fitted_parameter_relativeshifts(z);plus=fit.Compute_initial_state_vectors()
            fit.Set_fitted_parameter_relativeshifts(-z);minus=fit.Compute_initial_state_vectors()
            deriv.append((plus-minus)/(2*factor)) # derivative w.r.t. abs_steps coordinates
        jac.append(deriv[0]);half.append(deriv[1])
    fit.Set_fitted_parameter_relativeshifts(np.zeros(28))
    restored=fit.Compute_initial_state_vectors()
    assert np.array_equal(states,restored)
    assert np.isfinite(res).all() and np.isfinite(jac).all()
    np.savez(OUT/(label+'-initialization.npz'),state=states,res=res,jacobian=np.array(jac),half_jacobian=np.array(half),steps=h)
    report=dict(label=label,toas=len(res),initial_state_max_fd_halfstep_change=float(np.max(abs(np.array(jac)-half))),mapping_derivative_certified=False,fixed_turns=True,restoration_exact=True)
    (OUT/(label+'-initialization.json')).write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':
    {'mp-build':multiprecision_build,'mp-client':mp_client,'mp-run':mp_run,'init-build':init_build,'init-launch':init_launch,'init-evaluate':init_evaluate,'state-refit':state_refit,'readout-export':readout_export,'spin-calibrate':spin_calibrate,'reload':reload_check,'maintain':maintain,'seal':seal,'verify':verify}[sys.argv[1]]()

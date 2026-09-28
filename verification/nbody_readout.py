"""Request 16: 격리된 다체 관측식 수정·실행·검증. WSL에서 실행."""
import difflib
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/nbody-readout16'
RUN=Path.home()/'work/nutimo_pilot'
SOURCE=RUN/'nutimo_request15_initialization/src'
TARGET=RUN/'nutimo_request16/src'
TIM='0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def build():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'build.json').exists(), 'Do not overwrite a completed build'
    if not TARGET.exists():shutil.copytree(SOURCE,TARGET)
    header=(SOURCE/'Delay_brut.h').read_text()
    begin=header.index('void Delays_Brut_nogeometric(')
    old_signature=header[begin:header.index(';',begin)+1]
    member_signature=old_signature.replace('Delays_Brut_nogeometric','Nbody_nogeometric_delays')
    extra=', const std::vector<std::vector<value_type>>* full_states, const value_type* extra_masses, int nextra, value_type state_length'
    new_signature=old_signature.replace('value_type * spinaxis = NULL','value_type * spinaxis').replace('\n                ) ;',extra+'\n                ) ;')
    patches=[]
    for name in ['Delay_brut.h','Delay_brut.cpp','Fittriple.h','Fittriple-compute.cpp','Fittriple-init.cpp']:
        path=TARGET/name;old=(SOURCE/name).read_text();new=old
        if name=='Delay_brut.h':
            new=new.replace('#include "Constants.h"','#include "Constants.h"\n#include <vector>')
            new=new.replace(old_signature,new_signature)
        elif name=='Fittriple.h':
            new=new.replace('    void Initialize() ;','    void Initialize() ;\n\n    '+member_signature)
        elif name=='Delay_brut.cpp':
            start=new.index('void Delays_Brut_nogeometric(');end=new.index('\n}',start)
            part=new[start:end]
            sigend=part.index('{')
            # The original implementation repeats its default spin-axis argument.
            part=new_signature[:-1]+'\n'+part[sigend:]
            anchor='    nanflag = false; // set to true only if a delay is nan in the tests below'
            assert part.count(anchor)==1
            part=part.replace(anchor,anchor+r'''
    // Request 16: common Einstein/Shapiro summation for all additional bodies.
    // No class fields are added, preserving the existing Python extension ABI.
    std::vector<value_type> extra_u(ntis,0),extra_s(ntis,0);
    const char* scale_env=getenv("TIMING16_EXTRA_SCALE");
    char* scale_end=NULL;
    const value_type scale=scale_env?strtold(scale_env,&scale_end):1.L;
    if(!std::isfinite(scale)||scale<0||(scale_env&&(*scale_end||scale_end==scale_env))) {nanflag=true;return;}
    if(nextra<0 || (nextra && (!full_states||!extra_masses||full_states->size()!=static_cast<size_t>(ntis)))) {nanflag=true;return;}
    if(nextra && (einstein||shapiro)) {
        for(long int ti=0;ti<ntis;++ti) {
            if((*full_states)[ti].size()!=static_cast<size_t>(6*(3+nextra))) {nanflag=true;return;}
            for(int b=0;b<nextra;++b) {
                const value_type mass=extra_masses[b]*scale;
                if(!std::isfinite(mass)||mass<0) {nanflag=true;return;}
                if(mass==0)continue;
                value_type d[3];
                for(int k=0;k<3;++k)d[k]=sp[ti][k]-(*full_states)[ti][9+3*b+k]*state_length;
                const value_type r=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
                const value_type arg=(r-dotprod3d<value_type>(d,SSB_to_PSB))/clight;
                if(!(r>0)||!std::isfinite(r)||(shapiro&&(!(arg>0)||!std::isfinite(arg)))) {nanflag=true;return;}
                if(einstein)extra_u[ti]+=mass/r*GMsol/clight2;
                if(shapiro)extra_s[ti]+=-2.L*4.92521372097374e-6L/daysec*mass*log(arg);
            }
        }
    }
''')
            anchor='            deindt[i] += undemi * nvp2 / pow( clight, 2 ) ;'
            assert part.count(anchor)==2
            part=part.replace(anchor,anchor+'\n            deindt[i] += extra_u[i];')
            anchor='            delay[i] += Shap ;';assert part.count(anchor)==1
            part=part.replace(anchor,'            delay[i] += Shap + extra_s[i];')
            # Exact hexadecimal input/output samples, including the extra body,
            # at spaced grid points. Kept separate from frozen previous evidence.
            anchor='        return ;';assert part.count(anchor)==1
            part=part.replace(anchor,r'''
    if(const char* prefix=getenv("TIMING16_EXPORT")) if(nextra==1 && (einstein||shapiro)) {
        std::string filename=std::string(prefix)+(einstein?"-ein.hex":"-shap.hex");
        FILE* fp=fopen(filename.c_str(),"w");if(!fp)throw std::runtime_error("Request16 export failed");
        fprintf(fp,"%ld %La %La %La %La %La %La\n",ntis,t0,scale,extra_masses[0],GMsol,clight,daysec);
        for(long int ti=0;ti<ntis;ti+=std::max(1L,ntis/64)) {
            fprintf(fp,"%ld %La",ti,tis[ti]);
            for(int k=0;k<3;++k)fprintf(fp," %La %La %La",sp[ti][k],(*full_states)[ti][9+k]*state_length,SSB_to_PSB[k]);
            fprintf(fp," %La %La\n",extra_u[ti],extra_s[ti]);
        }
        fclose(fp);
    }
        return ;''')
            new=new[:start]+part+new[end:]
            new='#include <vector>\n#include <fstream>\n#include <stdexcept>\n'+new
        else:
            new=new.replace('Delays_Brut_nogeometric(', 'Nbody_nogeometric_delays(')
            if name=='Fittriple-compute.cpp':
                # Rotate retained native states on the cached RA/DEC path as well
                # as sp/si/so: extra-body getters otherwise return the old frame.
                anchor='                Rotate_SSB_to_PSB(sp[i], pRA, pDEC);';assert new.count(anchor)==1
                new=new.replace(anchor,'''                for(int body=0;body<3+parameters.nextra;++body) {
                    Rotate_SSB_to_PSB(&states[i][3*body],pRA,pDEC);
                    Rotate_PSB_to_SSB(&states[i][3*body],parameters.RA,parameters.DEC);
                    Rotate_SSB_to_PSB(&states[i][3*(3+parameters.nextra)+3*body],pRA,pDEC);
                    Rotate_PSB_to_SSB(&states[i][3*(3+parameters.nextra)+3*body],parameters.RA,parameters.DEC);
                }
'''+anchor)
                implementation=member_signature.replace('void Nbody_', 'void Fittriple::Nbody_').replace(' = NULL','')[:-1]
                new+='\n'+implementation+r'''
{
    ::Delays_Brut_nogeometric(tis,t0,ntis,nt0,sp,si,so,SSB_to_PSB,Mp,Mi,Mo,freq,
        einstein,shapiro,aberration,delay,truefreq,moydeindt,nanflag,spinaxis,
        &states,parameters.nextra?&int_M_extra[0]:NULL,parameters.nextra,length);
}
'''
        path.write_text(new)
        patches+=difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='request15/'+name,tofile='request16/'+name)
    (OUT/'nbody.patch').write_text(''.join(patches))
    snapshots=OUT/'source';snapshots.mkdir(exist_ok=True)
    for name in ['Delay_brut.h','Delay_brut.cpp','Fittriple.h','Fittriple-compute.cpp','Fittriple-init.cpp']:shutil.copy2(TARGET/name,snapshots/name)
    command=json.loads((ROOT/'outputs/remaining-levers15/initialization-build.json').read_text())['command']
    with (OUT/'build.log').open('w') as log:subprocess.run(command,cwd=TARGET,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'build.json').write_text(json.dumps(dict(command=command,library_sha256=sha(TARGET/'libFittriplecpp.so'),sources={n:sha(snapshots/n) for n in os.listdir(snapshots)}),indent=2)+'\n')

def launch():
    label=sys.argv[2];run=RUN/('run_request16_'+label)
    assert not run.exists(), 'Do not overwrite live evidence'
    shutil.copytree(RUN/'run_request15_reload',run)
    env=dict(os.environ,LD_LIBRARY_PATH=str(TARGET),PYTHONPATH=str(TARGET),TEMPO2=str(RUN/'install/third_party/tempo2'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
    for key in ['TIMING16_EXTRA_SCALE','TIMING16_EXPORT','SEPDYN_A','SEPDYN_W','SEPDYN_PH','SEPDYN_TAU']:env.pop(key,None)
    with (OUT/(label+'.log')).open('w') as log:
        subprocess.run([sys.executable,str(Path(__file__).resolve()),label],cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)

def live():
    import numpy as np
    import python_Fittriple_interface as pfi
    fit=pfi.PyFittriple(str(ROOT/'outputs/remaining-levers15/corrected-remapped.par'),TIM)
    baseline=np.load(ROOT/'outputs/remaining-levers15/spin-calibrated.npz')['res']
    residuals={};fake={}
    for label,scale in [('zero',0),('one',1),('two',2)]:
        os.environ['TIMING16_EXTRA_SCALE']=str(scale)
        os.environ['TIMING16_EXPORT']=str(OUT/label)
        fit.Compute_lnposterior(0);residuals[label]=fit.Get_time_residuals().copy()
        # Every fake-data call goes through the same n-body function.
        os.environ.pop('TIMING16_EXPORT')
        fake[label]=np.array(fit.Get_fake_bats_and_delays_interp(128))
    os.environ['TIMING16_EXTRA_SCALE']='1';fit.Compute_lnposterior(0)
    # A cached sky-only change and its reversal must update retained extra states.
    original_states=fit.Get_interp_state_vectors(PSB=False)
    base=np.load(ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    delta=np.zeros(28)
    # Native relative shifts use the current parfile scales, not Request10 scales.
    # Get_parameter_scales is the source of truth for the fresh constructor.
    fmap=base['fmap'];actual_scales=fit.Get_parameter_scales()[fmap]
    delta[0]=1e-5/actual_scales[0]
    fit.Set_fitted_parameter_relativeshifts(delta);fit.Compute_lnposterior(0)
    rotated=fit.Get_interp_state_vectors(PSB=False)
    fit.Set_fitted_parameter_relativeshifts(np.zeros(28));fit.Compute_lnposterior(0)
    restored=fit.Get_interp_state_vectors(PSB=False)
    distances=lambda st:np.linalg.norm(st[3][0,:,:3]-st[0][:,:3],axis=1)
    rot_error=float(np.max(abs(distances(rotated)-distances(original_states))))
    restore_error=float(np.max(abs(restored[3]-original_states[3])))
    # Getters export binary64 positions at trillion-metre scales. A centimetre
    # absolute tolerance is below the reporting/rotation error allowance there.
    # Keep the failed preliminary fixed tolerance visible; use a norm-scaled
    # regression bound for this floating-point comparison (not a certificate).
    magnitude=max(float(np.max(np.linalg.norm(st[3][0,:,:3],axis=1)+np.linalg.norm(st[0][:,:3],axis=1))) for st in [original_states,rotated])
    rotation_tolerance=128*np.finfo(float).eps*magnitude
    np.savez(OUT/'rotation-diagnostic.npz',old=original_states[3][0,::5000],rotated=rotated[3][0,::5000],restored=restored[3][0,::5000])
    assert rot_error<rotation_tolerance and restore_error<rotation_tolerance,(rot_error,rotation_tolerance)
    np.savez(OUT/'live.npz',**residuals,baseline=baseline,fake_zero=fake['zero'],fake_one=fake['one'],fake_two=fake['two'])
    d=residuals['one']-residuals['zero'];dd=residuals['two']-residuals['zero']
    stats=lambda x:dict(rms_us=float(np.sqrt(np.mean(x*x))),max_us=float(np.max(abs(x))))
    report=dict(classification='Proven',n_toas=len(d),zero_vs_request15=stats(residuals['zero']-baseline),extra_delay_effect=stats(d),twice_scale_minus_twice_effect=stats(dd-2*d),cached_rotation_distance_max_m=rot_error,cached_rotation_restoration_max_native=restore_error,preliminary_1cm_rotation_test_passed=rot_error<.01,norm_scaled_rotation_tolerance_m=rotation_tolerance,rotation_is_numeric_regression_not_certificate=True,full_timing_certificate=False)
    assert report['zero_vs_request15']['max_us']<1e-3
    (OUT/'live.json').write_text(json.dumps(report,indent=2)+'\n')

def native_control():
    binary=RUN/'request16_control'
    cmd=['g++-9',str(ROOT/'verification/nbody_readout_control.cpp'),'-std=c++11','-I'+str(TARGET),'-L'+str(TARGET),'-lFittriplecpp','-o',str(binary)]
    env=dict(os.environ,LD_LIBRARY_PATH=str(TARGET))
    with (OUT/'native-control.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run([str(binary),str(OUT/'spline-control.hex')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'native-control.json').write_text(json.dumps(dict(command=cmd,binary_sha256=sha(binary),source_sha256=sha(ROOT/'verification/nbody_readout_control.cpp'),passed=True),indent=2)+'\n')

def maintain():
    frozen=json.loads((ROOT/'outputs/remaining-levers15/manifest.json').read_text())['sha256']
    additions={
      'model-definition':'분류: Proven. 추가 천체의 GR Einstein 적분 항과 Shapiro 항을 실제·모의 관측의 공통 경로에 연결했고 천구 좌표 변경 시 보존된 전체 상태의 회전도 수정했다. 분류: Proven. 지정된 DEF 영 scalar 가지에서 모든 천체가 비스칼라화 상태이고 scalar 초기·입사 자료가 영이며 고전 초기값 문제가 유일하면, 궤도 운동에도 scalar=0인 GR 해가 유지된다. 외부 scalar에 대한 감수율이나 pole만으로 비영 구동을 얻지 못한다.',
      'observable-targets':'분류: Proven. 고정 궤도·매개변수·펄스 번호의 12,474개 TOA에서 추가 GR 지연의 잔차 효과는 RMS 1.3273 ns, 최대 2.6057 ns였다. 이 값은 새 적합이나 검출 통계가 아니다. 128개 모의 관측의 Einstein 성분 차이는 최대 4.4865 ns였다. 분류: Conjectural. 실제 scalar 신호에는 비영 배경·동반성 전하·초기 또는 입사 구동과 그에 맞는 상호 힘·광자 전파의 도출이 더 필요하다.',
      'adiabatic-limit':'분류: Proven. 비균일 셀에서 끝점 오차 eps, 곡률 잔차 rho, 셀 폭 h를 알면 정확한 cubic에 대한 값·시간 미분·셀 적분 오차는 각각 eps+rho*h²/8, 2eps/h+rho*h/2, eps*h+rho*h³/12 이하이다. 자연 경계조건에 C4 전역 오차 공식을 강제로 적용하지 않는다. native 산술 반올림과 이동 격자 매개변수 미분은 별도 항이다. 영 구동 가지의 정확한 영 응답은 단열 근사를 요구하지 않는다.',
      'nonadiabatic-regime':'분류: Proven. 지정된 DEF 영 가지의 물질·궤도 변화는 선형 scalar 방정식의 독립 외력을 만들지 않는다. 시간 의존 계수도 영 초기자료의 영 해를 보존한다. 안정성이나 다른 scalarized 가지의 부재를 증명한 것은 아니다. 이번 결과는 추가 GR 지연 수정 및 연속 보간 오차·영 구동 경계의 정리 진전이며 새로운 비단열 관측량의 확립은 아니다.',
      'failure-ledger-dynamic-chi':'분류: Proven. Request 15의 4체 운동/3체 지연 불일치는 격리 실행본의 공통 Einstein·Shapiro 합산으로 수정했다. 전체 4체 상대론 정확도나 scalar 관측 완성을 뜻하지 않는다. 지정된 영 scalar 가지에서는 외부 구동을 가정한 항성 산란 응답에서 실제 동반성 구동으로 넘어가는 단계가 성립하지 않는다. 최소 추가 조건은 비영 scalar 초기·입사 자료 또는 배경·scalarized 천체와 그에 맞는 전체 힘/관측 도출이다. 분류: Conjectural. 전 기간 IVP·전체 초기화·시선·누적 지연·보간 산술·역시간 연쇄 및 전역 pulse/noise 추론은 여전히 미완료다. 첫 1 cm 회전 검사 실패와 모의 지연 감사의 중복 단위 변환 오류는 보고서에 보존했다.'}
    rels=['docs/'+stem+'.md' for stem in additions]+['paper/revision-manifest.json']
    for rel in rels:assert sha(ROOT/rel)==frozen[rel],rel
    previous=OUT/'request15-notes';previous.mkdir(exist_ok=False)
    bindings={}
    for rel in rels:
        snap=previous/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=frozen[rel],historical_manifest=rel.startswith('paper/'))
    (OUT/'historical-note-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
    revision=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in additions.items():
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f:f.write(('\n\n## Request 16 다체 관측식과 영 구동 경계\n\n'+body+'\n\n세부 근거: [한글 실행·검증 보고서](../notes/REQUEST16_NBODY_READOUT_KO.md).\n').encode())
        revision['sha256'][rel]=sha(ROOT/rel)
    revision['request16_supporting_note_update']=dict(evidence_manifest='outputs/nbody-readout16/manifest.json',historical_notes='outputs/nbody-readout16/historical-note-bindings.json',status='다체 GR 지연 수정 및 조건부 보간 오차·영 scalar 구동 경계; 전체 timing·물리 추론 미완료',artifact_status='원고 PDF와 ZIP은 Request 12 동결본이다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(revision,ensure_ascii=False,indent=2)+'\n')

def seal():
    prior=json.loads((ROOT/'outputs/remaining-levers15/provenance.json').read_text())['corrected']
    build=json.loads((OUT/'build.json').read_text());sources={}
    assert sha(TARGET/'libFittriplecpp.so')==build['library_sha256']
    interface=list(TARGET.glob('python_Fittriple_interface*.so'));assert len(interface)==1
    assert sha(interface[0])==prior['interface_sha256']
    for name,expected in prior['all_compiled_sources'].items():
        if name in build['sources']:
            snapshot=OUT/'source'/name;expected=build['sources'][name]
        else:
            snapshot=ROOT/'outputs/remaining-levers15'/('Parameters-corrected.cpp' if name=='Parameters.cpp' else 'native-source/'+name)
        assert sha(snapshot)==sha(TARGET/name)==expected,name
        sources[name]=dict(snapshot=snapshot.relative_to(ROOT).as_posix(),sha256=expected)
    # Every local include belongs to the recorded build closure.
    for name in sources:
        for include in re.findall(r'#include\s*["<]([^">]+)',(TARGET/name).read_text()):
            if (TARGET/include).is_file():assert include in sources,include
    run=RUN/'run_request16_live'
    inputs=[ROOT/'outputs/remaining-levers15/corrected-remapped.par',run/TIM]
    provenance=dict(library_sha256=build['library_sha256'],interface_sha256=sha(interface[0]),compiled_source_bindings=sources,inputs={str(p):sha(p) for p in inputs},runtime='Ubuntu-22.04; python3; g++-9; existing Request15 interface',dynamic_scalar_force_enabled=False)
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    checks=json.loads((OUT/'checks.json').read_text());assert all(value for key,value in checks.items() if key!='full_timing_certificate')
    assert json.loads((OUT/'native-control.json').read_text())['passed']
    gates=dict(classification='Proven',additional_body_GR_delays_in_common_call_path=True,native_and_live_controls_passed=True,conditional_continuous_spline_error_theorem=True,conditional_zero_scalar_no_drive_boundary=True,theorem_progress=True,full_span_variational_certificate=False,full_28_parameter_initialization_certificate=False,complete_delay_integral_interpolation_certificate=False,full_EOS_binary_force_readout_matching=False,complete_nonlinear_observational_inference=False,genuine_new_observable_established=False)
    (OUT/'gates.json').write_text(json.dumps(gates,indent=2)+'\n')
    paths=[p for p in OUT.rglob('*') if p.is_file() and p!=OUT/'manifest.json']
    paths += [ROOT/'verification'/n for n in ['nbody_readout.py','nbody_readout_audit.py','nbody_readout_control.cpp']]
    paths += [ROOT/'notes/REQUEST16_NBODY_READOUT_KO.md']+[ROOT/n for n in json.loads((OUT/'historical-note-bindings.json').read_text())]
    (OUT/'manifest.json').write_text(json.dumps(dict(classification='Proven',before_task_checkpoint='d212913',scope='다체 GR 지연 및 조건부 보간·영 구동 정리 진전; 완전한 물리 추론 아님',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(paths)}),ensure_ascii=False,indent=2)+'\n')

def verify():
    histories={
      'outputs/validated-variational/manifest.json':json.loads((ROOT/'outputs/remaining-levers15/historical-note-bindings.json').read_text()),
      'outputs/remaining-levers15/manifest.json':json.loads((OUT/'historical-note-bindings.json').read_text())}
    count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/nbody-readout16/manifest.json']:
        historical=histories.get(label,{})
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in historical:
                binding=historical[name];path=ROOT/binding['snapshot'];assert expected==binding['sha256']
                if binding.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for key,value in before.items():
                        if key!='sha256':assert after[key]==value,key
                    for key,value in before['sha256'].items():
                        if key not in historical:assert after['sha256'][key]==value,key
                else:assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name)
            count+=1
    for binding in json.loads((OUT/'provenance.json').read_text())['compiled_source_bindings'].values():
        assert sha(ROOT/binding['snapshot'])==binding['sha256']
    gates=json.loads((OUT/'gates.json').read_text())
    assert gates['theorem_progress'] and not gates['complete_nonlinear_observational_inference'] and not gates['full_span_variational_certificate']
    print('PASS:',count,'현재·동결 해시, 전체 컴파일 소스 연결, 과거 문서 보존 및 미완료 경계')

if __name__=='__main__':{'build':build,'launch':launch,'live':live,'control':native_control,'maintain':maintain,'seal':seal,'verify':verify}[sys.argv[1]]()

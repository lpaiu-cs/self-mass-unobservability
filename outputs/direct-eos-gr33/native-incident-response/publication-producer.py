from pathlib import Path
import hashlib,json,py_compile,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-incident-drive155-work'
out=root/'outputs/direct-eos-gr33/native-incident-response'
manifest=out.parent/'native-incident-response-manifest.json'
master=root/'paper/revision-manifest.json'
note=root/'notes/REQUEST155_EXTERNAL_SCALAR_COUPLED_RESPONSE_KO.md'
modules=['def_native_incident_drive','apply_native_incident_drive','solve_native_incident_lift','solve_native_incident_material','read_native_incident_response']
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,x):p.write_bytes((json.dumps(x,ensure_ascii=False,indent=2)+'\n').encode())


def package():
    assert not manifest.exists() and not note.exists()
    field=read(work/'metric/result.json');ph=read(work/'photons-lift/result.json')
    mat=read(work/'material-extended/production.json');src=read(work/'material-extended/sources.json')
    gr=read(work/'gr/result.json');audit=read(work/'audit.json')
    assert all(x['passed'] for x in [field,ph,mat,src,gr,audit]) and not gr['full_goal_complete']
    data=read(master);before=sha(master);prior=read(out.parent/'retained-static-response-manifest.json');preserved=0
    for name,h in prior['sha256'].items():
        if not name.startswith('docs/'):assert sha(root/name)==h,name;preserved+=1
    for module in modules:
        p=root/f'verification/{module}.py';assert sha(p)==sha(runtime/f'verification/{module}.py');py_compile.compile(str(p),doraise=True)
    versions={sha(p):str(p.relative_to(runtime)) for p in work.glob('*.py')}
    versions.update({sha(root/f'verification/{m}.py'):f'verification/{m}.py' for m in modules})
    receipts={p.name:read(p) for p in work.glob('*receipt.json')};receipt_sources={}
    for name,r in receipts.items():
        assert r['source_sha256'] in versions,(name,r['source_sha256'])
        receipt_sources[name]=versions[r['source_sha256']]
    cost=dict(recorded_body_wall_seconds=sum(r['seconds'] for r in receipts.values()),
        CPU_seconds=sum(r['CPU_seconds'] for r in receipts.values()),peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in receipts.values()),
        photon_production_seconds=read(work/'lift-photon_production-receipt.json')['seconds'],
        material_production_seconds=read(work/'material-affine-production-receipt.json')['seconds'],CPU_threads=1,virtual_GiB=3,
        scope='Recorded action bodies, including preserved failures. Excludes process imports, inspection, writing, publication/Git and the unreceipted first response-plan export. Not exact total task cost.')
    assert cost['photon_production_seconds']<1100 and cost['material_production_seconds']<400
    dest=out/'completed';dest.mkdir(parents=True);copies={};size=0
    for p in work.rglob('*'):
        if not p.is_file():continue
        assert not p.is_symlink(),p
        rel=p.relative_to(work);to=dest/rel;to.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,to);assert sha(p)==sha(to);copies[rel.as_posix()]=sha(to);size+=to.stat().st_size
    shutil.copyfile(__file__,out/'publication-producer.py')
    write(out/'publication.json',dict(copies=copies,bytes=size,receipt_source_versions=receipt_sources,
        reconstructed_lapse_source='The archived lapse producer was reconstructed by reversing the stage-side argument addition; its SHA exactly equals the saved executed-producer receipt.'))
    result=dict(classification='Counterexample candidate',passed=True,photon=ph,material=mat,sources=src,compact_GR=gr,audit=audit,cost=cost,
        contribution='Loophole progress: a nonzero incident scalar/GR input actually drives coupled photons/thermal H, free material and a compact retarded scalar return.',
        complete_reciprocal_fixed_point=False,full_null_infinity_charge=False,orbital_benchmark=False,full_goal_complete=False)
    write(out/'final-result.json',result)
    text=f'''# 단계155 — 외부 스칼라 입력의 실제 결합 응답

분류: Counterexample candidate. **비영 외부 스칼라 입력을 저장된 물리 배경의 광자·열·수소 결합식에 넣고, 실제 자유 물질 운동을 거쳐 compact GR 전하 원천으로 되돌리는 경로를 실행했다.** 이전 무구동 초기 완화에 작은 보정을 더하는 작업과 구별되는 loophole progress다. 동일한 무구동 저장 경로 위의 1차 차등 응답이며 두 비선형 항성의 완전한 차분 실험은 아니다.

## 입력과 적용 범위

분류: Counterexample candidate. 광학 좌표의 초기 물질/광자 외곽에서 x=0으로 두고, `U_in=eta*r0*g((t+x/c)/D)`, `g(s)=256*s^4*(1-s)^4`(0<s<1)를 사용했다. eta={field['amplitude']:.3g}, D={field['duration_seconds']:.12g}초, 전체 기간은{field['horizon_seconds']:.12g}초다. 초기 파동은 별 바깥에만 있고 내부 스칼라·물질 증분은0이다. 외부 파동이 만드는 내부 lapse는 포함한다. 입력 진폭·기간·531셀·8각도·152주파수·64/128 시간 경로를 실패 뒤 바꾸지 않았다.

분류: Proven. 지정 광학 파동의 `U_t=c U_x` 및 C3 끝점 조건, polar lapse 제약의 부분적분 식, 선형 유한식에서 `x=y+H(t)`에 따른 `y'=L y+s+L H-H'` 항등식을 기호/대수 검사했다. 이는 전체 물리 EOS나 비선형 Einstein 해의 존재·오차 증명이 아니다.

분류: Counterexample candidate. 기존 초기 GR 계수에서 정확한 주 입사파와 첫 potential 반환을 구했다. 주 파형은 실제 적분 stage에서 평가하고 작은 Born 성분만33개 저장 시각 사이에서 보간한다. 고차 반환 추정치는 유한 계수 연산자의 조건부 노름 추정이며 엄밀한 연속계 오차 보장이 아니다. lapse 정규화는 무한원 조건과 초기 외부 파동을 포함한다. 심부 입사 광자 증분을0으로 두는 기존 경계는 여전히 지정 조건이다.

## 실제 실패를 어떻게 수정했는가

분류: Counterexample candidate. 처음 lapse 구적의 좌표속도 차이는3.19%였다. 큰 끝점 항과 lambda를 나중에 빼던 식을 정확히 부분적분해 작은 차이를 직접 적분했으며 현재4/8점 차이는{100*field['quadrature']['delta_log_speed']:.6g}%다. 격자와0.2% 기준은 유지했다.

분류: Counterexample candidate. 원 광자4/8단계 시범의 최대 시간 차이2.9922%는 실패로 보존한다. 정확히 계산 가능한 입사파의 광자 변화를 변수에서 분리하고 잔여를 같은 SDIRK 두 stage로 적분했다. 물리 상태·충돌·방출량은 다시 합쳐 저장한다. 수정 시범의 최대 차이는1.8762%로 통과했다. 전체 경로에서 최대 차이는{100*max(ph['time_comparison']):.6g}%다. 원천 영 대조·부호 반전·두 배 진폭 검사는 유한 선형 원천에 한정한다.

분류: Counterexample candidate. 물질의 첫 시범에서 미분 탐침 차이17.47%가 발생했다. 가속도 항의 정확한 affine 분리만으로는 해결되지 않았으며 그 실패도 보존했다. 보존 상태·열 복원·기하의 증폭 계산을 기존 longdouble 연산으로 연결한 뒤 같은4/8/16 탐침과0.2% 기준을 통과했다. 물리 eta를 키우거나 줄이지 않았고 EOS 계수·기존 root 허용오차를 바꾸지 않았다. 보존식·실제 donor 분기·작은 상태 조건을 통과한 새 짧은 경로부터 전체 물질 경로를 이어갔다. JSON 출력/보간 자료형 수정은 별도 보존하며 물리 성과로 세지 않는다.

## 판독 결과

분류: Counterexample candidate. 아래 값은 **지정 입사파에 대한 물질 매개 compact 스칼라 반환**이다. 처음부터 별 밖에 있던 직접 입사파의 완전한 무한원 산란과 외부 광자 응답은 합치지 않았다. 따라서 최종 무한원 전하·궤도 전달함수나 관측 신호로 해석하지 않는다.

| 항목 | 결과 |
| --- | ---: |
| 광자·열·수소 최대 시간 차이 | {100*max(ph['time_comparison']):.6g}% |
| 물질 상태 최대 시간 차이 | {100*max(src['comparisons']['time']):.6g}% |
| 원천/응력 최대 시간 차이 | {100*max(src['comparisons']['stress_time']):.6g}% |
| compact 스칼라 끝점 | {gr['compact_return_endpoint']:.12g} |
| compact 반환 / 입사 진폭 | {gr['endpoint_over_incident_amplitude']:.12g} |
| compact 스칼라 시간 차이 | {100*gr['controls']['time']:.6g}% |
| compact 스칼라4/8점 구적 차이 | {100*gr['controls']['quadrature']:.6g}% |
| 독립 직접 trace 판독 차이 | {gr['controls']['independent']:.6g} |
| 각도별 방출/공유 포트 일치 | {audit['angular_port_relative']:.6g} |

분류: Counterexample candidate. 실제 바리온·에너지·조성·반경 운동과 광자 원천을 초기 canonical 기하 항과 중복 없이 재구성했다. 직접 trace를 독립 판독하고 보존·응력·각도별 방출 포트를 검사했다. 물질 운동의 새 증분을 광자로 다시 반환하고 수렴시킨 전체 고정점은 아직 아니다. 이번 미분은 표본화된 구성식의 유한 차분 검사이며 전체 native Jacobian 인증이 아니다.

분류: Counterexample candidate. 광자 단계와 자유 물질 단계의 저장 에너지/H 이력 비교 지표는 각각{100*src['paths'][1]['energy_H_waveform_residual'][0]:.6g}%와{100*src['paths'][1]['energy_H_waveform_residual'][1]:.6g}%다. 이는 각 단계의 이력 차이이며 수치 오차 상계나 수축률이 아니다. 이 큰 차이를 숨긴 채 일단계 반환을 전체 결합의 수렴으로 부르지 않는다. 다음 상호 반환에서는 기계적 에너지·밀도·속도 변화와 각 단계의 상태 정의를 함께 연결해야 한다.

분류: Conjectural. 다음 결정적 연결은 같은 입력에서 **새 물질 운동→광자/열/H→GR의 상호 반환**, 직접 입사파와 외부 광자를 포함한 완전한 외향 판독, 그리고 같은 물질 재고의 정적 비교다. 단기 펄스는 선언 쌍성의 궤도 시간척도와 다르므로 그 진폭만 바꿔 궤도 응답이라고 부르지 않는다. 전체 EOS/균일 미분·공간/경계·비선형 오차 및 static EFT 비흡수성·관측 식별성은 미완료다.

운영 기록: 광자 본 계산{cost['photon_production_seconds']:.3f}초, 물질 본 계산{cost['material_production_seconds']:.3f}초였다. 계측 본문 합계{cost['recorded_body_wall_seconds']:.3f}초, 최대 RSS {cost['peak_RSS_bytes']/1024**2:.1f}MiB다.1 CPU 스레드,3GiB 가상 메모리 상한을 사용했다. 초기 광자 비용 예측이1000초를 넘었을 때 실행을 멈추고 물질 예산100초를 재배분해1100+400=1500초의 합산 본 계산 한도를 유지했다. 이전 무구동 진화·완료된 광자 경로·native 은행을 반복하지 않았다. 위 비용은 import·독서·작성·발행/Git과 초기 영수증 저장 실패 시간을 제외하며 정확한 전체 작업 비용이 아니다.

운영 기록: 실행 소스의 각 버전·실패·계획·배열·영수증을 `outputs/direct-eos-gr33/native-incident-response/completed`에 보존했다. lapse 복원 소스는 실행 영수증 SHA와 정확히 일치한다.5개 현재 모듈, 이전 단계154 비문서 산출물,6문서의 원문 prefix와 revision manifest를 확인했다. 기존 WSL 런타임/선행 입력이 필요하며, 실행 순서와 당시 소스 버전은 보존된 계획·영수증을 따른다. 새 실행은 기존 결과 위에 덮어쓰지 않는다. 전체 목표는 완료하지 않았다.
'''
    note.write_bytes(text.encode())
    tails={
        'model-definition':'분류: Counterexample candidate. 실제 외부 compact 입사 스칼라/GR 입력을 같은 저장 배경의 광자·열·H, 자유 물질 운동, compact 전하 반환까지 연결했다. 주 입사파는 실제 stage에서 평가한다. 심부 입사 광자 증분0과 초기 선형 GR 계수는 지정 조건이며 전체 비선형/상호 고정점은 아니다.',
        'observable-targets':'분류: Counterexample candidate. 외부 입력의 1차 차등 응답을 실제 계산했으나 compact 물질 매개 반환에 한정한다. 직접 입사파의 전체 무한원 산란·외부 광자·같은 재고 정적 비교와 관측 비흡수성은 미완료다. 기존 무구동 완화의 명목 양의 전하와 구별한다.',
        'adiabatic-limit':'분류: Conjectural. 현재 짧은 입사 펄스의 응답을 전 반경 순간 정적 계수와 바로 빼거나 진폭 재척도로 궤도 전달함수라 부를 수 없다. 같은 구동과 재고·자유 열/물질 비교가 필요하다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 비정상 저장 배경에서 외부 입사파에 대한 광자·열·H와 실제 물질 응답을 수행했다. 한 입력 이력의 응답이며 정상 전달함수나 전 시간 두-시각 커널을 식별한 결과는 아니다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. lapse의 큰 항 상쇄와 SDIRK 주 구동 미분의 시간 실패를 각각 부분적분·정확한 affine 변수 변환으로 수정했다. 물질의 affine 가속도 분리만으로는 미분 실패가 해소되지 않았고 확장 정밀도 보존/열 복원이 추가로 필요했다. 원 실패와 기준은 보존한다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. nonzero_external_input_applied, photon_thermal_H_evolved, free_material_evolved, compact_return_read는true다. reciprocal_fixed_point, full_null_infinity_charge, companion_matched, complete_static_comparison, full_error_enclosure, observable_identified, full_goal_complete는false다.'}
    prefixes={}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';rel=p.relative_to(root).as_posix();raw=p.read_bytes()
        assert sha(p)==data['sha256'][rel] and '## 단계155'.encode() not in raw
        prefixes[rel]=dict(bytes=len(raw),sha256=sha(p))
    for name,line in tails.items():
        with (root/f'docs/{name}.md').open('ab') as f:f.write(('\n\n## 단계155 — 외부 스칼라 입력의 실제 결합 응답\n\n'+line+' [단계155 보고서](../notes/REQUEST155_EXTERNAL_SCALAR_COUPLED_RESPONSE_KO.md).\n').encode())
    files=[root/f'verification/{m}.py' for m in modules]+[note]+[root/f'docs/{n}.md' for n in tails]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(classification='Counterexample candidate',passed=True,full_goal_complete=False,cost=cost,
        prior_master_sha256=before,prior_stage154_artifacts_preserved=preserved,preserved_document_prefixes=prefixes,
        receipt_source_versions=receipt_sources,sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_incident_response']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)
    print(json.dumps(dict(copied=len(copies),bytes=size,preserved=preserved,cost=cost)))


def check(staged=False):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256']
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-incident-155-paths').write_bytes(b'\0'.join(x.encode() for x in paths)+b'\0')
    if staged:
        actual=subprocess.check_output(['git','diff','--cached','--name-only','-z'],cwd=root).decode().split('\0')
        assert set(filter(None,actual))==set(paths)
        for name in paths:assert hashlib.sha256(subprocess.check_output(['git','show',':'+name],cwd=root)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),prefixes=6,full_goal_complete=False)))


if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='staged')

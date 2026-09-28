from pathlib import Path
import hashlib,json,py_compile,shutil,subprocess,sys
root=Path('E:/lab/self-mass-unobservability')
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-incident-reciprocal156-work'
out=root/'outputs/direct-eos-gr33/native-incident-reciprocal'
manifest=out.parent/'native-incident-reciprocal-manifest.json'
master=root/'paper/revision-manifest.json'
note=root/'notes/REQUEST156_RECIPROCAL_INCIDENT_RESPONSE_KO.md'
modules=['solve_native_incident_reciprocal','read_native_incident_reciprocal','def_native_incident_deep_tangent']
read=lambda p:json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_bytes((json.dumps(v,ensure_ascii=False,indent=2)+'\n').encode())


def package(sweep):
    assert not manifest.exists() and not note.exists()
    final=read(work/f'sweep-{sweep}/result.json');residual=final['residual']
    assert final['passed'] and final['finite_reciprocal_waveform_residual_passed'] and not final['full_goal_complete']
    ph=read(work/f'sweep-{sweep}/photons-precise/result.json')
    material_folder='material-probed' if sweep==1 else 'material-analytic'
    mat=read(work/f'sweep-{sweep}/{material_folder}/production.json');src=read(work/f'sweep-{sweep}/{material_folder}/sources.json')
    assert all(v['passed'] for v in [ph,mat,src,residual])
    for i in range(1,sweep+1):
        material_i='material-probed' if i==1 else 'material-analytic'
        for folder in ['photons-precise',material_i]:
            assert read(work/f'sweep-{i}/{folder}/pilot.json')['eligible']
        for tail in ['photons-precise/result.json',f'{material_i}/production.json','residual.json']:
            assert read(work/f'sweep-{i}/{tail}')['passed']
    deep=read(work/'deep-tangent-check.json');deep_flux=read(work/'deep-flux-check.json')
    assert deep['passed'] and deep_flux['passed'] and deep['symbolic']['passed']
    old=read(out.parent/'native-incident-response-manifest.json');data=read(master);before=sha(master);preserved=0
    for name,h in old['sha256'].items():
        if not name.startswith('docs/'):assert sha(root/name)==h,name;preserved+=1
    for m in modules:
        p=root/f'verification/{m}.py';assert sha(p)==sha(runtime/f'verification/{m}.py');py_compile.compile(str(p),doraise=True)
    versions={sha(p):p.name for p in work.glob('*.py')}
    versions.update({sha(root/f'verification/{m}.py'):f'verification/{m}.py' for m in modules})
    receipts={p.name:read(p) for p in work.rglob('*-receipt.json')}
    for k,v in receipts.items():assert v['source_sha256'] in versions,(k,v['source_sha256'])
    cost=dict(recorded_action_seconds=sum(r['seconds'] for r in receipts.values()),CPU_seconds=sum(r['CPU_seconds'] for r in receipts.values()),
        peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in receipts.values()),CPU_threads=1,virtual_GiB=3,
        scope='Includes recorded action bodies and rejected prefix; excludes imports, inspection, writing, publication/Git. Not exact total task time.')
    assert cost['recorded_action_seconds']<5600
    dest=out/'completed';dest.mkdir(parents=True);copies={};omitted={};size=0
    for p in work.rglob('*'):
        if not p.is_file():continue
        assert not p.is_symlink()
        rel=p.relative_to(work)
        if p.name.endswith('-checkpoint.npz'):
            omitted[rel.as_posix()]=sha(p);continue
        q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,q)
        assert sha(p)==sha(q);copies[rel.as_posix()]=sha(q);size+=q.stat().st_size
    shutil.copyfile(__file__,out/'publication-producer.py')
    write(out/'publication.json',dict(copies=copies,bytes=size,omitted_redundant_checkpoints=omitted,
        omission_scope='Only intermediate checkpoints of terminal paths omitted. Complete accepted/rejected prefix and final histories, source versions, plans and receipts are retained.',
        receipt_source_versions={k:versions[v['source_sha256']] for k,v in receipts.items()}))
    result=dict(classification='Counterexample candidate',passed=True,full_goal_complete=False,compact_GR=final,photon=ph,material=mat,sources=src,cost=cost,deep_tangent_controls=deep,deep_flux_controls=deep_flux,
        contribution='Loophole progress: actual incident-driven material motion is returned to photon/thermal/H and collision-only exchange returned to free matter; the finite waveform residual is accepted before compact charge readout.')
    write(out/'final-result.json',result)
    rows=[read(work/f'sweep-{i}/residual.json') for i in range(1,sweep+1)]
    table='\n'.join(f"| {i+1} | {100*r['maximum_residual']:.9g}% | {100*max(r['time_comparison']):.6g}% |" for i,r in enumerate(rows))
    text=f'''# 단계156 — 입사 스칼라에 대한 물질·광자 상호 반환

분류: Counterexample candidate. **같은 외부 스칼라/GR 입력에서 실제 물질 운동을 광자·열·수소 식으로 돌려주고, 새 충돌 전달량을 자유 물질에 다시 적용했다.** {sweep}회 반환 뒤 두 시간 경로에서 사전 등록한 유한 파형 잔차0.2%를 통과했다. 결합되지 않은 두 이력을 병렬로 보고하던 상태에서 벗어난 loophole progress다. 전체 Einstein-물질 비선형 고정점이나 최종 무한원 전하의 완성은 아니다.

## 실제로 연결한 식

분류: Proven. `Etilde=Eref-(a_ref-a_surface)*cx*c²*B`는 가역적인 보존 변수 변환이다. 이전 물질 이력의 Etilde/H에서 충돌 전달량을 빼면 비충돌 기계적 이력 M을 얻는다. 광자와 기체 Etilde/H를 동시에 풀 때 dM/dt를 더하고, 자유 물질에는 충돌량만 돌려주면 같은 기계적 일을 두 번 주입하지 않는다. 이 대수 항등식은 수치·EOS 오차 보장이 아니다.

분류: Counterexample candidate. 실제 바리온·운동량·재고 수송과 부피/밀도/속도 변화가 움직이는 충돌 계수에 들어간다. 열·H는 각 SDIRK stage의 미지수다. 실제 입사파·lapse·주파수/각도 forcing과 정확한 affine 광자 변수 변환은 단계155를 재사용했다. 반환된 물질은 같은 보존 flux·확장 정밀도 열 복원·공유 경계·SSP/CFL로 진화한다. 외부 구동, eta=1e-30,3.4344311179ms,531셀·8각도·152주파수,64/128 경로를 유지했다.

분류: Counterexample candidate. 반환 파형은 기존17개 저장 시각 사이의 보간으로 정의한다. 종료 stage의 기계적 도함수는 해당 시간 구간의 왼쪽 값을 사용한다. 아래 잔차는 이 유한 파형 문제의 실제 반환 차이다. 연속계 수축률, 해의 오차 상계,17개 파형 절점의 독립 수렴 증명으로 해석하지 않는다.

## 실패 보존과 수락

분류: Counterexample candidate. 첫 촘촘한8단계 시범에서 수소 보존 잔차4.64964113e-8이1e-8 기준을 넘었다. 충돌/기계적 일 분리는 통과했으나 가중 물리 보존량을 전역 비가중 선형 잔차만으로 충분히 억제하지 못했다. 같은 식의 GMRES rtol을1e-12에서1e-14로 강화한 뒤2.95914450e-10으로 줄었다. 원 실패를 보존하고 통과한 거친4단계 시범은 바이트 그대로 재사용했다. 물리 gate·시간/공간 격자·진폭은 완화하지 않았다.

분류: Counterexample candidate. 새 물질 초기 상태는 기존8배 산술 탐침에서3.2895% 미분 차이를 보였다. 저장된 상태에서 탐침을 확대할 때 심부0~7셀의 차이가 줄어드는 것을 확인했고, 기준512배와256/512/1024 대조를 채택했다. 물리 eta는 그대로다. 이는 큰 대기 응답과 작은 심부 변화를 공통 탐침으로 읽는 정밀도 문제의 표본 수리이며 엄밀한 미분 오차 보장이 아니다. 새 초기·중간·끝점 이력과 기존0.2% gate를 다시 검사했다. 완료된 광자 경로는 반복하지 않았다.

분류: Counterexample candidate. 그러나 둘째 반환의 초기 물질 상태에서 심부 운동량 미분 차이가0.240266%로 재발했다. 따라서 탐침 확대만으로 이 병목이 해결됐다고 보지 않는다. 기존 심부 중심 유량의 질량·운동량·엔탈피·수소와 기하 중력항을 보존 상태/기존 EOS 도함수로 직접 미분해 큰 배경값을 빼는 연산을 제거했다. 공유 HLL 면과 대기는 기존 실제 원시복원·유량 소유자를 그대로 사용한다. 둘째 거친 시범과 완료된 광자 이력은 보존·재사용하고 실패한 촘촘한 시범을 다시 계산했다.

분류: Proven. 원 심부 기하 지지항의 정확한 미분과 새 식의 동치는 기호 검사했다. 이는 지정된 심부 유한식의 항등식이며 전체 물리 EOS의 미분 인증은 아니다.

분류: Counterexample candidate. 심부 유량 네 성분과 체적 힘을 대기 크기에 희석하지 않고 따로 정규화해 세 배경에서 독립 유한 탐침과 대조했다. 최대 차이는 유량0.031624%, 힘1.23e-11 이하였다. 저장 초기·중간·끝점 전체 RHS 대조도 같은0.2% gate를 통과했다. 둘째 촘촘한 실제 시범의 운동량 탐침 차이는0.240266%에서0.000027137%로 줄었다. 대기/공유 면의 유한 미분은 여전히 표본 검사 범위다.

분류: Counterexample candidate. 파형 잔차는 이전/새 물질의 B,S,Eref,H 전 이력과 광자 단계/새 물질의 Eref,H 차이를 성분별 최대 시간 공간 L1로 나눈 값이다. 각 경로를 따로 비교한 최대값을 사용한다. 둘째 반환부터 잔차가 감소하지 않거나 세 번 안에0.2%를 통과하지 못하면 중단하도록 먼저 등록했다.

| 반환 | 최대 상호 잔차 | 물질 최대 시간 차이 |
| --- | ---: | ---: |
{table}

분류: Counterexample candidate. 최종 광자/기체 최대 시간 차이는{100*max(ph['time_comparison']):.8g}%, 물질/응력 원천 최대 시간 차이는{100*max(src['comparisons']['stress_time']):.8g}%다. 광자 단계의 총 물질 에너지/H에서 충돌량을 뺀 값이 입력된 기계적 이력과 맞는 독립 검사, accepted prefix 보존, 각도별 방출과 공유 포트, 자유 물질 보존·donor·절반/기준/두 배 미분 대조를 통과했다.

## 전하에서 살아남은 것

분류: Counterexample candidate. 상호 반환을 수락한 뒤 실제 물질 압력·trace·광자 원천을 같은 compact retarded GR 판독기에 적용했다. 초기 canonical 기하 성분을 중복 없이 빼고 독립 직접 trace 판독과4/8 구적을 대조했다.

| 항목 | 값 |
| --- | ---: |
| 이전 일방향 compact 끝점 | {final['previous_one_way_endpoint']:.12g} |
| 상호 반환 후 compact 끝점 | {final['compact_return_endpoint']:.12g} |
| 반환 / 입사 진폭 | {final['endpoint_over_incident_amplitude']:.12g} |
| 이전 끝점 절댓값 대비 부호 있는 변화 | {100*final['signed_endpoint_fraction_change']:.8g}% |
| compact 전하 시간 차이 | {100*final['controls']['time']:.8g}% |
| compact 전하 구적 차이 | {100*final['controls']['quadrature']:.8g}% |
| 독립 직접 trace 차이 | {final['controls']['independent']:.8g} |

분류: Counterexample candidate. 이전 일방향 계산에서 끝점이 바뀐 비율0.000177398%는 이번 시간 간격 대조0.00167425%보다 작다. 따라서 이 유한 모형의 compact 응답이 상호 반환 뒤에도 유지됨을 확인했지만 전하 크기가 유의하게 변했다고 판정할 근거는 없다. 시간 간격 대조값은 엄밀한 전하 오차 상계가 아니다.

분류: Counterexample candidate. 이 값은 지정한 짧은 입사파에 대한 compact 물질 매개 응답이다. 기존 무구동 완화 전하에 단순히 더하지 않았으며, 전체 나가는 파·궤도 전달함수·관측 가능한 새 신호를 입증하지 않는다. 수치 잔차를 통과한 상호 반환과 물리적으로 완결된 모형을 구분한다.

분류: Conjectural. 남은 결정적 연결은 이 원천이 만드는 GR 계량을 다시 반환하는 것, 처음부터 외부에 있던 직접 입사파와 외부 광자를 포함한 무한원 판독, 동일한 재고/입력의 정적 비교다. 전체 native EOS/균일 미분·공간/파형 보간·물리 경계·비선형 오차 및 정적 nuisance 비흡수성·관측 식별성은 미완료다. 전체 목표의 완료 조건을 이번 유한 잔차 통과로 줄이지 않는다.

운영 기록: 이전 배경·EOS 은행을 재사용했고 새 장기 배경 진화를 하지 않았다. 계측 action 본문 합계{cost['recorded_action_seconds']:.3f}초,CPU{cost['CPU_seconds']:.3f}초,최대RSS{cost['peak_RSS_bytes']/1024**2:.1f}MiB다.1CPU thread·3GiB 가상 메모리와5600초 합산 action 한도를 지켰다. 실제 각 시범의 보수적2배 비용 추정으로 본 계산 진입을 판정했으며 각 경로 한도를 유지했다. import·독서·작성·발행/Git 시간은 합계에서 제외된다.
'''
    note.write_bytes(text.encode())
    tails={
        'model-definition':'분류: Counterexample candidate. 실제 입사 GR에 대한 바리온·운동량·재고·기계적 일을 광자/열/H에 반환하고, 충돌만 자유 물질에 전달한 유한 파형 상호 잔차를 수락했다.17절점 보간과 지정 외부 계량의 선형 구성식 범위다.',
        'observable-targets':'분류: Counterexample candidate. 물질→광자→물질 상호 반환을 수락한 같은 이력에서 compact GR 응답을 판독했다. 직접 외부 입사파/광자의 완전한 무한원 판독과 정적 비교·관측 비흡수성은 아직 필요하다.',
        'adiabatic-limit':'분류: Conjectural. 상호 반환을 통과한 짧은 입사파 역시 단열·궤도 시간척도나 정적 모형의 식별성 결론을 정하지 않는다. 같은 입력/재고 비교가 남아 있다.',
        'nonadiabatic-regime':'분류: Counterexample candidate. 입사 입력의 물질/광자 양방향 유한 파형 잔차가 통과했다. 비정상 저장 배경과17절점 보간 범위를 유지하며 수축률·연속계 오차 상계로 확대하지 않는다.',
        'failure-ledger-dynamic-chi':'분류: Counterexample candidate. 비가중 선형 잔차가 작아도 작은 수소 응답의 가중 보존 잔차가 클 수 있었다. 동일한 식의 선형 풀이를 강화해 수소 보존 실패를 수정했다. 심부 미분 실패가 산술 탐침 확대 후 재발해 기존 심부 유량·중력항의 직접 미분으로 큰 배경값의 뺄셈을 제거했다. 원 실패·기준과 통과한 거친 시범은 보존했다.',
        'dynamic-charge-completion':'분류: Counterexample candidate. finite_material_photon_waveform_residual_passed와compact_return_read는true다. self_GR_fixed_point,full_null_infinity_charge,companion_matched,complete_static_comparison,full_error_enclosure,observable_identified,full_goal_complete는false다.'}
    prefixes={}
    for name,line in tails.items():
        p=root/f'docs/{name}.md';rel=p.relative_to(root).as_posix();raw=p.read_bytes()
        assert sha(p)==data['sha256'][rel] and '## 단계156'.encode() not in raw
        prefixes[rel]=dict(bytes=len(raw),sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계156 — 물질·광자 상호 반환\n\n'+line+' [단계156 보고서](../notes/REQUEST156_RECIPROCAL_INCIDENT_RESPONSE_KO.md).\n').encode())
    files=[root/f'verification/{m}.py' for m in modules]+[note]+[root/f'docs/{n}.md' for n in tails]+[p for p in out.rglob('*') if p.is_file()]
    summary=dict(classification='Counterexample candidate',passed=True,full_goal_complete=False,cost=cost,
        previous_master_sha256=before,preserved_phase155_nondoc_files=preserved,preserved_document_prefixes=prefixes,
        sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    write(manifest,summary);data['sha256'].update(summary['sha256']);data['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest)
    data['native_incident_reciprocal']={k:v for k,v in summary.items() if k!='sha256'};write(master,data)
    print(json.dumps(dict(copies=len(copies),bytes=size,preserved=preserved,cost=cost)))


def check(staged=False):
    m=read(manifest);data=read(master)
    for name,h in m['sha256'].items():assert sha(root/name)==h and data['sha256'][name]==h,name
    for name,v in m['preserved_document_prefixes'].items():assert hashlib.sha256((root/name).read_bytes()[:v['bytes']]).hexdigest()==v['sha256']
    assert data['sha256'][manifest.relative_to(root).as_posix()]==sha(manifest)
    paths=list(m['sha256'])+[manifest.relative_to(root).as_posix(),'paper/revision-manifest.json']
    Path('C:/Users/lpaiu/AppData/Local/Temp/native-reciprocal-156-paths').write_bytes(b'\0'.join(x.encode() for x in paths)+b'\0')
    if staged:
        actual=subprocess.check_output(['git','diff','--cached','--name-only','-z'],cwd=root).decode().split('\0')
        assert set(filter(None,actual))==set(paths)
        for name in paths:assert hashlib.sha256(subprocess.check_output(['git','show',':'+name],cwd=root)).hexdigest()==sha(root/name),name
    print(json.dumps(dict(bound_files=len(m['sha256']),exact_paths=len(paths),prefixes=6,full_goal_complete=False)))


if __name__=='__main__':
    if sys.argv[1]=='package':package(int(sys.argv[2]))
    check(sys.argv[1]=='staged')




"""Bind the actual dense GR return, its original failure and completed intervals."""
from pathlib import Path
from types import FunctionType
import importlib.util,shutil,sys
helper=Path('outputs/direct-eos-gr33/native-radau-pulse-gr/previous-publication-producer.py')
spec=importlib.util.spec_from_file_location('b',helper);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,runtime,master=b.root,b.runtime,b.master;read,write,sha=b.read,b.write,b.sha
out=root/'outputs/direct-eos-gr33/native-dense-gr-return';manifest=out.parent/'native-dense-gr-return-manifest.json'
note=root/'notes/REQUEST225_226_DENSE_GR_RETURN_KO.md'


def copy(src,dst):
    h=sha(src);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);assert sha(dst)==h==sha(src)


def package():
    old=runtime/'native-dense-return225-work';work=runtime/'native-interval-return226-work'
    row=read(work/'run-64.json');assert row['passed'] and read(work/'prefix-replay-identity.json')['passed']
    assert read(work/'coarse-receipt.json')['error'] is None
    failure=read(old/'coarse-receipt.json');assert 'Saved-stage anchor precision' in failure['error']
    previous=read(out.parent/'native-complete-radau-bridge-manifest.json');preserved={p:h for p,h in previous['sha256'].items() if not p.startswith('docs/')}
    for p,h in preserved.items():assert sha(root/p)==h,p
    assert not out.exists();out.mkdir();reused={}
    for src in old.rglob('*'):
        if not src.is_file():continue
        rel=src.relative_to(old)
        if rel.parts[0] in ['gr','sweep-0']:
            reused['225/'+rel.as_posix()]=dict(runtime_path=str(src),sha256=sha(src));continue
        copy(src,out/'original225'/rel)
    completed=[p for p in work.glob('*.json') if p.name not in ['controller-status.json','fine-receipt.json','audit-receipt.json','result.json']]
    completed += list(work.glob('expanded-*.py'))
    completed += list((work/'sweep-1/photons').glob('*64*'))
    for src in completed:copy(src,out/'coarse226'/src.relative_to(work))
    for src,name in [(root/'.phase225-followthrough.py','original-controller.py'),(root/'.phase226-followthrough.py','controller.py'),(Path(__file__),'publication-producer.py')]:copy(src,out/name)
    fine_complete=(work/'fine-receipt.json').exists() and read(work/'fine-receipt.json')['error'] is None
    paired=read(work/'result.json') if (work/'result.json').exists() else None
    if fine_complete:
        for src in [p for p in work.glob('*128*.json')]+list((work/'sweep-1/photons').glob('*128*'))+[work/'fine-receipt.json']:
            copy(src,out/'fine226'/src.relative_to(work))
    if paired:
        for name in ['result.json','audit-receipt.json']:copy(work/name,out/'paired226'/name)
    maximum=max(max(v['native_relative']) for v in row['anchor_checks'])
    note.write_text(f'''# 동일 결합 단계에 실제 밀집 GR 반환

분류: Counterexample candidate. **최종 전하의 결론은 미판정이다.** 수락된224의 같은T/8 물질·광자·에너지·출구 이력에서 만든 GR장을 실제 물질/광자 단계 방정식에 되먹였다. 현재 거친64시계15개 실제 단계가 원 수락 기준을 통과했다. 미세 경로 완료={fine_complete}. 짝 시간 판정={None if paired is None else paired['passed']}. 한 번의 반환이며 자기GR 고정점·완전 비선형 진화·전체 기간의 최종 전하 완료가 아니다.

분류: Counterexample candidate.225의 실제 metric은 원 시간/적분 차수/각도 광선 대조를 통과했다. 최초 물리 반환은T/16까지8단계를 통과한 뒤 다음 배경 구간에서 저장 native 반응 대조3.874713106e-12>1e-12로 실패했다.189.391초 이내의 오류였으며 예산 부족이 아니다. 원 실패·계획·실행 코드·수락된 체크포인트를 보존했다.

분류: Counterexample candidate. 원 물질 생산자가 각 canonical 구간에서 모델을 다시 만들었던 수명 규칙을 반환에도 적용했다. 같은225 방정식·정확 보존좌표 역변환·고정 GR 입력·격자·시계·물리 기간·수락 기준을 유지했다. 새 물리 모델이나 별도 전하를 더하지 않았다.225 체크포인트에는 누적 진단 최대값 일부가 없어 처음8단계를 한 번 재현했고, 모든 기존 체크포인트 배열과 정확히 같음을 확인했다. 이후 기존184 완전 재시작으로 상태·광자·물질·바닥 제거·단계율·각도/반경 출구 이력이 변경되지 않는지 검사하며 다음 구간을 실제로 진화했다.

분류: Counterexample candidate. 거친 경로의 최대 저장 native 대조는{maximum:.9e}, 실제 단계 잔차는{row['maximum_true_stage']:.9e}, 물리 모멘트 잔차는{row['maximum_true_physical_stage']:.9e}이다. 첫 구간 전환에서 발생한225 실패는 이 경로의 실제 속행에서 해결되었다. 미세 경로와 짝 시간 대조의 결과가 나오기 전에는 전체 반환을 수락하지 않는다. 이전188의T/64 국소 시간 실패는 보존하며, 더 긴 전역 노름의 성공으로 그 실패를 지우지 않는다.

분류: Proven. 같은 분기의 아핀 차분 항등식을 다시 기호 검사했다. 이는 정밀도·물리 폐쇄의 일반 정리가 아니다.

계산 예산은CPU1스레드·8GiB에서 거친30분/미세60분이며 원3Newton/8분기 반복과 정확도 기준을 유지한다. 재현 비용은 첫8단계의 기존172.40초를 기준으로 정했다. 고정된218/223/224 작업과 입력을 변경하지 않는다. 광자 공통 이력은 같은 해의 전체 출구 판독을 위해 계속 복원되고,224 소비자는 그 원 수락 후 자동 시작한다.

분류: Conjectural. EOS·미분·공간/외부 경계·전 기간·자기GR·완전 비선형·정적 비교·관측/무한대 전하를 닫아야 한다. 이번 실제 반환 성공을 최종 전하의 기존 부호나 결론 유지로 해석하지 않는다.
''',encoding='utf-8')
    prefixes={}
    for name in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi','dynamic-charge-completion']:
        p=root/f'docs/{name}.md';prefixes[p.relative_to(root).as_posix()]=dict(bytes=p.stat().st_size,sha256=sha(p))
        with p.open('ab') as f:f.write(('\n\n## 단계225–226 — 같은 실제 단계에 밀집 GR 반환\n\n분류: Counterexample candidate. 최종 전하는 미판정이다. 동일T/8 이력의 수락된 밀집 GR을 물질·광자 방정식에 실제 반환했다. 첫 배경 전환의 저장 native 재현 실패를 원 모델의 구간별 재생성으로 수정했고, 거친15단계를 원 기준으로 완료했다. 첫8단계와 재시작 이력의 정확 동일성을 확인했다. 이전 국소 시간 실패와 원 실패는 보존한다. 짝 시간·자기GR·전 기간·비선형·최종 전하의 전체 수락은 별도다. [근거](../notes/'+note.name+').\n').encode())
    final=dict(classification='Counterexample candidate',actual_same_joint_solution_GR_computed=True,actual_same_equation_GR_return_coarse_passed=True,coarse_result=row,fine_completed=fine_complete,paired_result=paired,original_failure_preserved=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(out/'result.json',final);write(out/'publication.json',dict(previous_nondoc=preserved,document_prefixes=prefixes,reused=reused))
    modules=[root/f'verification/{name}.py' for name in ['evolve_dense_same_history_return','continue_dense_gr_return']]
    for p in modules:assert sha(p)==sha(runtime/'verification'/p.name)
    files=modules+[note]+[root/p for p in prefixes]+[p for p in out.rglob('*') if p.is_file()]
    final.update(document_prefixes=prefixes,sha256={p.relative_to(root).as_posix():sha(p) for p in files});write(manifest,final)
    m=read(master);m['sha256'].update(final['sha256']);m['sha256'][manifest.relative_to(root).as_posix()]=sha(manifest);m['native_dense_GR_return']={k:v for k,v in final.items() if k!='sha256'};write(master,m)


check=FunctionType(b.check.__code__,dict(b.check.__globals__,out=out,manifest=manifest),argdefs=b.check.__defaults__)
if __name__=='__main__':
    if sys.argv[1]=='package':package()
    check(sys.argv[1]=='head')

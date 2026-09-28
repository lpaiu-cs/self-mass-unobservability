"""Save one completed pair and launch evidence while the frozen solve continues."""
from pathlib import Path
import importlib.util,shutil
spec=importlib.util.spec_from_file_location('admission',Path(__file__).with_name('.phase185-admission-publish.py'))
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
root,work,out=b.root,b.work,b.out;read,write,sha=b.read,b.write,b.sha
assert not (out/'first-continuation').exists()
pair=read(work/'comparison-02.json');cost=read(work/'cost-02.json');launch=read(work/'launch.json')
assert pair['passed'] and cost['eligible'] and pair['shared_interval']==2
assert all(r['accepted_prefix_preserved'] for r in pair['rows'])
assert sha(b.module)==launch['source_sha256'] and sha(work/'plan.json')==launch['plan_sha256']
folder=out/'first-continuation';folder.mkdir()
for name in ['comparison-02.json','cost-02.json','launch.json','launch-producer.ps1']:
    shutil.copyfile(work/name,folder/name);assert sha(work/name)==sha(folder/name)
shutil.copyfile(__file__,out/'continuation-publication-producer.py')
text='''

분류: Counterexample candidate. 전체 기간 실행을2026-09-24 14:51KST에 시작했다. 첫 추가 두 경로가0.429303890ms까지 실제 완료했고 원 시간 대조 최대0.243391%<2%, 네 물질 셀별 수지와 각도 출구, 이전 이력 보존을 통과했다. 두 새 경로 구간의 실제 비용은 합계378.811초였고 새 하위 단계는coarse7개·fine13개였다. 이 시각은 원 기간의12.5%이며 전체 연구 완성률이 아니다. 이후 구간은 같은 고정 계획으로 계속 진행하도록 했다. [첫 실제 연속 구간](../outputs/direct-eos-gr33/native-full-horizon-admission/first-continuation/comparison-02.json).

분류: Conjectural. 최초 실측 기반 예상은17:07KST전후, 실행 상한은19:21KST이다. 이는 물질·광자의 전체 기간 계산 일정이며 GR 반환과 최종 전하까지의 완료 예정시각이 아니다. 숨김 Windows WSL 부모 프로세스와 Linux실행 PID·시작 tick을 대조해 실행을 확인했다. 현재 상태는 런타임 `native-full-horizon185-work/status.json`과 실제 PID/시작 tick 및 종료 receipt를 함께 확인해야 한다. 원 코드·계획은 실행 중 변경하지 않는다.
'''
with b.note.open('ab') as f:f.write(text.encode('utf-8'))
manifest=read(b.manifest)
manifest['first_continuation']=dict(classification='Counterexample candidate',passed=True,
    same_horizon_seconds=pair['same_horizon_seconds'],photon_time=pair['photon_time'],material_time=pair['material_time'],
    first_pair_elapsed_seconds=pair['elapsed_seconds'],full_horizon_completed=False,final_charge_conclusion='unadjudicated')
manifest['scope']='Frozen admission snapshot plus first actual continuation. Later state requires the live process and runtime results; full-horizon success is not claimed.'
files=[b.note,out/'continuation-publication-producer.py']+list(folder.iterdir())
manifest['sha256'].update({p.relative_to(root).as_posix():sha(p) for p in files})
write(b.manifest,manifest)
master=read(b.master);master['sha256'].update(manifest['sha256']);master['sha256'][b.manifest.relative_to(root).as_posix()]=sha(b.manifest)
master['native_full_horizon_admission']={k:v for k,v in manifest.items() if k!='sha256'}
write(b.master,master);b.check()

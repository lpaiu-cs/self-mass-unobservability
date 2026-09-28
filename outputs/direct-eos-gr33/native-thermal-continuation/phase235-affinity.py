from pathlib import Path
import json,os,time
p=146514;root=Path(f'/proc/{p}');controller=json.loads(Path('native-thermal-continuation235-work/controller-start.json').read_text())
assert int(root.joinpath('stat').read_text().split()[3])==controller['pid']
assert b'verification/continue_thermal_native.py' in root.joinpath('cmdline').read_bytes()
assert 2 in os.sched_getaffinity(0)
threads=[int(v.name) for v in root.joinpath('task').iterdir()];before={t:sorted(os.sched_getaffinity(t)) for t in threads}
for tid in threads:os.sched_setaffinity(tid,{2})
after={t:sorted(os.sched_getaffinity(t)) for t in threads};assert all(v==[2] for v in after.values())
result=dict(pid=p,process_start_ticks=root.joinpath('stat').read_text().split()[21],boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),unix=time.time(),before=before,after=after,reason='232and235both pinned CPU0. Move235to CPU2 so the two independent one-CPU experiments can progress concurrently. Same one-CPU cap, no source, plan, arithmetic or acceptance change.')
Path('native-thermal-continuation235-work/resource-affinity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

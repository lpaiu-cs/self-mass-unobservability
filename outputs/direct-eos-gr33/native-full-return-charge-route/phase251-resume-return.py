from pathlib import Path
import json,os
out=Path('native-full-return249-work');root=Path('native-full-charge251-work/interruption-0605/249');root.mkdir()
read=lambda p:json.loads(Path(p).read_text())
entry=read(out/'controller-start.json');state=read(out/'controller-status.json')
assert entry['boot_id']!=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
assert state['state']=='waiting' and not state.get('completed')
assert not (out/'prepare-receipt.json').exists(),'Physical return already dispatched'
for name in ['controller-start.json','controller-status.json']:os.rename(out/name,root/name)

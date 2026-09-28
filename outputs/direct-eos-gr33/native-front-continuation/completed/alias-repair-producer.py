"""Restore the frozen linked gate from its SHA-bound published bytes."""
from pathlib import Path
import hashlib,json,shutil
root=Path(__file__).resolve().parent
runtime=Path('//wsl.localhost/Ubuntu-22.04/home/lpaiu/work/native-retained-tail-runtime')
work=runtime/'native-front-continuation184-work'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
binding=read(root/'outputs/direct-eos-gr33/native-material-front/publication.json')['reused']['check-result.json']
archive=root/binding['path'];assert sha(archive)==binding['sha256']
changed=work/'check-result.json';assert 'rows' in read(changed)
assert not (work/'restart-check.json').exists()
shutil.copyfile(changed,work/'restart-check.json')
before=sha(changed)
# Updating this shared inode restores every old runtime link to the immutable
# published gate. The new restart result now has its own inode and filename.
shutil.copyfile(archive,changed)
restored={}
for folder in ['native-full-incident181-work','native-full-time182-work','native-material-front183-work','native-front-continuation184-work']:
    p=runtime/folder/'check-result.json';assert sha(p)==binding['sha256'];restored[str(p.relative_to(runtime))]=sha(p)
shutil.copyfile(runtime/'verification/continue_full_material_front.py',work/'pre-alias-producer.py')
result=dict(classification='Counterexample candidate',
    failure='The inherited check-result.json was a hard link. Writing the new restart result at the same name changed its linked runtime gate. Detected before any new physical step or dispatch.',
    repair='Preserve the new result as restart-check.json on a separate inode. Restore the inherited runtime gates from their SHA-bound immutable published bytes. Future restart results use only the distinct filename.',
    overwritten_runtime_sha256=before,restored_runtime_files=restored,archive=binding,
    physical_arrays_unchanged=True,new_physical_steps=0)
(work/'alias-repair.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
shutil.copyfile(__file__,work/'alias-repair-producer.py')
print(json.dumps(result))

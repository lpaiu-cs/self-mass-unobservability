import os,json
from pathlib import Path
work=Path('retained-native-acoustic153-work');files=[];links={}
for directory,dirs,names in os.walk(work,followlinks=False):
    base=Path(directory)
    for name in dirs[:]:
        p=base/name
        if name=='__pycache__':dirs.remove(name)
        elif p.is_symlink():links[p.relative_to(work).as_posix()]=str(p.resolve());dirs.remove(name)
    for name in names:
        p=base/name
        if name.endswith('.pyc') or name=='publication-inventory.json':continue
        if p.is_symlink():links[p.relative_to(work).as_posix()]=str(p.resolve())
        else:files.append(p.relative_to(work).as_posix())
files.append('publication-inventory.json')
(work/'publication-inventory.json').write_text(json.dumps(dict(files=sorted(files),links=links),indent=2)+'\n')
print(json.dumps(dict(files=len(files),links=len(links),bytes=sum((work/p).stat().st_size for p in files))))

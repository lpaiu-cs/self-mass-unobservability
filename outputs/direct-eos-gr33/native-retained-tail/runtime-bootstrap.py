"""Byte-identical Linux code cache; no scientific state is changed."""
from pathlib import Path
import hashlib,json,sys,time,zipfile

ROOT=Path('E:/lab/self-mass-unobservability') if sys.platform=='win32' else Path('/mnt/e/lab/self-mass-unobservability')
OUT=ROOT/'outputs/direct-eos-gr33/native-retained-tail'
CACHE=Path('/home/lpaiu/work/native-retained-tail-runtime')

if sys.argv[1]=='pack':
    begin=time.monotonic();files=sorted((ROOT/'verification').glob('*.py'));hashes={}
    with zipfile.ZipFile(OUT/'runtime-source.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for p in files:
            content=p.read_bytes();rel=p.relative_to(ROOT).as_posix();hashes[rel]=hashlib.sha256(content).hexdigest();archive.writestr(rel,content)
    (OUT/'runtime-source.json').write_text(json.dumps(dict(files=hashes,seconds=time.monotonic()-begin),indent=2)+'\n')
    print(json.dumps(dict(files=len(files),seconds=time.monotonic()-begin)))
elif sys.argv[1]=='unpack':
    assert not CACHE.exists();CACHE.mkdir();begin=time.monotonic()
    with zipfile.ZipFile(OUT/'runtime-source.zip') as archive:archive.extractall(CACHE)
    for rel,h in json.loads((OUT/'runtime-source.json').read_text())['files'].items():
        assert hashlib.sha256((CACHE/rel).read_bytes()).hexdigest()==h
    (CACHE/'outputs').symlink_to(ROOT/'outputs',target_is_directory=True)
    print(json.dumps(dict(cache=str(CACHE),verified=True,seconds=time.monotonic()-begin)))
elif sys.argv[1]=='profile':
    import cProfile,pstats,io
    sys.path.insert(0,str(CACHE/'verification'));begin=time.monotonic();p=cProfile.Profile();p.enable()
    import def_native_cold_coupling
    p.disable();seconds=time.monotonic()-begin;stream=io.StringIO();pstats.Stats(p,stream=stream).sort_stats('cumulative').print_stats(20)
    (OUT/'runtime-import-profile.txt').write_text(stream.getvalue())
    result=dict(seconds=seconds,location=str(def_native_cold_coupling.__file__),new_native_states=0)
    (OUT/'runtime-import.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

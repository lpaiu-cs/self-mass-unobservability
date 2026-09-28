"""Preserve/recover large opacity text and frozen external input bytes."""
from pathlib import Path
import argparse
import gzip
import hashlib
import io
import json
import tarfile
import def_photon_actual_opacity as a


def seal():
    assert not (a.OUT/'archives.json').exists()
    names={r['path'] for r in json.loads((a.OUT/'source-download.json').read_text())['rows']}
    names|={r['path'] for r in json.loads((a.OUT/'metal-download.json').read_text())['files']}
    names|={r['path'] for r in json.loads((a.OUT/'build-balanced.json').read_text())['acquisition']['files']}
    names.add('selected-full.19')
    inputs=[]
    with tarfile.open(a.OUT/'inputs.tar.gz','w:gz') as tar:
        for name in sorted(names):
            raw=(a.CACHE/name).read_bytes()
            info=tarfile.TarInfo(name);info.size=len(raw);info.mtime=0
            tar.addfile(info,io.BytesIO(raw))
            inputs.append(dict(path=name,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    a.write(a.OUT/'inputs.json',dict(files=inputs,archive_sha256=a.digest(a.OUT/'inputs.tar.gz'),
        source_commit=a.COMMIT,attribution='SYNSPEC source MIT LICENSE included; atomic data and line-list provenance retained separately, not relicensed as MIT. The selected line-list bytes bind the public distribution receipt and selection receipt.'))
    rows=[]
    for p in sorted(a.OUT.rglob('*')):
        if not p.is_file() or p.stat().st_size<1_000_000 or p.name not in ['fort.27','fort.28','fort.29','fort.88','fort.89','opacity.txt','stdout.log']:continue
        assert p.resolve().is_relative_to(a.OUT.resolve())
        raw=p.read_bytes();target=p.with_name(p.name+'.gz')
        assert not target.exists();target.write_bytes(gzip.compress(raw,compresslevel=6,mtime=0))
        assert gzip.decompress(target.read_bytes())==raw
        rows.append(dict(path=str(p.relative_to(a.OUT)),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
            archive=str(target.relative_to(a.OUT)),archive_bytes=target.stat().st_size,archive_sha256=a.digest(target)))
        p.unlink()
    a.write(a.OUT/'archives.json',dict(files=rows,raw_bytes=sum(x['bytes'] for x in rows),archive_bytes=sum(x['archive_bytes'] for x in rows)))
    print('ARCHIVED',len(rows),sum(x['bytes'] for x in rows),sum(x['archive_bytes'] for x in rows))


def unpack():
    for row in json.loads((a.OUT/'archives.json').read_text())['files']:
        source=a.OUT/row['archive'];target=a.OUT/row['path']
        assert target.resolve().is_relative_to(a.OUT.resolve()) and a.digest(source)==row['archive_sha256']
        raw=gzip.decompress(source.read_bytes());assert hashlib.sha256(raw).hexdigest()==row['sha256']
        if target.exists():assert a.digest(target)==row['sha256']
        else:target.write_bytes(raw)
    info=json.loads((a.OUT/'inputs.json').read_text())
    assert a.digest(a.OUT/'inputs.tar.gz')==info['archive_sha256']
    rows={r['path']:r for r in info['files']}
    with tarfile.open(a.OUT/'inputs.tar.gz','r:gz') as tar:
        for member in tar:
            assert member.isfile() and member.name in rows
            target=a.CACHE/member.name;assert target.resolve().is_relative_to(a.CACHE.resolve())
            raw=tar.extractfile(member).read();assert hashlib.sha256(raw).hexdigest()==rows[member.name]['sha256']
            if target.exists():assert a.digest(target)==rows[member.name]['sha256']
            else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    print('Restored checked spectra and frozen external inputs.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['seal','unpack']);globals()[p.parse_args().action]()

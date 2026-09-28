"""Recover complete CRC-checked local ZIP entries; never infer missing data."""
from pathlib import Path
import binascii,hashlib,json,struct,zipfile,zlib

ROOT=Path('E:/lab/self-mass-unobservability')
OUT=ROOT/'outputs/direct-eos-gr33/native-retained-tail/supported-temperature'
EV=OUT/'evolution';source=EV/'checkpoint-128.tmp';data=source.read_bytes();pos=0;rows={};failure=None
try:
    while pos+30<=len(data) and data[pos:pos+4]==b'PK\x03\x04':
        _,version,flags,method,mt,md,crc,cs,us,nl,el=struct.unpack_from('<4s5H3I2H',data,pos)
        name=data[pos+30:pos+30+nl].decode();extra=data[pos+30+nl:pos+30+nl+el];i=0
        assert not flags&8, 'Unfinished streaming ZIP entry'
        while i+4<=len(extra):
            tag,size=struct.unpack_from('<HH',extra,i);value=extra[i+4:i+4+size];i+=4+size
            if tag==1:
                off=0
                if us==0xffffffff:us=struct.unpack_from('<Q',value,off)[0];off+=8
                if cs==0xffffffff:cs=struct.unpack_from('<Q',value,off)[0]
        at=pos+30+nl+el;assert at+cs<=len(data),('Truncated compressed entry',name)
        packed=data[at:at+cs];raw=zlib.decompress(packed,-15) if method==8 else packed
        assert method in [0,8] and len(raw)==us and (binascii.crc32(raw)&0xffffffff)==crc,('CRC',name)
        rows[name]=raw;pos=at+cs
    with zipfile.ZipFile(EV/'preproduction-checkpoint-128.npz') as z:expected=set(z.namelist())
    assert set(rows)==expected,('Missing fields',sorted(expected-set(rows)))
    # Every separately captured physical field must be byte-identical.
    with zipfile.ZipFile(EV/'state-128-32.npz') as z:
        checked=[k+'.npy' for k in ['U','I','u','theta','eta','Pi','h','j']]
        for name in checked:assert rows[name]==z.read(name),('Independent canonical state',name)
    import numpy as np,io
    assert int(np.load(io.BytesIO(rows['completed.npy'])))==32
    target=EV/'recovered-checkpoint-128-32.npz';assert not target.exists()
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED) as z:
        for name,raw in rows.items():z.writestr(name,raw)
    with zipfile.ZipFile(target) as z:assert z.testzip() is None and set(z.namelist())==expected
except Exception as exc:failure=repr(exc)
result=dict(passed=failure is None,error=failure,complete_entries=len(rows),bytes=len(data),consumed=pos,
    truncated_source_sha256=hashlib.sha256(data).hexdigest(),inferred_scientific_fields=0,
    recovery='Only complete local ZIP entries with matching CRC and full original field set are eligible; compare physical arrays to the independent32-step snapshot.',
    checkpoint_step=32 if failure is None else None)
(OUT/'checkpoint-recovery.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

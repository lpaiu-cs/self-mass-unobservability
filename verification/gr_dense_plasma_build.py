"""Compile the frozen derivative correction with its explicit 83-column line."""
import json,subprocess,sys
import gr_dense_plasma_defined as d


def run():
    d.bindings();out=d.OUT/'screening';before=json.loads((out/'build.json').read_text())
    assert before['returncode']==1 and not (out/'build-corrected.json').exists()
    source=out/'eos22-screening.f';lines=source.read_text().splitlines()
    long=[(i+1,s) for i,s in enumerate(lines) if s and s[0] not in '*Cc!' and len(s.split('!')[0])>72]
    assert len(long)==1 and long[0][0]==738 and 'H1*((2.d0/5.d0)' in long[0][1]
    original=(d.d.OUT/'potekhin-chabrier-eos22.f').read_text().splitlines()
    assert not [s for s in original if s and s[0] not in '*Cc!' and len(s.split('!')[0])>72]
    command=before['command'].copy();command.insert(1,'-ffixed-line-length-none')
    result=subprocess.run(command,capture_output=True,text=True)
    record=dict(classification='Counterexample candidate',command=command,returncode=result.returncode,
        stdout=result.stdout,stderr=result.stderr,only_extended_executable_line=long,
        bindings={p.relative_to(d.g.ROOT).as_posix():d.g.c.sha(p) for p in [
            d.g.ROOT/'verification/gr_dense_plasma_build.py',source,out/'build.json',d.OUT/'plan.json']})
    (out/'build-corrected.json').write_text(json.dumps(record,indent=2)+'\n');assert result.returncode==0
    runtime=dict(sha256={str(d.CACHE/'pc.so'):d.g.c.sha(d.CACHE/'pc.so')},
        inherited_runtime_sha256=json.loads((d.d.OUT/'runtime.json').read_text())['sha256'])
    (out/'runtime.json').write_text(json.dumps(runtime,indent=2)+'\n');verify()


def verify():
    record=json.loads((d.OUT/'screening/build-corrected.json').read_text());assert record['returncode']==0
    for rel,digest in record['bindings'].items():assert d.g.c.sha(d.g.ROOT/rel)==digest,rel
    print('PASS frozen screening source build; original 72-column compilation failure retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

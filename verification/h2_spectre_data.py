"""Build the sourced H2SPECTRE calculator and retain every finite grid result."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path, PurePosixPath
import hashlib, json, re, shutil, subprocess, sys, tarfile, time
import eos_spectral_build_audit as build_audit

g=build_audit.g;OUT=g.OUT/'gr-h2-spectre-audit';DATA=g.OUT/'gr-h2-spectre-source'
CACHE=g.CACHE/'h2-spectre';SOURCE=CACHE/'source'
NMAX=[31,30,28,27,25,23,22,20,18,16,14,12,10,7,4]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not (OUT/'plan.json').exists() and not (SOURCE/'H2spectre.exe').exists()
    OUT.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True);SOURCE.mkdir(parents=True,exist_ok=True);build_audit.verify()
    names=['h2spectre_7.4.tar.gz','h2-pachucki-codes.html','h2-napt2019.pdf',
        'h2-qed2011-metadata.json','h2-qed2011-supplement.pdf',
        'h2-partition2021-data-metadata.json','h2-partition2021-data.zip',
        'h2-partition2021-supplement-metadata.json','h2-partition2021-supplement.pdf']
    for name in names:
        original=g.ROOT/'outputs'/name;target=DATA/name
        if target.exists():assert target.read_bytes()==original.read_bytes()
        else:shutil.copy2(original,target)
    for name,md5 in [('h2-qed2011-supplement.pdf','dfc5819c6bb1ed6da19c52cc3e063a4e'),
        ('h2-partition2021-data.zip','347374f212921eab40b7a913dd2907a0'),
        ('h2-partition2021-supplement.pdf','b5b949b93d21b2c78822840db8764bf9')]:
        assert hashlib.md5((DATA/name).read_bytes()).hexdigest()==md5
    members=[]
    with tarfile.open(DATA/'h2spectre_7.4.tar.gz') as archive:
        for member in archive.getmembers():
            rel=PurePosixPath(member.name)
            assert not rel.is_absolute() and '..' not in rel.parts and '\\' not in member.name
            target=(SOURCE/str(rel)).resolve();assert target.is_relative_to(SOURCE.resolve())
            row=dict(name=member.name,type=member.type.decode(),size=member.size,link=member.linkname)
            if member.isdir():target.mkdir(exist_ok=True,parents=True);row['action']='directory'
            elif member.islnk():
                assert member.linkname==member.name and any(r['name']==member.name and r['type']=='0' for r in members)
                assert target.is_file() or target.suffix.lower() in ['.o','.a','.mod','.exe']
                row['action']='skip duplicate self-hardlink'
            elif member.issym():
                assert member.name=='h2spectr_pot.f90' and member.linkname=='data/h2spectr_pot.f90'
                row['action']='materialize verified in-archive potential source after regular files'
            else:
                assert member.isfile()
                if target.suffix.lower() in ['.o','.a','.mod','.exe']:
                    row['action']='skip upstream prebuilt object'
                else:
                    raw=archive.extractfile(member).read();target.parent.mkdir(exist_ok=True,parents=True)
                    if target.exists():assert target.read_bytes()==raw
                    else:target.write_bytes(raw)
                    row['action']='extract regular data/source'
            members.append(row)
    shutil.copy2(SOURCE/'data/h2spectr_pot.f90',SOURCE/'h2spectr_pot.f90')
    save('archive-members.json',members)
    main=(SOURCE/'h2spectre.f90').read_text()
    match=re.search(r'H2vJ=\(/([\d,]+)/\)',main);assert match and list(map(int,match[1].split(',')))==NMAX
    assert sum(n+1 for n in NMAX)==302
    for path in SOURCE.rglob('*'):
        if path.suffix.lower() in ['.f90','.f']:
            code='\n'.join(line.split('!')[0] for line in path.read_text().splitlines()
                if not (path.suffix.lower()=='.f' and line[:1] in ['*','c','C']))
            assert not re.search(r'\bcall\s+(system|execute_command_line)\b',code,re.I),path
    save('source-tree.json',{p.relative_to(SOURCE).as_posix():g.c.sha(p) for p in SOURCE.rglob('*') if p.is_file()})
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3d70a3b',
        source=dict(classification='Imported from prior work',version='7.4, upstream date 12.10.2022',
            archive='https://www.fuw.edu.pl/~krp/h2spectre_7.4.tar.gz',
            author_index='https://www.fuw.edu.pl/~krp/codes.html',
            methods='Komasa et al. 2019, Phys. Rev. A 100,032519, DOI 10.1103/PhysRevA.100.032519',
            methods_pdf='https://www.fuw.edu.pl/~krp/papers/H2spectr19.pdf',
            additional_references=['10.1021/ct200438t','10.1021/acs.jpca.1c06468'],
            support='The upstream check_vJ table specifies 302 H2 X-state bound keys, including the explicitly warned weak (14,4) state. This support is imported, not a proof of completeness or bound-state existence.',
            physical_uncertainty='Upstream componentwise truncation/potential error estimates are not rigorous bounds; README explicitly warns about high/weak states. CODATA versions in four-body inputs and current potentials differ as documented upstream.'),
        baseline_grids=[[200,10],[400,10],[400,20]],Nmax_by_v=NMAX,
        finite_energy_difference_gate_cm_inverse=1e-4,
        baseline_interpretation='Compare 200/10 versus 400/10 for spacing, and 400/10 versus 400/20 for range, separately. Retain nonpositive energies and all failures. No failed baseline level is silently dropped from support.',
        weak_state_grids=[[800,40],[1000,50],[1200,60]],
        execution='Unmodified upstream Makefile and Fortran/data; skip all upstream prebuilt objects. Run native Levels mode, default mixed FULL/NAPT policy, -V 1 -P 12. Each output binds exact stdin, command, elapsed time and return code. Four workers on spare CPUs 12-15. No forced -Dv validity override.',
        scope='Finite numerical spectral calculation. No interval eigensolver, physical EOS, continuous derivative, plasma occupation, continuum or excited-electronic completeness certificate.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/h2_spectre_data.py',OUT/'archive-members.json',OUT/'source-tree.json',*[DATA/n for n in names]]}))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'source-tree.json').read_text()).items():assert g.c.sha(SOURCE/rel)==digest,rel
    return plan


def build():
    bindings();assert not (OUT/'build.json').exists()
    result=subprocess.run(['make'],cwd=SOURCE,capture_output=True)
    (OUT/'build.log').write_bytes(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr.decode(errors='replace')
    command=[str(SOURCE/'H2spectre.exe'),'--PC']
    result=subprocess.run(command,cwd=SOURCE,capture_output=True);assert result.returncode==0
    (OUT/'physical-constants.txt').write_bytes(result.stdout+result.stderr)
    save('build.json',dict(command=['make'],executable=str(SOURCE/'H2spectre.exe'),
        sha256=g.c.sha(SOURCE/'H2spectre.exe'),compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
        build_log_sha256=g.c.sha(OUT/'build.log'),constants_sha256=g.c.sha(OUT/'physical-constants.txt')))


def parse(text,keys):
    parts=re.split(r"v'\s*=\s*(\d+)\s+J'\s*=\s*(\d+)",text);rows=[]
    assert len(parts)==1+3*len(keys),(len(parts),len(keys))
    for a in range(1,len(parts),3):
        v,J=map(int,parts[a:a+2]);body=parts[a+2]
        match=re.search(r'Total\s+(-?\d+\.\d+)\s+(\d\.\dE[+-]\d+)',body)
        method=re.search(r'E2\((\w+)\s*\)',body);assert match and method,body
        D,error=match.groups();assert float(error)>=0
        rows.append(dict(v=v,J=J,dissociation_cm_inverse=D,upstream_estimated_uncertainty_cm_inverse=error,
            E2_method=method[1],positive=float(D)>0,weak_warning='Weakly bound' in body))
    assert [(r['v'],r['J']) for r in rows]==keys
    return rows


def job(args):
    N,R,keys,label=args;build_record=json.loads((OUT/'build.json').read_text())
    assert g.c.sha(build_record['executable'])==build_record['sha256']
    assert not (OUT/(label+'.json')).exists() and not (OUT/(label+'.log')).exists()
    stdin='H2 Levels\n'+''.join(f'{v} {J}\n' for v,J in keys)
    (OUT/(label+'.input')).write_text(stdin)
    command=[build_record['executable'],'-N',str(N),'-R',str(R),'-V','1','-P','12']
    started=time.monotonic();result=subprocess.run(command,cwd=SOURCE,input=stdin.encode(),capture_output=True)
    raw=result.stdout+result.stderr;(OUT/(label+'.log')).write_bytes(raw)
    record=dict(classification='Counterexample candidate',command=command,elapsed_seconds=time.monotonic()-started,
        returncode=result.returncode,input_sha256=g.c.sha(OUT/(label+'.input')),raw_sha256=g.c.sha(OUT/(label+'.log')))
    save(label+'-execution.json',record);assert result.returncode==0
    rows=parse(raw.decode(),keys);save(label+'.json',dict(**record,rows=rows))
    print('H2 SPECTRE',label,len(rows),'positive',sum(r['positive'] for r in rows),round(record['elapsed_seconds'],2),flush=True)
    return rows


def run():
    plan=bindings();jobs=[]
    for N,R in plan['baseline_grids']:
        for v,n in enumerate(NMAX):jobs.append((N,R,[(v,J) for J in range(n+1)],f'grid-{N}-{R}-v{v:02}'))
    for N,R in plan['weak_state_grids']:jobs.append((N,R,[(14,4)],f'weak-{N}-{R}'))
    with ProcessPoolExecutor(max_workers=4) as pool:
        for done in as_completed([pool.submit(job,args) for args in jobs]):done.result()
    collect()


def collect():
    plan=bindings();grids=[]
    for N,R in plan['baseline_grids']:
        rows=[r for v in range(15) for r in json.loads((OUT/f'grid-{N}-{R}-v{v:02}.json').read_text())['rows']]
        assert len(rows)==302;grids.append(rows)
    rows=[]
    for a,b,c in zip(*grids):
        gap1=abs(float(a['dissociation_cm_inverse'])-float(b['dissociation_cm_inverse']))
        gap2=abs(float(b['dissociation_cm_inverse'])-float(c['dissociation_cm_inverse']))
        rows.append(dict(v=c['v'],J=c['J'],spacing_difference_cm_inverse=gap1,range_difference_cm_inverse=gap2,
            baseline_positive=c['positive'],baseline_passed=c['positive'] and max(gap1,gap2)<plan['finite_energy_difference_gate_cm_inverse']))
    save('baseline-result.json',dict(classification='Counterexample candidate',keys=302,
        passed_levels=sum(r['baseline_passed'] for r in rows),failed_levels=[r for r in rows if not r['baseline_passed']],
        maximum_spacing_difference=max(r['spacing_difference_cm_inverse'] for r in rows),
        maximum_range_difference=max(r['range_difference_cm_inverse'] for r in rows),rows=rows,
        complete_physical_spectrum_certified=False))
    print('H2 BASELINE',sum(r['baseline_passed'] for r in rows),'/ 302',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

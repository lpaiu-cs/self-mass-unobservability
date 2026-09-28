"""Request 14 producers. Run in WSL; never mutate earlier evidence/builds."""
import difflib
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/validated-variational'
RUNTIME = Path.home() / 'work/nutimo_pilot'
SOURCE = RUNTIME / 'nutimo_request13_stable/src'
TARGET = RUNTIME / 'nutimo_request14_export/src'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not TARGET.exists(), 'Do not overwrite a previous build'
    shutil.copytree(SOURCE, TARGET)
    path = TARGET / 'AllTheories3Bodies.cpp'
    old = path.read_text()
    anchor = '// Do the backward integration\n    Retro_Integre();'
    assert old.count(anchor) == 1
    new = old.replace(anchor, '''// Request 14: read-only IVP export before either integration direction.
    if (const char* path14 = getenv("TIMING14_EXPORT")) {
        FILE* f14 = fopen(path14, "w");
        if (!f14) throw runtime_error("Cannot open Request 14 IVP export");
        fprintf(f14, "%d\\n", nbody_plus_extra);
        // Hexadecimal long doubles serialize the binary values exactly.
        fprintf(f14, "%La %La %La %La %La %La %La\\n", t0, ts.front(), ts.back(), length, timescale, clightAd, clightAd2);
        for (auto v : x0) fprintf(f14, "%La\\n", v);
        for (int i=0;i<nbody_plus_extra;++i) fprintf(f14, "%La\\n", Ms[i]);
        for (int i=0;i<nbody_plus_extra;++i) for (int j=0;j<nbody_plus_extra;++j) fprintf(f14, "%La %La\\n", Gg[i][j], gammabar[i][j]);
        for (int i=0;i<nbody_plus_extra;++i) for (int j=0;j<nbody_plus_extra;++j) for (int k=0;k<nbody_plus_extra;++k) fprintf(f14, "%La\\n", betabar[i][j][k]);
        fprintf(f14, "%La %La %La %La\\n", SEPdyn_A, SEPdyn_w, SEPdyn_ph, SEPdyn_tau);
        fclose(f14);
    }
''' + anchor)
    begin = new.index('void Integrateur::rhs_GR_nbody(')
    end = new.index(' #undef GGEFF', begin)
    part = new[begin:end]
    anchor = '\n\n };'
    assert part.count(anchor) == 1
    part = part.replace(anchor, '''
    static int samples14 = 0;
    if (const char* path14 = getenv("TIMING14_RHS")) {
        if (samples14 < 64) {
            FILE* f14 = fopen(path14, samples14 == 0 ? "w" : "a");
            if (!f14) throw runtime_error("Cannot open Request 14 RHS export");
            fprintf(f14, "%La", t);
            for (auto v : xvect) fprintf(f14, " %La", v);
            for (auto v : dxdtvect) fprintf(f14, " %La", v);
            fprintf(f14, "\\n"); fclose(f14); ++samples14;
        }
    }
''' + anchor)
    new = new[:begin] + part + new[end:]
    path.write_text(new)
    (OUT/'export.patch').write_text(''.join(difflib.unified_diff(old.splitlines(True), new.splitlines(True), fromfile='request13_stable/AllTheories3Bodies.cpp', tofile='request14_export/AllTheories3Bodies.cpp')))
    command = json.loads((ROOT/'outputs/research-remediation/stable-build.json').read_text())['command']
    with (OUT/'export-build.log').open('w') as log:
        subprocess.run(command, cwd=TARGET, stdout=log, stderr=subprocess.STDOUT, check=True)
    run = RUNTIME/'run_request14_export'
    assert not run.exists()
    shutil.copytree(RUNTIME/'run_request13_stable', run)
    env = dict(os.environ, LD_LIBRARY_PATH=str(TARGET), PYTHONPATH=str(TARGET),
               TEMPO2=str(RUNTIME/'install/third_party/tempo2'),
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2',
               TIMING14_EXPORT=str(OUT/'ivp.hex'), TIMING14_RHS=str(OUT/'native-rhs.hex'))
    for key in ['SEPDYN_A','SEPDYN_W','SEPDYN_PH','SEPDYN_TAU']:
        env.pop(key, None)
    code = "import python_Fittriple_interface as p; f=p.PyFittriple('parfile-planetGR-max-bestfit','0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'); print('EXPORT_COMPLETE')"
    with (OUT/'export-run.log').open('w') as log:
        subprocess.run([sys.executable, '-c', code], cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    (OUT/'export.json').write_text(json.dumps(dict(source_sha256=sha(SOURCE/'AllTheories3Bodies.cpp'), export_source_sha256=sha(path), library_sha256=sha(TARGET/'libFittriplecpp.so'), command=command, ivp_sha256=sha(OUT/'ivp.hex'), rhs_sha256=sha(OUT/'native-rhs.hex'),par_sha256=sha(run/'parfile-planetGR-max-bestfit'),tim_sha256=sha(run/'0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim')), indent=2)+'\n')


def build():
    import shlex
    import math
    from fractions import Fraction
    controls={}
    for name,odd,bounds in [('cos',0,('0x1.14a280fb5068bp-1','0x1.14a280fb5068dp-1')),('sin',1,('0x1.aed548f090cedp-1','0x1.aed548f090cefp-1'))]:
        partial=sum((Fraction((-1)**k,math.factorial(2*k+odd)) for k in range(20)),Fraction(0))
        other=partial+Fraction(1,math.factorial(40+odd))
        lo,hi=map(float.fromhex,bounds)
        assert Fraction(lo)<=min(partial,other)<=max(partial,other)<=Fraction(hi)
        controls[name]=dict(hex_bounds=bounds,exact_series_low=str(partial),exact_series_high=str(other))
    (OUT/'positive-control.json').write_text(json.dumps(controls,indent=2)+'\n')
    capd = RUNTIME/'request14_capd'
    flags = subprocess.check_output([str(capd/'build/bin/capd-config'),'--cflags','--libs'],text=True).strip()
    binary = RUNTIME/'request14_variational'
    command=['g++',str(ROOT/'verification/validated_variational.cpp'),*shlex.split(flags),'-o',str(binary)]
    with (OUT/'capd-client-build.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    (OUT/'capd-build.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=capd,text=True).strip(),flags=flags,command=command,binary_sha256=sha(binary),library_sha256=sha(capd/'build/libcapd.a'),source_sha256=sha(ROOT/'verification/validated_variational.cpp')),indent=2)+'\n')
    subprocess.run([str(binary)],check=True)


def seal():
    source=OUT/'sources';source.mkdir(exist_ok=True)
    capd=RUNTIME/'request14_capd'
    for path in [SOURCE/'AllTheories3Bodies.cpp',SOURCE/'Constants.h',capd/'build/CMakeFiles/capd.dir/flags.make',capd/'build/CMakeCache.txt',capd/'COPYING']:
        if path.exists():shutil.copy2(path,source/path.name)
    shutil.copy2(ROOT/'verification/validated_variational.cpp',source/'final-client.cpp')
    runs={}
    for name,horizon,mode,build,client in [
        ('forward',.07,'rect2','rect2','initial-rect2'),
        ('backward',-.068488,'rect2','rect2','initial-rect2'),
        ('full-forward',124.162279,'rect2','rect2','initial-rect2'),
        ('ho-forward',124.162279,'ho','ho','rect2-ho'),
        ('scaled-forward',124.162279,'scaled','scaled','scaled-client'),
        ('scaled-local',.07,'scaled','scaled','scaled-client'),
        ('scaled-backward',-.068488,'scaled','jacobi-adaptive','jacobi-adaptive'),
        ('jacobi-forward',124.162279,'jacobi','jacobi-adaptive','jacobi-adaptive'),
        ('jacobi-fixed-forward',124.162279,'jacobi','jacobi-fixed','jacobi-fixed'),
        ('jacobi-translation-forward',124.162279,'jacobi','jacobi-translation','jacobi-translation'),
        ('jacobi-local',.07,'jacobi','jacobi-local','jacobi-local'),
        ('mass-local',.07,'mass',None,'final-client')]:
        p=OUT/('capd-build'+('-'+build if build else '')+'.json')
        record=json.loads(p.read_text());client_path=source/(client+'.cpp')
        assert sha(client_path)==record['source_sha256'],name
        runs[name]=dict(horizon_internal=horizon,mode=mode,build_record=p.name,source_snapshot=str(client_path.relative_to(OUT)),source_sha256=sha(client_path),binary_sha256=record['binary_sha256'])
    (OUT/'run-index.json').write_text(json.dumps(runs,indent=2)+'\n')
    record=json.loads((OUT/'export.json').read_text());run=RUNTIME/'run_request14_export'
    record.update(par_sha256=sha(run/'parfile-planetGR-max-bestfit'),tim_sha256=sha(run/'0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'))
    (OUT/'export.json').write_text(json.dumps(record,indent=2)+'\n')
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    paths += [ROOT/'verification'/n for n in ['validated_variational.cpp','validated_variational.py','variational_audit.py']]
    paths += [ROOT/'notes/REQUEST14_VALIDATED_VARIATIONAL.md']
    documents=['docs/'+n+'.md' for n in ['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi']]
    paper=ROOT/'paper/revision-manifest.json';revision=json.loads(paper.read_text())
    for name in documents:revision['sha256'][name]=sha(ROOT/name)
    revision['request14_supporting_note_update']=dict(paths=documents,evidence_manifest='outputs/validated-variational/manifest.json',status='Conditional GR IVP/initial-state variational enclosures; full timing D2 remains false',artifact_status='The manuscript PDF/source ZIP remain the frozen Request 12 artifacts; they do not contain a completed Request 14 certificate.')
    paper.write_text(json.dumps(revision,indent=2,ensure_ascii=False)+'\n')
    paths += [ROOT/n for n in documents]+[paper]
    result=dict(baseline='874d7c9',before_task_checkpoint='b62a460',scope='Conditional validated flow enclosures; full D2 remains false',sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)})
    (OUT/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print('sealed',len(paths),'files')


def archive_initial_client():
    text=(OUT/'sources/rect2-ho.cpp').read_text()
    start=text.index('    capd::IMap nonlinear(')
    end=text.index('\n}',start)
    text=text[:start]+text[end+1:]
    text=text.replace('if(argc!=5&&argc!=6) throw std::runtime_error("args: ivp.hex native-rhs.hex horizon output.jsonl [ho]");\n        const bool ho=argc==6&&std::string(argv[5])=="ho";\n        if(argc==6&&!ho)throw std::runtime_error("unknown set representation");','if(argc!=5) throw std::runtime_error("args: ivp.hex native-rhs.hex horizon output.jsonl");')
    text=text.replace('capd::ITimeMap tm(solver);tm.stopAfterStep(true);','capd::ITimeMap tm(solver);tm.stopAfterStep(true);capd::C1Rect2Set s(u);')
    text=text.replace('auto integrate=[&](auto& s) {','try {')
    text=text.replace('        };\n        try {\n            if(ho){capd::C1HORect2Set s(u);integrate(s);}\n            else {capd::C1Rect2Set s(u);integrate(s);}\n','')
    expected=json.loads((OUT/'capd-build-rect2.json').read_text())['source_sha256']
    assert hashlib.sha256(text.encode()).hexdigest()==expected
    (OUT/'sources/initial-rect2.cpp').write_text(text)


def check():
    import tempfile
    binary=RUNTIME/'request14_variational'
    base=(OUT/'ivp.hex').read_text().split()
    first=(OUT/'native-rhs.hex').read_text().splitlines()[0].split()
    assert first[0]==base[1] and first[1:25]==base[8:32], 'RHS samples must start at the exported IVP'
    results={}
    with tempfile.TemporaryDirectory(prefix='request14-') as folder:
        folder=Path(folder)
        for label,index,value in [('non_GR_pair',38,'0x1.0001p+0'),('dynamic_coupling',len(base)-4,'0x1p-20'),('nonfinite_state',8,'nan')]:
            tokens=base.copy();tokens[index]=value
            p=folder/'invalid.hex';p.write_text(' '.join(tokens))
            run=subprocess.run([str(binary),str(p),str(OUT/'native-rhs.hex'),'0',str(folder/'test.jsonl')],capture_output=True,text=True)
            assert run.returncode!=0 and 'requested_horizon_completed' not in run.stdout,label
            results[label]=dict(rejected=True,reason=run.stderr.strip())
        p=folder/'bad-rhs.hex';lines=(OUT/'native-rhs.hex').read_text().splitlines();tokens=lines[0].split();tokens[25]='0x1p+20';lines[0]=' '.join(tokens);p.write_text('\n'.join(lines)+'\n')
        run=subprocess.run([str(binary),str(OUT/'ivp.hex'),str(p),'0',str(folder/'test.jsonl')],capture_output=True,text=True)
        assert run.returncode!=0 and 'native RHS regression failed' in run.stderr
        results['corrupted_native_rhs']=dict(rejected=True,reason=run.stderr.strip())
    (OUT/'negative-controls.json').write_text(json.dumps(results,indent=2)+'\n')
    print('PASS: non-GR, nonzero signal, nonfinite input and corrupted native RHS are rejected')


def report():
    from decimal import Decimal, localcontext, ROUND_CEILING
    def upper(value):
        with localcontext() as ctx:
            ctx.prec=8;ctx.rounding=ROUND_CEILING
            return str(+Decimal.from_float(value))
    audit=json.loads((OUT/'audit.json').read_text())
    scale=float.fromhex((OUT/'ivp.hex').read_text().split()[5])/86400
    local=audit['runs']['jacobi-local']['checkpoints'][-1]
    text='\n## Recorded results\n\n'
    text+=f"Status: Proven. At the binary64 internal epoch 0.07 (approximately {0.07*scale:.8f} days), the Jacobi run encloses the 24-state flow with midpoint state error norm at most {upper(local['state_midpoint_error_l2_upper'])} and initial-state Jacobian operator error at most {upper(local['initial_state_Jacobian_midpoint_error_operator_upper'])}, in the original dimensionless coordinates. The geometric-delay box width is at most {upper(local['geometric_delay_box_width_us_upper'])} microseconds. These displayed upper bounds are rounded upward. Raw IVP Jacobian norms cannot be compared directly to singular values of the normalized timing nuisance matrix.\n\n"
    text+='Status: Proven. The Cartesian, scaled and Jacobi state/Jacobian intervals at this common epoch intersect component by component. The augmented mass calculation also agrees with the scaled state block. The four independent fractional-mass column error norms, in pulsar/inner/outer/extra order, are bounded by '+', '.join(upper(v) for v in audit['independent_fractional_mass_columns']['error_l2_upper'])+'.\n\n'
    text+='Status: Imported from prior work. The full observation interval extends to approximately '+f'{124.16227849383201*scale:.3f}'+' days after the reference epoch. Every long-run representation below stopped at the declared Jacobian-width ceiling, without restarting its uncertainty. The listed last epochs are not useful-accuracy guarantees over the whole preceding interval.\n\n'
    text+='| Representation | Last returned epoch (days) | Full horizon completed |\n|---|---:|---|\n'
    for name in ['full-forward','ho-forward','scaled-forward','jacobi-translation-forward']:
        run=audit['runs'][name];last=run['checkpoints'][-1]
        text+=f"| {name} | {sum(last['time_internal'])/2*scale:.8f} | {run['requested_horizon_completed']} |\n"
    text+='\nStatus: Conjectural. The next dependency is a long-span variational representation that avoids this interval overestimation, with the same inclusion test retained. In parallel mathematical preparation, the timing-parameter-to-IVP chain and full delay/inverse-time map must be defined. Extra numerical precision or Kepler-based coordinates are candidate remedies, not established solutions. D2 and full physical inference remain open.\n'
    path=ROOT/'notes/REQUEST14_VALIDATED_VARIATIONAL.md'
    before=path.read_text().split('\n## Recorded results\n')[0]
    path.write_text(before+text)
    print(text)


def verify():
    count=0
    for path in [ROOT/'outputs/research-remediation/manifest.json',OUT/'manifest.json',ROOT/'paper/revision-manifest.json']:
        for name,expected in json.loads(path.read_text())['sha256'].items():
            assert sha(ROOT/name)==expected,name
            count+=1
    audit=json.loads((OUT/'audit.json').read_text())
    assert audit['D2']['pass_gate'] is False
    assert audit['D2']['local_independent_mass_partials_certified'] is True
    assert audit['coordinate_enclosures_intersect_at_local_epoch'] is True
    print('PASS:',count,'manifest bindings; local certificates retained, full D2 not promoted')


if __name__ == '__main__':
    {'export': export, 'build':build, 'seal':seal, 'archive-initial':archive_initial_client, 'check':check, 'report':report, 'verify':verify}[sys.argv[1]]()

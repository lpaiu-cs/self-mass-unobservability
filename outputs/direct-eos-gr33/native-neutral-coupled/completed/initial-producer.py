"""Counterexample candidate: close native neutral transport inside radiation.

Reuse the frozen175 histories and native tangent. Before any new evolution,
identify the actual missing neutral rate and test a sparse stage representation.
"""
from pathlib import Path
import json,resource,sys,time
import numpy as np
import sympy as sp
import return_native_pressure_reciprocity as previous

OUT=Path('native-neutral-coupled177-work');BEFORE=previous.OUT
read,write,sha=previous.read,previous.write,previous.sha
LD=np.longdouble;AMP=previous.AMP
CAPS=dict(prepare=10,localize=30)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failed=read(BEFORE/'block-result.json');assert not failed['passed']
    files=[Path(__file__),Path(previous.__file__),Path(previous.engine.__file__),Path(previous.engine.face.__file__),
        BEFORE/'block-result.json',BEFORE/'result.json']
    files.extend(p/f'steps-128-reference-128.npz' for s in [0,1] for p in previous.paths(s))
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='04bf81e3d',
        claim='Resolve the failed native noncollisional H feedback in the same photon/material equation before any GR charge readout.',
        failure='175 full histories pass their individual time/conservation tests but M_H and dM_H differ20.7percent and42.1percent; preserve the0.2percent rejection.',
        method='Reuse the native analytic conserved-to-primitive/HLL/donor tangent. On three saved states decompose the actual neutral flux difference and test additivity before selecting a sparse coupled-stage representation. A separately evolved gas H must not be overwritten to fake agreement.',
        decision='Only a representation that reproduces the native neutral flux, including branch effects, may enter an actual joint stage solve. Failure requires a source-level repair or explicitly nonlinear stage; no automatic further waveform sweep.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=3,new_photon_or_material_steps=0,
        scope='Source-level preparation for a real coupled solve; neither final charge progress by itself nor a uniform derivative certificate.',
        gates=dict(native_representation=.002,conservation=1e-12),
        bindings={str(p):sha(p) for p in files}))
    H,C,M,T=sp.symbols('H C M T',cls=sp.Function);t=sp.symbols('t')
    assert sp.simplify(sp.diff(H(t)-C(t),t).subs(sp.diff(H(t),t),sp.diff(C(t),t)+T(t))-T(t))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='For Hdot=Cdot+T_H, the noncollisional inventory M=H-C satisfies Mdot=T_H. This algebra does not establish a sparse native derivative or a coupled numerical solution.'))


def localize():
    previous.initialize();m=previous.run.c.Material(128,128);rows=[];states=[]
    for s in [0,1]:
        d=np.load(previous.paths(s)[1]/'steps-128-reference-128.npz')
        ids=[np.argmin(abs(d['t']-t)) for t in m.t];states.append(d['history_scaled'][ids].copy())
    norm=lambda a:np.sum(abs(a),dtype=LD)
    for k in [1,8,16]:
        old,new=states[0][k],states[1][k]
        def transport(z):return -np.diff(previous.engine.tangent(m,k,z,1.)[0][3])
        oldrate,newrate=transport(old),transport(new);difference=newrate-oldrate
        changes=[]
        for j in range(4):
            z=old.copy();z[j]=new[j];changes.append(transport(z)-oldrate)
        zero=np.zeros_like(new);parts=[]
        for j in range(4):
            z=zero.copy();z[j]=new[j];parts.append(transport(z))
        scale=max(norm(newrate),LD('1e-290'))
        rows.append(dict(k=k,t=float(m.t[k]),deep_cells=m.nb,largest_changed_rate_cell=int(np.argmax(abs(difference))),
            neutral_rate_L1_physical=float(AMP*norm(newrate)),difference_L1_physical=float(AMP*norm(difference)),
            changed_component_L1_physical=[float(AMP*norm(v)) for v in changes],component_order=['B','S','Eref','H'],
            difference_decomposition_relative=float(norm(sum(changes)-difference)/max(norm(difference),LD('1e-290'))),
            zero_based_additivity_relative=float(norm(sum(parts)-newrate)/scale),
            neutral_only_rate_L1_physical=float(AMP*norm(parts[3])),
            conservation_relative=float(abs(np.sum(newrate,dtype=LD))/scale)))
    write(OUT/'localization.json',dict(classification='Counterexample candidate',rows=rows,no_evolution_steps=True))
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))

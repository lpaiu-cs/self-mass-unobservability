"""Counterexample candidate: compact charge of reciprocal incident response."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,resource,sys,time
import numpy as np
import solve_native_incident_reciprocal as solve
import read_native_incident_response as prior

read,write,sha=solve.read,solve.write,solve.sha


def run(sweep):
    root=solve.OUT;out=root/f'sweep-{sweep}';photon,material=solve.paths(sweep);gr=out/'gr'
    residual=read(out/'residual.json')
    assert residual['passed'] and residual['finite_reciprocal_residual_passed']
    partition=[]
    for n in [64,128]:
        p=np.load(photon/f'steps-{n}-reference-128.npz')
        pp,mm=solve.paths(sweep-1)
        oldp=np.load(pp/f'steps-{n}-reference-128.npz');oldm=np.load(mm/f'steps-{n}-reference-128.npz')
        ids=[int(np.argmin(abs(oldm['t']-t))) for t in p['t']]
        a=p['moments'][:,[1,2]]-p['collision_transfer'].transpose(0,2,1)
        b=oldm['history_scaled'][ids][:,[2,3]]*solve.AMP-oldp['collision_transfer'].transpose(0,2,1)
        err=solve.relative(a,b);assert max(err)<1e-8;partition.append(err)
        prefix=np.load(photon/f'pilot-{n}.npz')
        for key in ['t','moments','material_history','collision_transfer','radial_ports','accepted_angular_times','accepted_angular_luminosity']:
            assert np.array_equal(p[key][:len(prefix[key])],prefix[key]),key
    assert not gr.exists();gr.mkdir()
    files=[Path(__file__),Path(solve.__file__),Path(prior.__file__),out/'residual.json',root/'plan.json',root/'linear-repair-plan.json']
    files += [p/f'steps-{n}-reference-128.npz' for p in [photon,material] for n in [64,128]]
    write(out/'readout-plan.json',dict(classification='Counterexample candidate',
        claim='Read the compact retarded scalar response after the returned matter/photon waveform passes its0.2percent residual. Reuse the direct primitive, canonical geometry subtraction and independent GR readout.',
        scope='Prescribed external metric and17-knot retained waveform equations. Residual is not an error bound or nonlinear/native/Einstein closure. Full direct incident scattering and exterior photons at infinity remain outside this readout.',
        gates=dict(time=.02,quadrature=.002,pressure=.002,independent=1e-9,identities=1e-12),
        seconds=180,bindings={str(p):sha(p) for p in files}))
    def initialize():
        solve.initialize(sweep);solve.base.Material=solve.Material
    namespace=dict(prior.compact.__globals__,OUT=out,PHOTON=photon,MATERIAL=material,GR=gr,
                   fixed=SimpleNamespace(initialize=initialize,write=write))
    FunctionType(prior.compact.__code__,namespace)()
    maximum=0.;LD=np.longdouble
    for n in [64,128]:
        d=np.load(gr/f'source-{n}-reference-128.npz');s=np.load(material/f'stress-{n}-reference-128.npz')['material']
        rest=d['baryon_g'].astype(LD)*LD(d['cx'])*LD(solve.base.C)**2;total=rest+d['gas_nonrest_energy_erg']
        errors=[total-s[:,0],d['nonrest_trace_erg']+rest-(s[:,0]-s[:,1]-2*s[:,3]),
                d['nonrest_stress_erg']+rest-(s[:,0]-s[:,1]),d['pressure_volume_erg']-s[:,3],
                d['metric_stress_erg']-(total+d['photon_energy_erg']-s[:,1]-d['photon_radial_pressure_erg'])]
        maximum=max(maximum,float(max(np.max(abs(e)) for e in errors)/max(np.max(abs(s)),1e-290)))
    assert maximum<1e-12
    result=read(gr/'result.json');old=read(solve.prior.OUT/'gr/result.json')
    result.update(finite_reciprocal_waveform_residual_passed=True,sweeps=sweep,residual=residual,mechanical_collision_partition_relative=partition,accepted_prefixes_preserved=True,
                  source_identity_relative=maximum,previous_one_way_endpoint=old['compact_return_endpoint'],
                  signed_endpoint_fraction_change=(result['compact_return_endpoint']-old['compact_return_endpoint'])/abs(old['compact_return_endpoint']),
                  complete_reciprocal_fixed_point=False,self_GR_returned=False,
                  scope='Accepted finite material-photon waveform residual under prescribed external GR; not full Einstein-matter fixed point, continuum error certificate or final infinity charge.')
    write(out/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    sweep=int(sys.argv[1]);resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    solve.base.drive.native.deadline(180);start=time.monotonic();cpu=time.process_time();error=None
    try:
        spent=sum(read(p)['seconds'] for p in solve.OUT.rglob('*-receipt.json'))
        assert spent+180<=read(solve.OUT/'plan.json')['budget']['total_action_seconds']
        run(sweep)
    except Exception as exc:error=repr(exc);raise
    finally:
        p=solve.OUT/f'compact-{sweep}-receipt.json';assert not p.exists()
        write(p,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
                     peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))

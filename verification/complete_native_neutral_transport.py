"""Counterexample candidate: return the actual177 prefix before a full run.

Use only completed photon inputs, with NaNs beyond their horizon. Keep the
native free-material equation; an independent H overwrite is forbidden.
"""
from pathlib import Path
import gc,json,resource,sys,time
import numpy as np
import couple_native_neutral_transport as coupled

previous=coupled.previous;OUT=coupled.OUT;LD=coupled.LD;AMP=coupled.AMP
read,write,sha=coupled.read,coupled.write,coupled.sha
CAP=60


def run(retry=False):
    assert read(OUT/'pilot-result.json')['passed']
    if retry:
        failure=read(OUT/'material-prefix-receipt.json')
        assert 'Future photon transfer' in failure['error'] and failure['seconds']+CAP<60
        write(OUT/'material-prefix-guard-repair.json',dict(classification='Counterexample candidate',failure=failure,
            repair='Place the finite-transfer guard at actual rhs consumption, not fields used for background-only activity masks. At the prefix endpoint use the existing closing-stage left derivative for the endpoint probe. Future transfers remain NaN; no future collision input is used.',
            source_sha256=sha(__file__),original_producer_sha256=sha(OUT/'material-prefix-pre-guard.py'),remaining_seconds=CAP))
    old=read(coupled.BEFORE/'sweep-1/material/pilot.json')
    write(OUT/'material-prefix-plan.json',dict(classification='Conjectural',checkpoint='9e32f1896',
        claim='Apply the actually completed177 photon prefix to the unchanged free-material equation. Test whether its noncollisional H inventory and native transport match the photon/H joint solve before spending on a full horizon.',
        reuse='No new photon, background or EOS trajectory. Construct the original17-knot transfer array with ONLY measured prefix values; all future entries are NaN. Reject any use beyond the prefix. Keep the original full-horizon macro clock and stop after4/8steps.',
        gates=dict(time=.02,conservation=1e-8,directional=.002,owner=1e-8,branch=.01,
                   paired=.002,mechanical_H=.002,native_H_rate=.002),
        comparison='Compare free H-C to the actual Radau integral of native H flux; compare native flux evaluated on the free material trajectory at every actual photon stage. These are the original0.2percent transport gates on the newly implicit H equation, not the obsolete lagged H input.',
        forecast=dict(previous_same_material_pair_seconds=old['seconds'],
            estimate_range_seconds=[old['seconds'],4*old['seconds']],
            assumption='Extra flux comparisons and setup may change timing;60s hard cap,3GiB,one CPU thread. No full horizon or refinement on failure.'),
        seconds=CAP,full_horizon_authorized=False,
        bindings={str(p):sha(p) for p in [Path(__file__),Path(coupled.__file__),OUT/'pilot-result.json',
            OUT/'pilot_retry-receipt.json']+[coupled.paths(1)[0]/f'pilot-{n}.npz' for n in [64,128]]}))
    previous.OUT=OUT;previous.paths=coupled.paths;previous.initialize();owner=previous.run.c
    original=owner.transfers
    def prefix_transfer(m,path,n):
        if Path(path)!=coupled.paths(1)[0]:return original(m,path,n)
        with np.load(path/f'pilot-{n}.npz') as p:
            count=len(p['t']);assert np.array_equal(p['t'],m.t[:count])
            c=p['collision_transfer'].astype(LD)
            m.transfer=np.full((len(m.t),4,m.n),np.nan,dtype=LD)
            m.transfer[:count]=np.stack([np.zeros_like(c[:,:,0]),p['moments'][:,3]/m.a,c[:,:,0],c[:,:,1]],axis=1)/LD(AMP)
            m.prefix_end=float(p['t'][-1])
    owner.transfers=prefix_transfer;rows=[];histories=[]
    for n in [64,128]:
        mark=time.monotonic();m=owner.Material(128,n);fields=m.fields;rhs=m.rhs
        def limited(t):
            assert 0<=t<=m.prefix_end+1e-18,('Uncomputed photon input',t,m.prefix_end)
            return fields(t)
        def consumed(t,z,probe=1.):
            j=limited(t)[0]
            assert np.isfinite(m.transfer[j:j+2]).all(),('Future photon transfer',j,t)
            return rhs(t,z,probe)
        m.fields=limited;m.rhs=consumed;label=f'prefix-{n}'
        row=m.run(n,label,n//16)
        path=coupled.paths(1)[1]/f'{label}.npz';d=dict(np.load(path))
        p=dict(np.load(coupled.paths(1)[0]/f'pilot-{n}.npz'))
        assert d['t'][-1]==m.prefix_end
        ids=[int(np.argmin(abs(d['t']-t))) for t in p['t']]
        assert np.max(abs(d['t'][ids]-p['t']))<1e-18
        history=d['history_scaled'][ids].copy();histories.append(history)
        balance=float(np.max(abs(np.sum(d['history_scaled'],axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        actual=m.transfer[:len(p['t'])]
        MH=history[-1,3]-actual[-1,3]
        joint=np.sum(p['native_neutral_stage_weights'][:,None].astype(LD)*p['native_neutral_stage_rates'],axis=0,dtype=LD)/LD(AMP)
        mismatch=float(np.sum(abs(MH-joint),dtype=LD)/max(np.sum(abs(MH),dtype=LD),LD('1e-290')))
        differences=[];norms=[]
        for t,rate in zip(p['native_neutral_stage_times'],p['native_neutral_stage_rates']):
            j=np.clip(np.searchsorted(d['t'],t,side='left')-1,0,len(d['t'])-2)
            w=(t-d['t'][j])/(d['t'][j+1]-d['t'][j]);z=(1-w)*d['history_scaled'][j]+w*d['history_scaled'][j+1]
            k,v,_,_=m.fields(t)
            F=sum(a*coupled.flux(m,i,z) for i,a in [(k,1-v),(k+1,v)] if a)
            native=-np.diff(F)*AMP
            differences.append(float(np.sum(abs(native-rate),dtype=LD)));norms.append(float(np.sum(abs(native),dtype=LD)))
        rate_error=max(differences)/max(max(norms),1e-290)
        m.closing_stage=True
        rates=[m.rhs(m.prefix_end,d['delta_scaled'],v)[0] for v in [.5,1.,2.]]
        m.closing_stage=False
        probe=[(np.sum(abs(v-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).astype(float).tolist() for v in [rates[0],rates[2]]]
        paired=owner.relative(p['moments'][:,[1,2]]/AMP,history[:,[2,3]])
        row.update(balance=balance,mechanical_H_relative=mismatch,native_H_rate_relative=rate_error,
            paired_E_H=paired,endpoint_probe_half_nominal_double=probe,physical_branch_ratio=m.physical_branch_ratio,
            maximum_owner_error=max(v['owner_error'] for v in m.cache.values()),worker_seconds=time.monotonic()-mark)
        row['physical_material_passed']=bool(row['passed'] and balance<1e-8 and np.max(probe)<.002 and row['maximum_owner_error']<1e-8 and m.physical_branch_ratio<.01)
        row['transport_matches_joint_H']=bool(max([mismatch,rate_error]+paired)<.002)
        rows.append(row);write(path.with_suffix('.json'),row);print(json.dumps(row),flush=True)
        np.savez_compressed(OUT/f'prefix-transport-{n}.npz',mechanical_H_free_scaled=MH,
            mechanical_H_joint_scaled=joint,stage_rate_difference_L1=differences,stage_rate_norm_L1=norms)
        assert row['physical_material_passed'],row
        del m;gc.collect()
    errors=owner.relative(*histories)
    result=dict(classification='Counterexample candidate',passed=max(errors)<.02 and all(v['transport_matches_joint_H'] for v in rows),
        rows=rows,time_comparison=errors,new_photon_steps=0,original_transport_gate=.002,
        full_horizon_authorized=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'material-prefix-result.json',result);print(json.dumps(result),flush=True)
    assert result['passed'],('Actual prefix material transport',result)


if __name__=='__main__':
    retry=len(sys.argv)>1 and sys.argv[1]=='retry'
    if retry:CAP=50
    receipt=OUT/('material-prefix-retry-receipt.json' if retry else 'material-prefix-receipt.json');assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));previous.original.inf.incident.native.deadline(CAP)
    start=time.monotonic();cpu=time.process_time();error=None
    try:run(retry)
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))

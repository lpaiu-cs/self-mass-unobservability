"""Preserve the failed double-precision defect; reuse packet owner stably."""
from pathlib import Path
from types import FunctionType
import inspect,json,signal,sys,textwrap,time
import numpy as np
import def_native_collision_defect_response as response

forcing=response.forcing; ORIGINAL=forcing.OUT; OUT=ORIGINAL/'paired-forcing'
write=forcing.write;sha=forcing.sha


def configure():
    response.configure();forcing.OUT=OUT


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(response.OUT/'pilot-failure.json',dict(classification='Counterexample candidate',
        error='Native source number pairing: knot0 relative1.835233536714013e-10 exceeds1e-12',
        failed_cell=134,deep_max=5.211841994585709e-18,physical_steps=0,
        exact_elapsed_unavailable=True,charged_pilot_seconds=90))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Repair cancellation in native-minus-table atmospheric scattering before applying the actual finite response. Preserve the same remapping, velocities, opacity and frequency exits.',
        method='At identical velocities apply the signed opacity-product difference through the existing linear scattering owner. For unequal velocities evaluate the same owner twice in extended precision and retain the tiny velocity difference. Preserve each old forcing file and reuse its bound-free/deep fields.',
        correction='The earlier j diagnosis was incorrect: runtime velocity owner is def_native_material_join.Coupled.velocity and reads Pi. Restoring j is harmless but is not a physical operator repair. The regenerated coefficient bank is retained; no missing-velocity claim is supported.',
        gates=dict(forcing_number=1e-12,arithmetic_change=.002,positive=.002,integrated_net=.02),
        budgets=dict(repair_seconds=60,retry_pilot_seconds=60,production_seconds=900,CPU_threads=1,virtual_GiB=3,new_native_calls=0),
        reallocation='Charge the failed preflight the full90s pilot cap. Use60s of72.946806533s unused bank allowance for one retry pilot. Repair uses60s from166.228125535s unused collision production allowance. Original aggregate budgets unchanged; completed bank reused.',
        stop='No threshold relaxation or further automatic retries; preserve a failure and assess the actual owner.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(response.__file__),Path(forcing.__file__),ORIGINAL/'production.json',response.OUT/'plan.json',response.OUT/'bank-result.json',response.OUT/'pilot-failure.json']}))


def repair():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(60)
    data=forcing.inputs();m,ds,ats,z,by=data;f=m.flow;nb=m.bulk.n
    owner=m.scattering.__func__;source=textwrap.dedent(inspect.getsource(owner))
    assert source.count('escape=np.zeros((n,3))')==1
    source=source.replace('escape=np.zeros((n,3))','escape=np.zeros((n,3),dtype=I.dtype)')
    source=source.replace('loss.sum((1,2)),1e-250','abs(loss).sum((1,2)),1e-250')
    scope=dict(owner.__globals__);exec(compile(source,__file__,'exec'),scope);scatter=scope['scattering']
    (OUT/'expanded-scattering.py').write_text(source)
    assert m.velocity.__func__.__module__=='def_native_material_join'
    rows=[];changed=0
    for k in range(17):
        d=ats[k];ids=d['active'];v=d['v'][ids];rho=d['rho'][ids]*f.eos.rho0
        f.eos.y=d['y'];kap=f.eos(d['rho'],d['lt'])[4][ids]
        roots=[by['atmosphere',k,int(j)] for j in ids]
        rn=np.array([r['rho'] for r in roots])*f.eos.rho0;vn=np.array([r['v'] for r in roots])
        kn=np.array([r['raw'][13]/r['raw'][0]*6.6524587321e-25/forcing.AMU for r in roots])
        I=z['snapshot_I'][k].sum(0)[ids];a=m.m.a[ids];vol=(4*np.pi*m.m.RJ**2*m.m.vol)[ids]
        original=dict(np.load(ORIGINAL/f'point-{k}.npz'));ph=original['photon'].copy();exit=original['escape'].copy()
        saved_error=m.scatter_number_error
        old,oe=m.scattering(I,rho,v,kap,a);new,ne=m.scattering(I,rn,vn,kn,a)
        ld=np.longdouble;same=vn==v;changed+=int((~same).sum())
        delta=np.zeros(I.shape,dtype=ld);ex=np.zeros((len(ids),3),dtype=ld)
        # All operations stay in the same packet owner; signed rates need an
        # absolute-loss diagnostic, never the positive-rate denominator.
        if same.any():
            product=rn[same].astype(ld)*kn[same]-rho[same].astype(ld)*kap[same]
            delta[same],ex[same]=scatter(m,I[same].astype(ld),np.ones(same.sum(),dtype=ld),v[same].astype(ld),product,a[same].astype(ld))
        if (~same).any():
            j=~same
            p,e=scatter(m,I[j].astype(ld),rn[j].astype(ld),vn[j].astype(ld),kn[j].astype(ld),a[j].astype(ld))
            p0,e0=scatter(m,I[j].astype(ld),rho[j].astype(ld),v[j].astype(ld),kap[j].astype(ld),a[j].astype(ld))
            delta[j]=p-p0;ex[j]=e-e0
        m.scatter_number_error=saved_error
        # Recompose from bound-free instead of subtracting rounded old packets.
        ph[nb+ids]=original['bound'][nb+ids]+delta
        exit[:,nb+ids]=(ex*vol[:,None]).T
        N=m.energy_weight[ids]/m.E
        residual=np.sum(delta*N,axis=(1,2))+ex[:,0]*vol
        norm=np.maximum(np.sum(abs(delta)*N,axis=(1,2))+abs(ex[:,0])*vol,1.)
        number=float(np.max(abs(residual)/norm))
        arithmetic=float(np.sum(abs(delta-(new-old))*m.energy_weight[ids])/max(np.sum(abs(delta)*m.energy_weight[ids]),1.))
        # The tiny velocity part can have a worse relative cancellation ratio;
        # the actual exported paired source is checked independently below.
        number_source=np.sum((ph-original['bound'])*np.r_[m.bulk.photon_energy_weight/m.bulk.d['Einf'],m.energy_weight/m.E],axis=(1,2))+exit[0]
        source_scale=np.maximum(np.sum((abs(ph)+abs(original['bound']))*np.r_[m.bulk.photon_energy_weight/m.bulk.d['Einf'],m.energy_weight/m.E],axis=(1,2))+abs(exit[0]),1.)
        paired=float(np.max(abs(number_source)/source_scale))
        assert arithmetic<.002 and paired<1e-12,(k,arithmetic,paired)
        dm=original['defect_moments'].copy();w=m.energy_weight[ids]
        dm[0,nb+ids]=np.sum(ph[nb+ids]*w,axis=(1,2))+exit[1,nb+ids]
        dm[2,nb+ids]=-(np.sum(ph[nb+ids]*w*m.mu[None,:,None],axis=(1,2))+exit[2,nb+ids])/a
        original.update(photon=ph,escape=exit,defect_moments=dm)
        np.savez_compressed(OUT/f'point-{k}.npz',**original)
        rows.append(dict(k=k,scattering_number=number,export_number=paired,arithmetic_relative=arithmetic,changed_velocity_cells=int((~same).sum())))
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,changed_velocity_cells=changed,
                seconds=time.monotonic()-start,new_native_calls=0,old_velocity_diagnosis_retracted=True)
    write(OUT/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def pilot():
    assert json.loads((OUT/'result.json').read_text())['passed'];configure()
    source=inspect.getsource(response.pilot).replace('signal.alarm(90)','signal.alarm(60)')
    scope=dict(vars(response));exec(compile(source,__file__,'exec'),scope);scope['pilot']()


def production():configure();response.production()


if __name__=='__main__':
    cap=3*1024**3;forcing.resource.setrlimit(forcing.resource.RLIMIT_AS,(cap,cap))
    signal.signal(signal.SIGALRM,forcing.history.flow.old.optical.timeout)
    action=sys.argv[1]
    try:globals()[action]()
    except Exception as exc:write(OUT/f'{action}-failure.json',dict(error=repr(exc)));raise

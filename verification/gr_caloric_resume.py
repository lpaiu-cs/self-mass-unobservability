"""Resume accepted caloric endpoints without changing the scientific algorithm."""
import ast,inspect,json,shutil,sys
from pathlib import Path
import numpy as np
import gr_caloric_increment as base
import gr_caloric_refinement as refinement

g=base.g;OUT=refinement.OUT;RECOVERY=OUT/'recovery-1'


def save(name,value): (RECOVERY/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def initialize_path(steps,n,duration,mass,initial):
    ld=np.longdouble;shift=np.zeros(n,dtype=ld);old_a=initial.copy()
    total_energy=np.zeros(n,dtype=ld);exchanged=ld(0);history=[];dt=duration/steps
    path=OUT/f'path-{steps}-progress.json'
    if path.exists():
        history=json.loads(path.read_text())['rows']
        assert [row['step'] for row in history]==list(range(len(history))) and len(history)<=steps
        for step,row in enumerate(history):
            data=dict(np.load(OUT/f'endpoint-{steps}-{step}.npz'))
            assert np.array_equal(data['old_shift'],shift),(steps,step,'shift continuity')
            U=data['caloric_increment'];L=data['interior_Linf'];a=data['eos']
            energy=mass.astype(ld)*U;target=ld(dt)*(np.r_[L,ld(0)]-np.r_[ld(0),L])
            residual=energy-target
            assert np.array_equal(residual,data['energy_residual']),(steps,step,'residual replay')
            norm=float(abs(residual/(mass.astype(ld)*a[:,10])).max());exchange=abs(target).sum()
            global_error=float(abs(residual.sum())/max(ld(1),exchange))
            assert norm==row['local_scaled_energy_residual'] and global_error==row['global_energy_relative_to_exchange']
            assert row['finite_quadrature_passed'] and row['caloric_chart_entropy_change_erg_K']>=0
            assert row['native_endpoint_difference_within_32ulp_cells']==n
            shift+=data['step_shift'].astype(ld);total_energy+=energy;exchanged+=exchange;old_a=a.copy()
    print('RESUME ACCEPTED CALORIC ENDPOINTS',steps,len(history),flush=True)
    return shift,old_a,total_energy,exchanged,history,dt,len(history)


def save_path(path,**arrays):
    if path.exists():
        old=dict(np.load(path));assert old.keys()==arrays.keys()
        assert all(np.array_equal(old[k],value) for k,value in arrays.items()),path
    else:np.savez_compressed(path,**arrays)


def transformed_source():
    source=inspect.getsource(base.run)
    replacements={
        'shift=np.zeros(n,dtype=ld);old_a=base.copy();total_energy=np.zeros(n,dtype=ld);exchanged=ld(0);history=[];dt=duration/steps':
        'shift,old_a,total_energy,exchanged,history,dt,start_step=initialize_path(steps,n,duration,mass,base)',
        'for step in range(steps):':'for step in range(start_step,steps):',
        "np.savez_compressed(OUT/f'path-{steps}.npz',":"save_path(OUT/f'path-{steps}.npz',"}
    changed=source
    for old,new in replacements.items():
        assert changed.count(old)==1,old;changed=changed.replace(old,new)
    reverse=changed
    for old,new in replacements.items():reverse=reverse.replace(new,old)
    assert reverse==source and ast.dump(ast.parse(reverse))==ast.dump(ast.parse(source))
    return source,changed,replacements


def prepare():
    assert not RECOVERY.exists();RECOVERY.mkdir()
    source,changed,replacements=transformed_source()
    (RECOVERY/'original-run.py').write_text(source)
    (RECOVERY/'resumed-run.py').write_text(changed)
    paths=[g.ROOT/'verification/gr_caloric_resume.py',g.ROOT/'verification/gr_caloric_increment.py',
        g.ROOT/'verification/gr_caloric_refinement.py',OUT/'plan.json',OUT/'path-8-progress.json',
        *sorted(OUT.glob('endpoint-8-*.npz'))]
    shutil.copy2(OUT/'path-8-progress.json',RECOVERY/'accepted-path-8-progress.json')
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f161e59',
        reason='Authoritative process listing found all prior scientific PIDs absent; /proc/uptime showed a newly restarted WSL instance. Four accepted path-8 endpoints were retained. No endpoint-8-4 was accepted.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        changes=replacements,algorithm='Restore shift, endpoint EOS, accumulated energy and exchange with the exact original ordering. Skip only accepted steps. Original Newton, EOS, opacity, Gauss, energy and time tolerances remain byte-for-byte unchanged after reversing the three explicit source overlays.',
        physical_EOS_certified=False,full_GR_evolution=False))
    state,aux=base.old.micro.inputs();mass=state['dm']*np.exp(state['nu'])
    duration=json.loads((g.OUT/'gr-nonlinear-thermal/duration.json').read_text())['coordinate_seconds']
    restored=initialize_path(8,len(mass),duration,mass,aux['eos']);assert restored[-1]==4
    np.savez_compressed(RECOVERY/'restored-state.npz',lnT_shift=restored[0],eos=restored[1],
        total_energy_increment=restored[2],exchange=restored[3])
    save('control.json',dict(classification='Proven',passed=True,accepted_steps=4,
        stored_shift_chain_bitwise=True,stored_residual_replay_bitwise=True,algorithm_overlay_reversible=True,
        scope='Saved endpoint arithmetic and exact source preservation only; no new EOS/time result.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in RECOVERY.iterdir() if p.is_file()}))
    print('PASS accepted endpoint reconstruction and unchanged solver overlay',flush=True)


def run():
    plan=json.loads((RECOVERY/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((RECOVERY/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    source,changed,_=transformed_source();assert changed==(RECOVERY/'resumed-run.py').read_text()
    base.OUT=OUT
    namespace=dict(base.run.__globals__,OUT=OUT,initialize_path=initialize_path,save_path=save_path)
    exec(compile(changed,str(RECOVERY/'resumed-run.py'),'exec'),namespace)
    namespace['run']()
    for steps in [8,16,32]:refinement.audit(steps)


if __name__=='__main__':globals()[sys.argv[1]]()

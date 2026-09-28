"""Preserve native excitation scaling before diagnosing tiny derivatives."""
import ctypes, shutil, sys
import numpy as np
import excitation_curvature as e

ORIGINAL_OUT=e.OUT; ORIGINAL_CACHE=e.CACHE; ORIGINAL_NAME=e.NAME
e.OUT=e.g.OUT/'excitation-scaled'; e.CACHE=e.g.CACHE/'excitation-scaled'
e.NAME='free_eos_direct24_excitation_scaled'
e.LIB=e.CACHE/'build/src'/('lib'+e.NAME+'.so.1.0.0')
OUT=e.OUT; g=e.g; p=e.p; Base=e.EOS


def prepare():
    assert not OUT.exists() and not e.CACHE.exists(); OUT.mkdir(); e.CACHE.mkdir()
    source=e.CACHE/'source'; shutil.copytree(ORIGINAL_CACHE/'source',source); changes={}
    path=source/'src/mod_excitation_block.f90'; text=path.read_text()
    text=text.replace('  integer, save, public :: extrace_count',
        '  real(fp_kind), save, public :: extrace_scale(318)=1._fp_kind\n  integer, save, public :: extrace_count',1)
    path.write_text(text)
    path=source/'src/excitation_sum.f90'; text=path.read_text()
    text=text.replace('extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht',
        'extrace_value, extrace_grad, extrace_hess, extrace_hf, extrace_ht, extrace_scale')
    head,tail=text.split('  extrace_count=extrace_count+1;',1)
    trace,tail=tail.split('  nuvarsum_scale = nuvar/qratio_scale',1)
    assert trace.count('/qratio_scale')==7
    trace=trace.replace('/qratio_scale','').replace('  extrace_ids(:,trace_slot)=extrace_tag',
        '  extrace_ids(:,trace_slot)=extrace_tag\n  extrace_scale(trace_slot)=qratio_scale')
    path.write_text(head+'  extrace_count=extrace_count+1;'+trace+'  nuvarsum_scale = nuvar/qratio_scale'+tail)
    path=source/'src/CMakeLists.txt'; path.write_text(path.read_text().replace('OUTPUT_NAME '+ORIGINAL_NAME,'OUTPUT_NAME '+e.NAME))
    for name in ['mod_excitation_block.f90','excitation_sum.f90','CMakeLists.txt']:
        shutil.copy2(source/'src'/name,OUT/name)
        changes[name]=dict(before=g.c.sha(ORIGINAL_CACHE/'source/src'/name),after=g.c.sha(source/'src'/name))
    shutil.copy2(ORIGINAL_OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    plan=e.read(ORIGINAL_OUT/'plan.json')
    plan.update(checkpoint='f8066a3',source_changes=changes,control_cells=[0,1607,1612,1907,1916,2972,5734],
        method='Keep native exp(600)-scaled component values and all derivatives plus the actual scaling constant. Convert in extended exponent range before dimensionless products. Physics and original validation thresholds are unchanged. Preserve the failed original trace.',
        bindings={rel:g.c.sha(g.ROOT/rel) for rel in [
            'verification/excitation_scaled.py','verification/excitation_curvature.py',
            'outputs/direct-eos-gr33/excitation-curvature/failure-manifest.json',
            'outputs/direct-eos-gr33/pressure-ionization/manifest.json','outputs/direct-eos-gr33/reference-state.npz']})
    e.save('plan.json',plan)


class EOS(Base):
    def trace(self):
        result=super().trace()
        scale=self.array('mod_excitation_block','extrace_scale',318)[:result['count']]
        raw_scale=np.ones(5); raw_scale[:len(scale)]=scale
        assert np.all(np.isfinite(raw_scale)) and np.all(raw_scale>0)
        result['raw_scale']=raw_scale
        for name in ['value','grad','hess','hf','ht']:
            raw=result[name]; result['raw_'+name]=raw.copy()
            converted=raw.astype(np.longdouble); divisor=raw_scale.astype(np.longdouble)
            if name=='value': converted[:,3:]/=divisor[:,None]
            elif name=='hess': converted/=divisor[:,None,None]
            else: converted/=divisor[:,None]
            result[name]=converted
        return result


def control():
    plan=e.read(OUT/'plan.json'); state=dict(np.load(g.OUT/'reference-state.npz')); eos=EOS(); rows=[]; e.symbolic()
    for i in plan['control_cells']:
        start=i//128*128; j=i-start
        electron=dict(np.load(p.ex.OUT/f'block-{start}.npz'))['values'][j]
        base,trace,G,H,Hnon,row=eos.full_excitation(state['lnd'][i],state['lnT'][i],state['X'][i],electron)
        snap,a,meta,extra,rn,ri,G0,H0,computed,native,_=base
        old=dict(np.load(p.OUT/f'block-{start}.npz'))
        for k,v in dict(state=meta,extra=extra,neutral=rn,ion3=ri,Gram=G0,Hessian=H0,computed=computed,native=native).items():
            assert np.array_equal(v,old[k][j]),(i,k)
        row.update(cell=i,passed=e.gates(row,plan)); rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**trace,Gram=G,Hessian=H,nonideal_Hessian=Hnon)
        print('SCALED EXCITATION CONTROL',row,flush=True)
    e.save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows))
    assert all(r['passed'] for r in rows)


e.EOS=EOS
build=e.build; run=e.run
if __name__=='__main__': globals()[sys.argv[1]]()

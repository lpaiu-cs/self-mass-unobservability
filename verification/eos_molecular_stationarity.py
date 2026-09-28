"""Trace molecular reaction constants and their free-energy gradients."""
import ctypes,shutil,sys
import numpy as np
import eos_atomic_stationarity_defined as defined

a=defined.a;e=a.e;p=a.p;g=a.g;Base=e.EOS;atomic_gates=e.gates
PREVIOUS_OUT=e.OUT;PREVIOUS_CACHE=e.CACHE;PREVIOUS_NAME=e.NAME
e.OUT=g.OUT/'molecular-stationarity';e.CACHE=g.CACHE/'molecular-stationarity'
e.NAME='free_eos_direct24_molecular_stationarity'
e.LIB=e.CACHE/'build/src'/('lib'+e.NAME+'.so.1.0.0');OUT=e.OUT


def prepare():
    assert not OUT.exists() and not e.CACHE.exists();OUT.mkdir();e.CACHE.mkdir()
    source=e.CACHE/'source';shutil.copytree(PREVIOUS_CACHE/'source',source);changes={}
    replacements={
        'mod_nuvar.f90':[('  implicit none','''  implicit none
  real(fp_kind), save, public :: mol_eq(2)=0._fp_kind, mol_dv(2)=0._fp_kind
  real(fp_kind), save, public :: mol_pl(2)=0._fp_kind, mol_linear(2)=0._fp_kind
  real(fp_kind), save, public :: mol_logs(3)=0._fp_kind
  integer, save, public :: mol_flags(4)=0''')],
        'eos_calc.f90':[
            ('  use mod_nuvar, only: nuvar','''  use mod_nuvar, only: nuvar, mol_eq, mol_dv, mol_pl, mol_linear, mol_logs, mol_flags'''),
            ('  nions = nionsp2 - 2','''  nions = nionsp2 - 2
  mol_eq=0._fp_kind;mol_dv=0._fp_kind;mol_pl=0._fp_kind
  mol_linear=0._fp_kind;mol_logs=0._fp_kind
  mol_flags=[ifh2,ifh2plus,ifnr,0]'''),
            ('     h2equil = logalpha - log(4._fp_kind) + qh2 - 1.5_fp_kind*tl + h2diss*tc2 + h2_log + dv(nions+1)',
             '''     h2equil = logalpha - log(4._fp_kind) + qh2 - 1.5_fp_kind*tl + h2diss*tc2 + h2_log + dv(nions+1)
     mol_flags(4)=1;mol_eq(1)=h2equil;mol_dv(1)=dv(nions+1)
     mol_linear(1)=logalpha-log(4._fp_kind)+qh2-1.5_fp_kind*tl+h2diss*tc2+h2_log
     if(ifpl.eq.1) mol_pl(1)=plop(nions+1)'''),
            ('        h2plusequil = h2equil + qh2plus - qh2 - bi(nions+1)*tc2 + dv(nions+2)',
             '''        h2plusequil = h2equil + qh2plus - qh2 - bi(nions+1)*tc2 + dv(nions+2)
        mol_eq(2)=h2plusequil;mol_dv(2)=dv(nions+2)
        mol_linear(2)=mol_linear(1)+qh2plus-qh2-bi(nions+1)*tc2
        if(ifpl.eq.1) mol_pl(2)=plop(nions+2)'''),
            ('     ! Transform nh2, nh2plus, nh_neutral, and nh_ion from ln (n/avogadro) form to n form',
             '''     mol_logs(1)=nh_neutral;mol_logs(2)=nh2
     if(ifh2plus.gt.0) mol_logs(3)=nh2plus
     ! Transform nh2, nh2plus, nh_neutral, and nh_ion from ln (n/avogadro) form to n form''')],
        'CMakeLists.txt':[('OUTPUT_NAME '+PREVIOUS_NAME,'OUTPUT_NAME '+e.NAME)]}
    for name,pairs in replacements.items():
        path=source/'src'/name;text=path.read_text()
        for old,new in pairs:
            assert text.count(old)==1,(name,old);text=text.replace(old,new)
        path.write_text(text);shutil.copy2(path,OUT/name)
        changes[name]=dict(before=g.c.sha(PREVIOUS_CACHE/'source/src'/name),after=g.c.sha(path))
    shutil.copy2(PREVIOUS_OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    molecules=np.concatenate([dict(np.load(p.c.s.OUT/f'block-{start}.npz'))['molecular_H_fractions'] for start in range(0,5735,128)])
    controls=sorted(set([0,2972,5734,*map(int,np.argmax(molecules,axis=0))]))
    plan=e.read(PREVIOUS_OUT/'plan.json')
    plan.update(checkpoint='092a341',source_changes=changes,control_cells=controls,
        molecular_gradient_absolute_tolerance=1e-8,molecular_log_population_tolerance=1e-8,
        molecular_linear_replay_tolerance=1e-8,
        method='At the actual final eos_calc call record H2 and H2+ constants, their explicit temperature/partition linear pieces, both native nonideal dv values, source-defined molecular PL terms and pre-underflow log(n/N_A). Compare with independent molecular free-energy gradients and positive saved populations. Native dv for H2+ is relative to H2, so the formation constant relative to 2H uses dv(H2)+dv(H2+).',
        scope='Frozen molecular reactions and existing atomic support only. Continuous primitive/partition/branch/root errors and physical EOS remain separate.',
        bindings={rel:g.c.sha(g.ROOT/rel) for rel in [
            'verification/eos_molecular_stationarity.py','verification/eos_atomic_stationarity.py',
            'verification/eos_atomic_stationarity_defined.py','verification/excitation_scaled.py','verification/excitation_curvature.py',
            'outputs/direct-eos-gr33/atomic-stationarity-defined/manifest.json','outputs/direct-eos-gr33/reference-state.npz']})
    e.save('plan.json',plan)


class EOS(Base):
    def full_excitation(self,r,t,x,electron):
        base,trace,G,H,Hnon,row=super().full_excitation(r,t,x,electron)
        snap,_,meta,*_=base;prefix='__mod_nuvar_MOD_'
        flags=np.ctypeslib.as_array((ctypes.c_int*4).in_dll(self.inventory_lib,prefix+'mol_flags')).copy()
        assert flags[2]==0 and flags[3] in [0,1]
        mol={name:self.array('mod_nuvar','mol_'+name,size) for name,size in [('eq',2),('dv',2),('pl',2),('linear',2),('logs',3)]}
        assert all(np.all(np.isfinite(v)) for v in mol.values())
        mu=trace['station_chemical_potentials'].astype(np.longdouble).copy();mu[-2:]-=mol['pl']
        predicted=np.array([2*mu[0]-mu[-2],mu[-2]-mu[-1]],dtype=float)
        ym=(x/g.c.A)@self.mapping;eps=ym/float(ym@self.weights)
        molecular=eps[0]*snap['molecular_H_fractions']/2;neutral=snap['number_fractions'][0,0]
        gradient_error=0.;population_error=0.;linear_error=0.;directions=0
        formation=np.cumsum(mol['dv']);chemical=np.cumsum(predicted)
        if flags[3]:
            linear_error=float(abs(mol['eq'][0]-mol['linear'][0]-formation[0]))
            if flags[1]>0:linear_error=max(linear_error,float(abs(mol['eq'][1]-mol['linear'][1]-formation[1])))
            for j in range(2):
                if j==1 and flags[1]<=0:continue
                if neutral>0 and molecular[j]>0:
                    logratio=float(np.log(molecular[j])-2*np.log(neutral)-np.log(meta[1]))
                    gradient_error=max(gradient_error,abs(chemical[j]-formation[j]))
                    population_error=max(population_error,abs(logratio-mol['eq'][j]));directions+=1
        else:assert np.all(molecular==0)
        trace.update(**{'mol_'+name:v for name,v in mol.items()},mol_flags=flags,mol_predicted=predicted,
            mol_populations=molecular,mol_neutral=np.float64(neutral))
        row.update(molecular_gradient_absolute_error=float(gradient_error),molecular_log_population_error=float(population_error),
            molecular_linear_replay_error=float(linear_error),positive_molecular_directions=directions)
        return base,trace,G,H,Hnon,row


def gates(row,plan):
    return atomic_gates(row,plan) and row['molecular_gradient_absolute_error']<plan['molecular_gradient_absolute_tolerance'] and row['molecular_log_population_error']<plan['molecular_log_population_tolerance'] and row['molecular_linear_replay_error']<plan['molecular_linear_replay_tolerance']


def control():
    plan=e.read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[]
    for i in plan['control_cells']:
        start=i//128*128;j=i-start;electron=dict(np.load(p.ex.OUT/f'block-{start}.npz'))['values'][j]
        base,trace,G,H,Hnon,row=eos.full_excitation(state['lnd'][i],state['lnT'][i],state['X'][i],electron)
        old=dict(np.load(PREVIOUS_OUT/f'block-{start}.npz'));current=dict(**trace,Gram=G,Hessian=H,nonideal_Hessian=Hnon)
        for k,v in old.items():assert np.array_equal(v[j],current[k]),(i,k)
        assert np.array_equal(base[0]['eos'],dict(np.load(p.c.s.OUT/f'block-{start}.npz'))['eos'][j])
        row.update(cell=i,passed=gates(row,plan));rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**current);print('MOLECULAR STATIONARITY CONTROL',row,flush=True)
    e.save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows))
    assert sum(r['positive_molecular_directions'] for r in rows)>0 and all(r['passed'] for r in rows)


def run():
    e.run();rows=[]
    for path in OUT.glob('block-*.json'):
        record=e.read(path);rows+=record['rows'];new=dict(np.load(path.with_suffix('.npz')))
        old=dict(np.load(PREVIOUS_OUT/path.with_suffix('.npz').name))
        for k,v in old.items():assert np.array_equal(v,new[k]),(record['start'],k)
    result=e.read(OUT/'result.json');result.update(
        molecular_gradient_absolute_error=max(r['molecular_gradient_absolute_error'] for r in rows),
        molecular_log_population_error=max(r['molecular_log_population_error'] for r in rows),
        molecular_linear_replay_error=max(r['molecular_linear_replay_error'] for r in rows),
        positive_molecular_directions=sum(r['positive_molecular_directions'] for r in rows),
        previous_atomic_trace_bitwise=True,continuous_stationarity_certified=False)
    e.save('result.json',result);print('MOLECULAR STATIONARITY COMPLETE',result,flush=True)


e.EOS=EOS;e.gates=gates;build=e.build
if __name__=='__main__':globals()[sys.argv[1]]()

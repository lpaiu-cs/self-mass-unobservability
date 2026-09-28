"""Compare final atomic equilibrium constants with free-energy gradients."""
import ctypes, shutil, sys
import numpy as np
import excitation_scaled as scaled

e=scaled.e; p=e.p; g=e.g; Base=e.EOS; original_gates=e.gates
PREVIOUS_OUT=e.OUT; PREVIOUS_CACHE=e.CACHE; PREVIOUS_NAME=e.NAME
e.OUT=g.OUT/'atomic-stationarity'; e.CACHE=g.CACHE/'atomic-stationarity'
e.NAME='free_eos_direct24_atomic_stationarity'
e.LIB=e.CACHE/'build/src'/('lib'+e.NAME+'.so.1.0.0'); OUT=e.OUT


def prepare():
    assert not OUT.exists() and not e.CACHE.exists(); OUT.mkdir(); e.CACHE.mkdir()
    source=e.CACHE/'source'; shutil.copytree(PREVIOUS_CACHE/'source',source); changes={}
    replacements={
        'mod_nuvar.f90':[('  implicit none','''  implicit none
  real(fp_kind), save, public :: station_dv(316)=0._fp_kind, station_zero(24)=0._fp_kind
  real(fp_kind), save, public :: station_ce(316)=0._fp_kind, station_binding(316)=0._fp_kind
  real(fp_kind), save, public :: station_plop(318)=0._fp_kind, station_tc2=0._fp_kind
  integer, save, public :: station_flags(3)=0, station_active(24)=0''')],
        'ionize.f90':[
            ('  use mod_nuvar, only:', '''  use mod_nuvar, only: station_dv, station_zero, station_ce, station_binding, &
       station_plop, station_tc2, station_flags, station_active
  use mod_nuvar, only:'''),
            ('  n_partial_elements = size(partial_elements) - 2','''  station_dv=0._fp_kind;station_zero=0._fp_kind
  station_ce=0._fp_kind;station_binding=0._fp_kind;station_active=0
  station_plop=0._fp_kind
  if(ifpl.eq.1) station_plop=plop
  station_tc2=tc2;station_flags=[ifnr,ifpl,ifionized]
  n_partial_elements = size(partial_elements) - 2'''),
            ('        ln_fractmax = -dvzero(ielement)','''        station_active(ielement)=1;station_zero(ielement)=dvzero(ielement)
        station_dv(ion0+1:ion0+iatomic_number(ielement))=dv(ion0+1:ion0+iatomic_number(ielement))
        station_ce(ion0+1:ion0+iatomic_number(ielement))=ce(ion0+1:ion0+iatomic_number(ielement))
        station_binding(ion0+1:ion0+iatomic_number(ielement))=bi_ref(ion0+1:ion0+iatomic_number(ielement))
        ln_fractmax = -dvzero(ielement)''')],
        'CMakeLists.txt':[('OUTPUT_NAME '+PREVIOUS_NAME,'OUTPUT_NAME '+e.NAME)]}
    for name,pairs in replacements.items():
        path=source/'src'/name; text=path.read_text()
        for old,new in pairs:
            assert text.count(old)==1,(name,old); text=text.replace(old,new)
        path.write_text(text); shutil.copy2(path,OUT/name)
        changes[name]=dict(before=g.c.sha(PREVIOUS_CACHE/'source/src'/name),after=g.c.sha(path))
    shutil.copy2(PREVIOUS_OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    plan=e.read(PREVIOUS_OUT/'plan.json'); plan.update(checkpoint='b22682b',source_changes=changes,
        control_cells=[0,2972,5734],atomic_gradient_absolute_tolerance=1e-8,
        atomic_log_population_tolerance=1e-8,
        method='Capture actual final atomic dv, reference shifts, statistical/binding constants, and optional temperature-only PL weights. Reconstruct chemical-potential differences independently from electron/exchange, Coulomb, MDH and scaled excitation gradients at the saved final state. Compare final positive atomic log ratios and largest-population reference differences; do not require a numerically absent neutral atom as the residual reference.',
        scope='Finite final-state atomic redistribution on current positive support. Molecular dissociation/ionization equilibrium, continuous root enclosures and physical EOS error are separate.',
        bindings={rel:g.c.sha(g.ROOT/rel) for rel in [
            'verification/eos_atomic_stationarity.py','verification/excitation_scaled.py','verification/excitation_curvature.py',
            'outputs/direct-eos-gr33/excitation-scaled/manifest.json','outputs/direct-eos-gr33/reference-state.npz']})
    e.save('plan.json',plan)


class EOS(Base):
    def full_excitation(self,r,t,x,electron):
        base,trace,G,H,Hnon,row=super().full_excitation(r,t,x,electron)
        snap,coulomb,meta,extra,rn,ri,*_=base
        module='mod_nuvar'; prefix='__'+module+'_MOD_'
        flags=np.ctypeslib.as_array((ctypes.c_int*3).in_dll(self.inventory_lib,prefix+'station_flags')).copy()
        active=np.ctypeslib.as_array((ctypes.c_int*24).in_dll(self.inventory_lib,prefix+'station_active')).copy()
        assert flags[0]==0 and flags[1] in [0,1]
        raw={name:self.array(module,'station_'+name,size) for name,size in [('dv',316),('zero',24),('ce',316),('binding',316),('plop',318)]}
        tc2=ctypes.c_double.in_dll(self.inventory_lib,prefix+'station_tc2').value
        assert all(np.all(np.isfinite(v)) for v in raw.values()) and tc2>0
        # Include zero-population reference species for gradient differences only.
        allpop=np.zeros((24,29))
        for el,Z in enumerate(g.d.CHARGES): allpop[el,:Z+1]=1.
        _,B,_,_=p.species_coordinates(allpop,np.ones(2),rn,ri)
        keys=[(el+1,z) for el,Z in enumerate(g.d.CHARGES) for z in range(Z+1)]+[(25,0),(26,1)]
        S=meta[1]*self.constants[0]; units=np.r_[p.L**np.arange(7),1.,p.L**3]
        u=extra[0]/units; kappa=S*(4*np.pi/3)*p.L**3
        _,_,grad2,grad3,_,_=[np.asarray(fun(*u),float) for fun in self.poly]
        gradient=kappa*grad2.ravel()+meta[11]*kappa**2*grad3.ravel()
        mu=np.asarray(B[:3].T@(coulomb[1:4]/coulomb[15])+B[3:].T@gradient,np.longdouble)
        mu+=B[0]*float(electron[10]+electron[1])
        factor=S*np.array([1.,p.L,p.L*p.L,1.]); moment_gradient=np.zeros(4,dtype=np.longdouble)
        L=np.zeros(len(keys),dtype=np.longdouble)
        for j in range(trace['count']):
            L[keys.index(tuple(trace['ids'][j]))]=trace['value'][j,3]
            moment_gradient+=trace['value'][j,0]*factor*trace['grad'][j]
        mu-=L+B[[3,4,5,10]].T@moment_gradient
        offset=0
        for el,Z in enumerate(g.d.CHARGES):
            for z in range(Z): mu[keys.index((el+1,z))]-=raw['plop'][offset+z]
            offset+=Z
        expected=np.zeros(316); actual=np.zeros(316); chemical_error=0.; population_error=0.; comparisons=0;offset=0
        for el,Z in enumerate(g.d.CHARGES):
            if active[el]:
                neutral=keys.index((el+1,0)); absolute=np.r_[0.,raw['dv'][offset:offset+Z]+raw['zero'][el]]
                candidate=np.array([float(mu[neutral]-mu[keys.index((el+1,z))]) for z in range(Z+1)])
                expected[offset:offset+Z]=candidate[1:]; actual[offset:offset+Z]=absolute[1:]
                positive=np.flatnonzero(snap['number_fractions'][el,:Z+1]>0)
                reference=int(positive[np.argmax(snap['number_fractions'][el,positive])])
                logw=np.r_[0.,raw['ce'][offset:offset+Z]-tc2*raw['binding'][offset:offset+Z]]+absolute
                for z in positive:
                    if z==reference: continue
                    chemical_error=max(chemical_error,abs((candidate[z]-candidate[reference])-(absolute[z]-absolute[reference])))
                    ratio=float(np.log(snap['number_fractions'][el,z])-np.log(snap['number_fractions'][el,reference]))
                    population_error=max(population_error,abs(ratio-(logw[z]-logw[reference])))
                    comparisons+=1
            offset+=Z
        trace.update(**{'station_'+name:value for name,value in raw.items()},station_tc2=np.float64(tc2),
            station_flags=flags,station_active=active,station_expected=expected,station_actual=actual,
            station_chemical_potentials=np.asarray(mu,float))
        row.update(atomic_gradient_absolute_error=float(chemical_error),atomic_log_population_error=float(population_error),
            positive_atomic_comparisons=comparisons)
        return base,trace,G,H,Hnon,row


def gates(row,plan):
    return original_gates(row,plan) and row['atomic_gradient_absolute_error']<plan['atomic_gradient_absolute_tolerance'] and row['atomic_log_population_error']<plan['atomic_log_population_tolerance']


def control():
    plan=e.read(OUT/'plan.json'); state=dict(np.load(g.OUT/'reference-state.npz')); eos=EOS(); rows=[]
    for i in plan['control_cells']:
        start=i//128*128; j=i-start; electron=dict(np.load(p.ex.OUT/f'block-{start}.npz'))['values'][j]
        base,trace,G,H,Hnon,row=eos.full_excitation(state['lnd'][i],state['lnT'][i],state['X'][i],electron)
        previous=dict(np.load(PREVIOUS_OUT/f'control-{i}.npz'))
        current=dict(**trace,Gram=G,Hessian=H,nonideal_Hessian=Hnon)
        for k,v in previous.items(): assert np.array_equal(v,current[k]),(i,k)
        row.update(cell=i,passed=gates(row,plan)); rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**current)
        print('ATOMIC STATIONARITY CONTROL',row,flush=True)
    e.save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows))
    assert all(r['passed'] for r in rows)


def run():
    e.run(); rows=[]
    for path in OUT.glob('block-*.json'):
        record=e.read(path); rows+=record['rows']; start=record['start']
        new=dict(np.load(path.with_suffix('.npz'))); old=dict(np.load(PREVIOUS_OUT/path.with_suffix('.npz').name))
        for k,v in old.items(): assert np.array_equal(v,new[k]),(start,k)
    result=e.read(OUT/'result.json'); result.update(atomic_gradient_absolute_error=max(r['atomic_gradient_absolute_error'] for r in rows),
        atomic_log_population_error=max(r['atomic_log_population_error'] for r in rows),
        positive_atomic_comparisons=sum(r['positive_atomic_comparisons'] for r in rows),
        previous_excitation_trace_bitwise=True,molecular_stationarity_certified=False,continuous_stationarity_certified=False)
    e.save('result.json',result); print('ATOMIC STATIONARITY COMPLETE',result,flush=True)


e.EOS=EOS; e.gates=gates; build=e.build
if __name__=='__main__': globals()[sys.argv[1]]()

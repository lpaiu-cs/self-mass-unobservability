"""Independent rational gradients and interval residuals on frozen atoms."""
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
from mpmath import iv
import eos_atomic_stationarity_defined as defined
import verify_pressure_ionization as pi
from interval_records import exact_endpoint

a=defined.a;e=a.e;p=a.p;g=a.g;OUT=e.OUT;rat=pi.rat;read=e.read


def interval(value): return iv.mpf(value.numerator)/value.denominator
def upper(value): return max(map(abs,map(exact_endpoint,value._mpi_)))
def upper_text(value):
    with localcontext() as context:
        context.prec=40;context.rounding=ROUND_CEILING
        text=str(Decimal(value.numerator)/Decimal(value.denominator))
    assert F(text)>=value
    return text


def run():
    assert not (OUT/'audit.json').exists();iv.dps=60
    plan=read(OUT/'plan.json');result=read(OUT/'result.json');bridge=read(OUT/'bridge.json')
    assert result['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(e.LIB)==bridge['library_sha256'] and g.c.sha(e.CACHE/'excitation.so')==bridge['bridge_sha256']
    assert g.c.sha(p.c.s.LIB)==plan['original_library_sha256']
    state=dict(np.load(g.OUT/'reference-state.npz'));eos=a.EOS();bounds=[]
    max_mu=F(0);max_equil=F(0);count=0;comparisons=0;cells_without_directions=0
    units=np.r_[p.L**np.arange(7),1.,p.L**3]
    for start in range(0,plan['cells'],128):
        record=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert record['start']==count and g.c.sha(path)==record['output_sha256']
        assert g.c.sha(OUT/'plan.json')==record['plan_sha256'] and all(a.gates(r,plan) for r in record['rows'])
        data=dict(np.load(path));previous=[]
        for folder in [p.OUT,p.c.OUT,p.ex.OUT,p.c.s.OUT,a.PREVIOUS_OUT]:
            path=folder/f'block-{start}.npz';assert g.c.sha(path)==read(folder/f'block-{start}.json')['output_sha256']
            previous.append(dict(np.load(path)))
        mdh,coulomb,electron,species,excitation=previous
        assert all(np.array_equal(data[k],v) for k,v in excitation.items())
        for cell in plan['control_cells']:
            if start<=cell<record['stop']:
                control=dict(np.load(OUT/f'control-{cell}.npz'))
                assert all(np.array_equal(data[k][cell-start],v) for k,v in control.items())
        for j in range(record['stop']-start):
            i=start+j
            # Fresh native output replay is separate from the gradient arithmetic.
            snap=eos.snapshot(state['lnd'][i],state['lnT'][i],state['X'][i])
            assert np.array_equal(snap['eos'],species['eos'][j]),i
            assert np.array_equal(snap['number_fractions'],species['number_fractions'][j]),i
            pop=species['number_fractions'][j];active=data['station_active'][j]
            total_pairs=sum(max(0,np.count_nonzero(pop[el,:Z+1])-1) for el,Z in enumerate(g.d.CHARGES) if active[el])
            assert total_pairs==record['rows'][j]['positive_atomic_comparisons']
            if total_pairs==0:
                bounds.append(dict(cell=i,comparisons=0,maximum_atomic_residual_upper='0',Fisher_squared_residual_upper='0'))
                cells_without_directions+=1;continue
            meta=mdh['state'][j];S=float(meta[1]*eos.constants[0]);k=rat(float(S*(4*np.pi/3)*p.L**3));q=rat(meta[11])
            _,g2,_=pi.monomial(pi.P2,list(map(rat,mdh['extra'][j,0]/units)))
            _,g3,_=pi.monomial(pi.P3,list(map(rat,mdh['extra'][j,0]/units)))
            gradient=[k*v+q*k*k*w for v,w in zip(g2,g3)]
            allpop=np.zeros((24,29))
            for el,Z in enumerate(g.d.CHARGES):allpop[el,:Z+1]=1.
            _,B,_,_=p.species_coordinates(allpop,np.ones(2),mdh['neutral'][j],mdh['ion3'][j])
            keys=[(el+1,z) for el,Z in enumerate(g.d.CHARGES) for z in range(Z+1)]+[(25,0),(26,1)]
            B=[[rat(v) for v in row] for row in B];factors=list(map(rat,S*np.array([1.,p.L,p.L*p.L,1.])))
            L={};gx=[F(0)]*4
            for component in range(int(data['count'][j])):
                scale=rat(data['raw_scale'][j,component]);nu=rat(data['raw_value'][j,component,0])
                tag=tuple(data['ids'][j,component]);L[tag]=rat(data['raw_value'][j,component,3])/scale
                for z in range(4):gx[z]+=nu*factors[z]*rat(data['raw_grad'][j,component,z])/scale
            v=coulomb['values'][j];muC=[rat(v[k])/rat(v[15]) for k in [1,2,3]]
            psie=rat(electron['values'][j,10])+rat(electron['values'][j,1])
            tc2=rat(data['station_tc2'][j]);offset=0;cell_bound=F(0);fisher=F(0)
            for el,Z in enumerate(g.d.CHARGES):
                positive=np.flatnonzero(pop[el,:Z+1]>0)
                if not active[el] or len(positive)<2:offset+=Z;continue
                reference=int(positive[np.argmax(pop[el,positive])]);ir=keys.index((el+1,reference))
                def primitive(name,z):return F(0) if z==0 else rat(data['station_'+name][j,offset+z-1])
                def dv(z):return F(0) if z==0 else primitive('dv',z)+rat(data['station_zero'][j,el])
                def pl(z):return rat(data['station_plop'][j,offset+z]) if z<Z else F(0)
                for z in positive:
                    z=int(z)
                    if z==reference:continue
                    iz=keys.index((el+1,z));dB=[row[iz]-row[ir] for row in B]
                    mu=sum(dB[k]*muC[k] for k in range(3))+dB[0]*psie
                    mu+=sum(dB[3+k]*gradient[k] for k in range(9))
                    mu-=L.get((el+1,z),F(0))-L.get((el+1,reference),F(0))
                    mu-=sum(dB[k]*gx[l] for l,k in enumerate([3,4,5,10]))
                    mu-=pl(z)-pl(reference)
                    # This is a free-energy gradient difference, not a replay of dv.
                    max_mu=max(max_mu,abs(mu+dv(z)-dv(reference)))
                    chemical=mu-primitive('ce',z)+primitive('ce',reference)+tc2*(primitive('binding',z)-primitive('binding',reference))
                    residual=iv.log(interval(rat(pop[el,z]))/interval(rat(pop[el,reference])))+interval(chemical)
                    bound=upper(residual);cell_bound=max(cell_bound,bound)
                    fisher+=rat(pop[el,z])*bound*bound;comparisons+=1
                offset+=Z
            max_equil=max(max_equil,cell_bound)
            bounds.append(dict(cell=i,comparisons=total_pairs,maximum_atomic_residual_upper=upper_text(cell_bound),
                Fisher_squared_residual_upper=upper_text(fisher)))
        count=record['stop'];print('ATOMIC STATIONARITY AUDIT',count,'/',plan['cells'],flush=True)
    assert count==5735 and comparisons==result['positive_atomic_comparisons']
    assert max_mu<F(str(plan['atomic_gradient_absolute_tolerance']))
    # Sum the two preregistered discrepancy budgets; do not re-fit a threshold.
    combined_budget=F(str(plan['atomic_gradient_absolute_tolerance']))+F(str(plan['atomic_log_population_tolerance']))
    assert max_equil<combined_budget
    e.save('frozen-residual-bounds.json',dict(classification='Proven',passed=True,rows=bounds,
        maximum_absolute_atomic_residual_upper=upper_text(max_equil),
        maximum_Fisher_squared_residual_upper=max((v['Fisher_squared_residual_upper'] for v in bounds),key=F),
        proof='For each positive atomic redistribution, evaluate the frozen free-energy gradient difference by rational component arithmetic and outward interval logarithms of the frozen positive populations. The summed nu_j times squared bounds upper-bounds the dual ideal-species Fisher norm after removal of the arbitrary per-element constant, since minimizing over that constant cannot increase the norm.',
        scope='Frozen primitive coefficients and populations, with molecular populations held fixed. Zero-comparison cells have no tested atomic tangent direction; their zero is not a certificate of complete equilibrium. Native primitive evaluation, continuous roots, omitted populations, molecular reactions and physical EOS error are not enclosed.'))
    e.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,positive_atomic_comparisons=comparisons,
        cells_without_tested_atomic_directions=cells_without_directions,original_21_EOS_outputs_and_atomic_populations_bitwise=True,
        previous_excitation_trace_bitwise=True,independent_gradient_constant_residual_upper=upper_text(max_mu),
        molecular_stationarity_certified=False,continuous_stationarity_certified=False,physical_EOS_certified=False))
    paths=[v for v in OUT.rglob('*') if v.is_file()]+[g.ROOT/'verification/eos_atomic_stationarity.py',g.ROOT/'verification/eos_atomic_stationarity_defined.py',g.ROOT/'verification/verify_atomic_stationarity.py']
    e.save('manifest.json',dict(classification='Counterexample candidate',sha256={v.relative_to(g.ROOT).as_posix():g.c.sha(v) for v in paths},
        runtime={str(e.LIB):g.c.sha(e.LIB),str(e.CACHE/'excitation.so'):g.c.sha(e.CACHE/'excitation.so')},
        original_library_sha256=plan['original_library_sha256']))
    print('PASS ATOMIC STATIONARITY AUDIT',count,'pairs',comparisons,'residual',upper_text(max_equil),flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items():assert g.c.sha(Path(path))==digest,path
    assert g.c.sha(p.c.s.LIB)==manifest['original_library_sha256']
    print('PASS ATOMIC STATIONARITY',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

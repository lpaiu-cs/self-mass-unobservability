"""Independent interval residuals for the full positive hydrogen inventory."""
from fractions import Fraction as F
from pathlib import Path
import sys
import numpy as np
from mpmath import iv
import eos_molecular_stationarity as m
import verify_atomic_stationarity as v
import verify_pressure_ionization as pi
from interval_records import exact_endpoint

e=m.e;p=m.p;g=m.g;OUT=e.OUT;rat=pi.rat;read=e.read


def run():
    assert not (OUT/'audit.json').exists();iv.dps=60
    plan=read(OUT/'plan.json');result=read(OUT/'result.json');bridge=read(OUT/'bridge.json')
    assert result['passed']
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(e.LIB)==bridge['library_sha256'] and g.c.sha(e.CACHE/'excitation.so')==bridge['bridge_sha256']
    assert g.c.sha(p.c.s.LIB)==plan['original_library_sha256']
    atomic=read(m.PREVIOUS_OUT/'frozen-residual-bounds.json')
    state=dict(np.load(g.OUT/'reference-state.npz'));eos=m.EOS();rows=[];count=0
    max_gradient=F(0);max_residual=F(0);max_fisher=F(0);directions=0
    zero_logs=[];zero_molecules=0;disabled_molecular_cells=0;without_H_directions=0
    units=np.r_[p.L**np.arange(7),1.,p.L**3];charges=[0,1,0,1];nuclei=[1,1,2,2]
    keys=[(1,0),(1,1),(25,0),(26,1)]
    for start in range(0,plan['cells'],128):
        record=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert record['start']==count and record['stop']==min(count+128,plan['cells'])
        assert g.c.sha(path)==record['output_sha256'] and g.c.sha(OUT/'plan.json')==record['plan_sha256']
        assert all(m.gates(r,plan) for r in record['rows'])
        data=dict(np.load(path));prior=[]
        for folder in [p.OUT,p.c.OUT,p.ex.OUT,p.c.s.OUT,m.PREVIOUS_OUT]:
            path=folder/f'block-{start}.npz';assert g.c.sha(path)==read(folder/f'block-{start}.json')['output_sha256']
            prior.append(dict(np.load(path)))
        mdh,coulomb,electron,species,old=prior
        assert all(np.array_equal(data[k],value) for k,value in old.items())
        for cell in plan['control_cells']:
            if start<=cell<record['stop']:
                control=dict(np.load(OUT/f'control-{cell}.npz'))
                assert all(np.array_equal(data[k][cell-start],value) for k,value in control.items())
        for j in range(record['stop']-start):
            i=start+j
            snap=eos.snapshot(state['lnd'][i],state['lnT'][i],state['X'][i])
            assert np.array_equal(snap['eos'],species['eos'][j]),i
            assert np.array_equal(snap['number_fractions'],species['number_fractions'][j]),i
            assert np.array_equal(snap['molecular_H_fractions'],species['molecular_H_fractions'][j]),i
            population=list(map(rat,[*species['number_fractions'][j,0,:2],*data['mol_populations'][j]]))
            active=[z for z,n in enumerate(population) if n>0];flags=data['mol_flags'][j]
            if not flags[3]:
                assert all(n==0 for n in population[2:]);disabled_molecular_cells+=1
            log_rho=iv.log(v.interval(rat(mdh['state'][j,1])))
            for z in range(2):
                if flags[3] and (z==0 or flags[1]>0) and population[2+z]==0:
                    # Saved pre-underflow log(n/N_A), not log of a returned zero.
                    log_nu=v.interval(rat(data['mol_logs'][j,z+1]))-log_rho
                    zero_logs.append(exact_endpoint(log_nu._mpi_[1]));zero_molecules+=1
            if len(active)<2:
                without_H_directions+=1;cell_bound=F(0);fisher=F(0)
            else:
                meta=mdh['state'][j];S=float(meta[1]*eos.constants[0]);k=rat(float(S*(4*np.pi/3)*p.L**3));quad=rat(meta[11])
                _,g2,_=pi.monomial(pi.P2,list(map(rat,mdh['extra'][j,0]/units)))
                _,g3,_=pi.monomial(pi.P3,list(map(rat,mdh['extra'][j,0]/units)))
                gradient=[k*x+quad*k*k*y for x,y in zip(g2,g3)]
                support=np.zeros((24,29));support[0,:2]=1.
                _,B,_,_=p.species_coordinates(support,np.ones(2),mdh['neutral'][j],mdh['ion3'][j])
                assert B.shape==(12,4);B=[[rat(x) for x in row] for row in B]
                factors=list(map(rat,S*np.array([1.,p.L,p.L*p.L,1.])))
                L={};gx=[F(0)]*4
                for component in range(int(data['count'][j])):
                    scale=rat(data['raw_scale'][j,component]);nu=rat(data['raw_value'][j,component,0])
                    tag=tuple(data['ids'][j,component]);L[tag]=rat(data['raw_value'][j,component,3])/scale
                    for z in range(4):gx[z]+=nu*factors[z]*rat(data['raw_grad'][j,component,z])/scale
                c=coulomb['values'][j];muC=[rat(c[z])/rat(c[15]) for z in [1,2,3]]
                psie=rat(electron['values'][j,10])+rat(electron['values'][j,1])
                PL=[rat(data['station_plop'][j,0]),F(0),*map(rat,data['mol_pl'][j])]
                mu=[]
                for z in range(4):
                    value=sum(B[l][z]*muC[l] for l in range(3))+charges[z]*psie
                    value+=sum(B[3+l][z]*gradient[l] for l in range(9))
                    value-=L.get(keys[z],F(0))+sum(B[l][z]*gx[q] for q,l in enumerate([3,4,5,10]))+PL[z]
                    mu.append(value)
                constants=[F(0),-rat(data['station_ce'][j,0])+rat(data['station_tc2'][j])*rat(data['station_binding'][j,0]),
                    -rat(data['mol_linear'][j,0]),-rat(data['mol_linear'][j,1])]
                density_coeff=[0,0,-1,-1]
                if flags[3]:
                    if population[0]>0 and population[2]>0:
                        max_gradient=max(max_gradient,abs(2*mu[0]-mu[2]-rat(data['mol_dv'][j,0])))
                    if flags[1]>0 and population[0]>0 and population[3]>0:
                        max_gradient=max(max_gradient,abs(2*mu[0]-mu[3]-sum(map(rat,data['mol_dv'][j]))))
                # One gauge for every positive H species, including molecular-only support.
                reference=max(active,key=lambda z:population[z]*nuclei[z]**2)
                cell_bound=F(0);fisher=F(0)
                for z in active:
                    if z==reference:continue
                    q=F(nuclei[z],nuclei[reference])
                    residual=iv.log(v.interval(population[z]))-v.interval(q)*iv.log(v.interval(population[reference]))
                    residual+=v.interval(constants[z]-q*constants[reference]+mu[z]-q*mu[reference])
                    residual+=v.interval(F(density_coeff[z])-q*density_coeff[reference])*log_rho
                    bound=v.upper(residual);cell_bound=max(cell_bound,bound);fisher+=population[z]*bound*bound;directions+=1
            # The prior total includes the old H contribution. Keeping it is a conservative overcount.
            combined_fisher=F(atomic['rows'][i]['Fisher_squared_residual_upper'])+fisher
            max_residual=max(max_residual,cell_bound);max_fisher=max(max_fisher,combined_fisher)
            rows.append(dict(cell=i,positive_H_species=len(active),H_directions=max(0,len(active)-1),
                H_residual_upper=v.upper_text(cell_bound),H_Fisher_squared_upper=v.upper_text(fisher),
                combined_positive_support_Fisher_squared_upper=v.upper_text(combined_fisher)))
        count=record['stop'];print('MOLECULAR STATIONARITY AUDIT',count,'/',plan['cells'],flush=True)
    assert count==5735
    assert max_gradient<F(str(plan['molecular_gradient_absolute_tolerance']))
    # A molecular reaction plus up to two atomic reference changes gives this predeclared-budget sum.
    budget=sum(F(str(plan[key]))*factor for key,factor in [
        ('molecular_gradient_absolute_tolerance',1),('molecular_log_population_tolerance',1),
        ('atomic_gradient_absolute_tolerance',2),('atomic_log_population_tolerance',2)])
    assert max_residual<budget
    omitted=None
    if zero_logs:
        log_bound=max(zero_logs);exponent=max(log_bound,F(-1000))
        omitted=dict(classification='Proven',zero_molecule_slots=zero_molecules,
            maximum_frozen_log_nu_upper=v.upper_text(log_bound),exponent_used=v.upper_text(exponent),
            frozen_preunderflow_nu_upper=v.upper_text(v.upper(iv.exp(v.interval(exponent)))),
            scope='Unreturned populations inferred from saved pre-underflow molecular logarithms on executed molecular branches only. The exponent is capped toward zero at -1000 to avoid enormous exact denominators. Native log evaluation, forced full-ionization branches and physical population error are not certified.')
    e.save('frozen-residual-bounds.json',dict(classification='Proven',passed=True,rows=rows,
        maximum_H_residual_upper=v.upper_text(max_residual),
        maximum_combined_positive_support_Fisher_squared_upper=v.upper_text(max_fisher),
        atomic_residual_bound=atomic['maximum_absolute_atomic_residual_upper'],omitted_molecules=omitted,
        proof='Assign nuclear counts A=(1,1,2,2) to H,H+,H2,H2+. Independently reconstruct each frozen chemical potential, including electron charge, MDH, Coulomb, scaled excitation and source-defined PL. Set c=(0,-ce_H+tc2*E_H,-linear_H2,-linear_H2+) and density coefficients d=(0,0,-1,-1). Outward intervals enclose g_i=log(nu_i)+c_i+d_i*log(rho_native)+mu_i. For any positive reference b, r_i=g_i-(A_i/A_b)*g_b spans every allowed H redistribution. Sum nu_i*r_i^2 bounds the minimized nuclear-constrained dual Fisher norm. Add the prior atomic bound conservatively, retaining its old H contribution as an overcount.',
        scope='The declared frozen primitive coefficients and current positive support, now including molecular redistribution. No uniform Hessian neighborhood, root enclosure, native primitive error, omitted-species derivative error or physical EOS certification is asserted.'))
    e.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,H_reaction_directions=directions,
        cells_without_H_directions=without_H_directions,cells_without_executed_molecular_branch=disabled_molecular_cells,
        original_21_EOS_outputs_atomic_and_molecular_populations_bitwise=True,previous_atomic_trace_bitwise=True,
        independent_molecular_gradient_constant_residual_upper=v.upper_text(max_gradient),
        continuous_stationarity_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    paths=[x for x in OUT.rglob('*') if x.is_file()]+[g.ROOT/'verification/eos_molecular_stationarity.py',g.ROOT/'verification/verify_molecular_stationarity.py']
    e.save('manifest.json',dict(classification='Counterexample candidate',sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in paths},
        runtime={str(e.LIB):g.c.sha(e.LIB),str(e.CACHE/'excitation.so'):g.c.sha(e.CACHE/'excitation.so')},
        original_library_sha256=plan['original_library_sha256']))
    print('PASS MOLECULAR STATIONARITY AUDIT',count,'H directions',directions,'residual',v.upper_text(max_residual),flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items():assert g.c.sha(Path(path))==digest,path
    assert g.c.sha(p.c.s.LIB)==manifest['original_library_sha256']
    print('PASS MOLECULAR STATIONARITY',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

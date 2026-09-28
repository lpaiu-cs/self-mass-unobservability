"""Independent constrained projection and rational frozen excitation bounds."""
from fractions import Fraction as F
import sys
import numpy as np
import excitation_scaled as scaled
import verify_pressure_ionization as pi

e=scaled.e
p=e.p; g=e.g; OUT=e.OUT; rat=pi.rat; read=e.read


def positive_scale(value):
    assert value>0
    exponent=(value.numerator.bit_length()-value.denominator.bit_length())//2
    return F(2)**exponent


def diagonal(n, B, atoms, elements):
    n=list(map(rat,n)); B=[[rat(x) for x in row] for row in B]
    result=[F(0)]*len(B)
    for element in set(elements):
        ix=[i for i,v in enumerate(elements) if v==element]
        den=sum(n[i]*int(atoms[i])**2 for i in ix)
        for k in range(len(B)):
            first=sum(n[i]*int(atoms[i])*B[k][i] for i in ix)
            result[k]+=sum(n[i]*B[k][i]**2 for i in ix)-first**2/den
    assert all(x>=0 for x in result)
    return result


def run():
    assert not (OUT/'audit.json').exists()
    plan=read(OUT/'plan.json'); bridge=read(OUT/'bridge.json')
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(e.LIB)==bridge['library_sha256']
    assert g.c.sha(e.CACHE/'excitation.so')==bridge['bridge_sha256']
    assert g.c.sha(p.c.s.LIB)==plan['original_library_sha256']
    for value in [F(1,10**1200),F(3,7),F(10**1200)]:
        scale=positive_scale(value); assert scale>0 and F(1,4)<value/(scale*scale)<4
    eos=p.EOS(); state=dict(np.load(g.OUT/'reference-state.npz'))
    units=np.r_[p.L**np.arange(7),1.,p.L**3]
    bounds=[]; spectra=[]; count=0
    max_H=max_G=max_spectrum=max_pair=0.
    for start in range(0,plan['cells'],128):
        record=read(OUT/f'block-{start}.json'); path=OUT/f'block-{start}.npz'
        assert record['start']==count and record['stop']==min(count+128,plan['cells'])
        assert g.c.sha(path)==record['output_sha256']
        assert g.c.sha(OUT/'plan.json')==record['plan_sha256']
        assert all(e.gates(r,plan) for r in record['rows'])
        data=dict(np.load(path)); previous=[]
        for folder in [p.OUT,p.c.OUT,p.ex.OUT,p.c.s.OUT]:
            path=folder/f'block-{start}.npz'
            assert g.c.sha(path)==read(folder/f'block-{start}.json')['output_sha256']
            previous.append(dict(np.load(path)))
        mdh,coulomb,electron,species=previous
        for cell in plan['control_cells']:
            if start<=cell<record['stop']:
                control=dict(np.load(OUT/f'control-{cell}.npz'))
                assert all(np.array_equal(data[k][cell-start],v) for k,v in control.items())
        for j in range(record['stop']-start):
            meta=mdh['state'][j]; v=coulomb['values'][j]
            S=float(meta[1]*eos.constants[0]); k=rat(float(S*(4*np.pi/3)*p.L**3)); q=rat(meta[11])
            u=list(map(rat,mdh['extra'][j,0]/units))
            _,_,h2=pi.monomial(pi.P2,u); _,_,h3=pi.monomial(pi.P3,u)
            H=[[F(0) for _ in range(16)] for _ in range(16)]
            for a,row in enumerate([[4,5,6],[5,7,8],[6,8,9]]):
                for b,index in enumerate(row): H[a][b]=rat(S)*rat(v[index])/rat(v[15])
            for a in range(9):
                for b in range(9): H[a+3][b+3]=k*h2[a][b]+q*k*k*h3[a][b]
            population=species['number_fractions'][j]
            x=state['X'][start+j]; ym=(x/g.c.A)@eos.mapping; eps=ym/float(ym@eos.weights)
            mol=eps[0]*species['molecular_H_fractions'][j]/2
            n,B,atoms,elements=p.species_coordinates(population,mol,mdh['neutral'][j],mdh['ion3'][j])
            keys=[(el+1,z) for el,Z in enumerate(g.d.CHARGES) for z in range(Z+1) if population[el,z]>0]
            keys += [(25+z,z) for z,vv in enumerate(mol) if vv>0]
            assert len(keys)==len(n)
            factors=S*np.array([1.,p.L,p.L*p.L,1.])
            D=np.zeros((4,len(n))); W=[[F(0) for _ in range(4)] for _ in range(4)]
            active_components=int(data['count'][j])
            for component in range(active_components):
                tag=tuple(data['ids'][j,component]); nu=rat(data['value'][j,component,0])
                if tag in keys:
                    D[:,keys.index(tag)]=factors*data['grad'][j,component]
                for a in range(4):
                    for b in range(4):
                        W[a][b]+=nu*rat(factors[a])*rat(data['raw_hess'][j,component,a,b])*rat(factors[b])/rat(data['raw_scale'][j,component])
            # Native coordinate order is alpha0, alpha1, alpha2, beta.
            selected=[3,4,5,10]
            for a,i in enumerate(selected):
                H[i][12+a]=H[12+a][i]=F(-1)
                for b,l in enumerate(selected): H[i][l]-=W[a][b]
            B=np.vstack([B,D])
            Hnon=np.array([[float(z) for z in row] for row in H])
            max_H=max(max_H,float(abs(Hnon-data['nonideal_Hessian'][j]).max()/max(1.,abs(Hnon).max())))
            allH=Hnon.copy(); numerator=rat(v[13])+rat(electron['values'][j,2])
            assert numerator>0 and v[10]>0 and v[12]>0
            allH[0,0]+=float(rat(S)*numerator/(rat(v[10])*rat(v[12])))
            max_H=max(max_H,float(abs(allH-data['Hessian'][j]).max()/max(1.,abs(allH).max())))
            # Compare explicit pairwise species derivatives with the low-rank pullback.
            C=B[selected]; direct=np.zeros((len(n),len(n)))
            for component in range(active_components):
                tag=tuple(data['ids'][j,component]); value=data['value'][j,component]
                gradient=factors*data['grad'][j,component]
                hh=factors[:,None]*data['hess'][j,component]*factors[None,:]
                direct-=value[0]*(C.T@hh@C)
                if tag in keys:
                    index=keys.index(tag); row=gradient@C
                    direct[index,:]-=row; direct[:,index]-=row
            ww=np.array([[float(vv) for vv in row] for row in W])
            lowrank=-(D.T@C+C.T@D+C.T@ww@C)
            max_pair=max(max_pair,float(abs(direct-lowrank).max()/max(1.,abs(direct).max())))
            A=np.array([np.where(elements==el,atoms,0) for el in sorted(set(elements))],dtype=np.longdouble)
            nn=n.astype(np.longdouble); bb=B.astype(np.longdouble)
            first=(bb*nn)@bb.T; cross=(bb*nn)@A.T; den=np.sum(A*A*nn,axis=1)
            G=first-(cross/den)@cross.T
            max_G=max(max_G,float(abs(G-data['Gram'][j]).max()/max(1.,abs(first).max())))
            spectrum=np.linalg.eigvals(np.asarray(G,float)@allH)
            assert abs(spectrum.imag).max()<1e-9
            margin=1+min(0.,float(spectrum.real.min())); spectra.append(margin)
            max_spectrum=max(max_spectrum,abs(margin-record['rows'][j]['combined_excitation_margin']))
            diag=diagonal(n,B,atoms,elements); active=[i for i,d in enumerate(diag) if d>0]
            scale={i:positive_scale(diag[i]) for i in active}
            assert all(vv>0 for vv in scale.values())
            loss=sum(abs(H[a][b])*(diag[a]*scale[b]/scale[a]+diag[b]*scale[a]/scale[b])/2 for a in active for b in active)
            bounds.append(pi.lower_text(1-loss))
        count=record['stop']; print('EXCITATION AUDIT',count,'/',plan['cells'],flush=True)
    assert count==5735 and max_H<1e-12 and max_G<1e-12 and max_pair<1e-12 and max_spectrum<1e-9
    result=read(OUT/'result.json'); assert abs(min(spectra)-result['minimum_combined_excitation_margin'])<1e-9
    minimum=min(map(F,bounds)); positive=minimum>0
    e.save('frozen-matrix-bounds.json',dict(classification='Proven',passed=True,
        minimum_lower_bound=pi.lower_text(minimum),positive_lower_bound=positive,lower_bounds=bounds,
        proof='Evaluate the constrained Gram diagonal exactly from frozen positive species and frozen coordinate coefficients. Sum absolute exact nonideal Hessian entries against weighted Young bounds on sqrt(G_ii G_jj). Add ideal-ion identity and the independently positive ideal-electron-plus-exchange correction. The MDH polynomial is differentiated by independent exact monomials; excitation Hessians use rational products of the saved native component derivatives.',
        scope='Finite reconstructed positive support, frozen MDH moments and coefficients, and frozen native excitation derivatives. This certifies the declared finite matrix, not native derivative accuracy, omitted populations, continuous curvature, equilibrium stationarity, physical EOS or GR evolution. A nonpositive lower bound is inconclusive, not proof of instability.'))
    e.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,
        independent_Hessian_scaled_error=max_H,independent_Gram_scaled_error=max_G,
        independent_pairwise_species_Hessian_scaled_error=max_pair,independent_spectrum_error=max_spectrum,
        exact_frozen_component_lower_bound_positive=positive,full_EOS_Hessian_certified=False,physical_EOS_certified=False))
    paths=[a for a in OUT.rglob('*') if a.is_file()]+[g.ROOT/'verification/excitation_curvature.py',g.ROOT/'verification/excitation_scaled.py',g.ROOT/'verification/verify_excitation_curvature.py']
    e.save('manifest.json',dict(classification='Counterexample candidate',sha256={a.relative_to(g.ROOT).as_posix():g.c.sha(a) for a in paths},
        runtime={str(e.LIB):g.c.sha(e.LIB),str(e.CACHE/'excitation.so'):g.c.sha(e.CACHE/'excitation.so')},original_library_sha256=plan['original_library_sha256']))
    print('PASS EXCITATION AUDIT',count,'frozen lower bound',pi.lower_text(minimum),'positive',positive,flush=True)


def verify():
    from pathlib import Path
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items(): assert g.c.sha(Path(path))==digest,path
    assert g.c.sha(p.c.s.LIB)==manifest['original_library_sha256']
    print('PASS EXCITATION',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

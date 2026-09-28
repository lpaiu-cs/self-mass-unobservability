"""Independent monomial arithmetic and exact frozen-component lower bounds."""
from fractions import Fraction as F
import json, math, sys
import numpy as np
import pressure_ionization_curvature as p
from verify_coulomb_curvature import lower_text

g=p.g;OUT=p.OUT
P2=[(6,(1,2)),(2,(0,3)),(1,(7,8))]
P3=[(3,(0,3,3)),(30,(1,2,3)),(9,(2,2,2)),(6,(0,2,4)),(9,(1,1,4)),(6,(0,1,5)),(1,(0,0,6))]


def read(path): return json.loads(path.read_text())
def rat(x): return F(float(x))


def monomial(terms,x):
    value=F(0);grad=[F(0)]*9;H=[[F(0) for _ in range(9)] for _ in range(9)]
    for coefficient,indices in terms:
        value+=coefficient*math.prod(x[k] for k in indices)
        for j,a in enumerate(indices):
            grad[a]+=coefficient*math.prod(x[k] for i,k in enumerate(indices) if i!=j)
            for l,b in enumerate(indices):
                if l!=j: H[a][b]+=coefficient*math.prod(x[k] for i,k in enumerate(indices) if i not in [j,l])
    return value,grad,H


def exact_diagonal(n,B,atoms,elements):
    n=list(map(rat,n));B=[[rat(v) for v in row] for row in B];diagonal=[F(0)]*12
    for element in set(elements):
        ix=[i for i,v in enumerate(elements) if v==element]
        den=sum(n[i]*int(atoms[i])**2 for i in ix)
        for k in range(12):
            first=sum(n[i]*int(atoms[i])*B[k][i] for i in ix)
            diagonal[k]+=sum(n[i]*B[k][i]**2 for i in ix)-first**2/den
    assert all(v>=0 for v in diagonal)
    return diagonal


def run():
    assert not (OUT/'audit.json').exists();plan=read(OUT/'plan.json');bridge=read(OUT/'bridge.json')
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(p.LIB)==bridge['library_sha256'] and g.c.sha(p.CACHE/'pi.so')==bridge['bridge_sha256']
    assert g.c.sha(p.c.s.LIB)==plan['original_library_sha256']
    eos=p.EOS();state=dict(np.load(g.OUT/'reference-state.npz'));count=0;bounds=[];spectra=[]
    max_H=0.;max_value=0.;max_G=0.;max_spectrum=0.
    units=np.r_[p.L**np.arange(7),1.,p.L**3]
    for start in range(0,plan['cells'],128):
        record=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert record['start']==count and record['stop']==min(count+128,plan['cells'])
        assert record['output_sha256']==g.c.sha(path) and record['plan_sha256']==g.c.sha(OUT/'plan.json')
        assert all(p.gates(r,plan) for r in record['rows'])
        a=dict(np.load(path));spfile=p.c.s.OUT/f'block-{start}.npz';species=dict(np.load(spfile))
        assert g.c.sha(spfile)==read(p.c.s.OUT/f'block-{start}.json')['output_sha256']
        cf=p.c.OUT/f'block-{start}.npz';coulomb=dict(np.load(cf));assert g.c.sha(cf)==read(p.c.OUT/f'block-{start}.json')['output_sha256']
        ef=p.ex.OUT/f'block-{start}.npz';electron=dict(np.load(ef));assert g.c.sha(ef)==read(p.ex.OUT/f'block-{start}.json')['output_sha256']
        for cell in plan['control_cells']:
            if start<=cell<record['stop']:
                control=dict(np.load(OUT/f'control-{cell}.npz'))
                assert all(np.array_equal(a[k][cell-start],v) for k,v in control.items())
        for j in range(record['stop']-start):
            meta=a['state'][j];raw=a['extra'][j];v=coulomb['values'][j]
            u=raw[0]/units;S=float(meta[1]*eos.constants[0]);kT=float(eos.constants[1]*meta[0])
            kappa=float(S*(4*np.pi/3)*p.L**3);k=rat(kappa);q=rat(meta[11])
            v2,g2,h2=monomial(P2,list(map(rat,u)));v3,g3,h3=monomial(P3,list(map(rat,u)))
            hpi=[[k*h2[i][l]+q*k*k*h3[i][l] for l in range(9)] for i in range(9)]
            f=rat(kT)*rat(S)*(k*v2+q*k*k*v3)
            max_value=max(max_value,float(abs(f-rat(a['computed'][j,0]))/max(F(1,10**300),abs(f))))
            H=[[F(0) for _ in range(12)] for _ in range(12)]
            for i,row in enumerate([[4,5,6],[5,7,8],[6,8,9]]):
                for l,index in enumerate(row): H[i][l]=rat(S)*rat(v[index])/rat(v[15])
            for i in range(9):
                for l in range(9): H[i+3][l+3]=hpi[i][l]
            allH=np.array([[float(z) for z in row] for row in H])
            numerator=rat(v[13])+rat(electron['values'][j,2]);assert numerator>0 and v[10]>0 and v[12]>0
            allH[0,0]+=float(rat(S)*numerator/(rat(v[10])*rat(v[12])))
            max_H=max(max_H,float(abs(allH-a['Hessian'][j]).max()/max(1.,abs(allH).max())))
            x=state['X'][start+j];ym=(x/g.c.A)@eos.mapping;eps=ym/float(ym@eos.weights)
            molecules=eps[0]*species['molecular_H_fractions'][j]/2
            n,B,atoms,elements=p.species_coordinates(species['number_fractions'][j],molecules,a['neutral'][j],a['ion3'][j])
            # Independent direct Schur projection; the driver uses centered sums.
            A=np.array([np.where(elements==el,atoms,0) for el in sorted(set(elements))],dtype=np.longdouble)
            nn=n.astype(np.longdouble);bb=B.astype(np.longdouble)
            first=(bb*nn)@bb.T;cross=(bb*nn)@A.T;den=np.sum(A*A*nn,axis=1)
            G=first-(cross/den)@cross.T
            max_G=max(max_G,float(abs(G-a['Gram'][j]).max()/max(1.,abs(first).max())))
            spectrum=np.linalg.eigvals(np.asarray(G,float)@allH);assert abs(spectrum.imag).max()<1e-9
            margin=1+min(0.,float(spectrum.real.min()));spectra.append(margin)
            max_spectrum=max(max_spectrum,abs(margin-record['rows'][j]['combined_component_margin']))
            diagonal=exact_diagonal(n,B,atoms,elements);active=[i for i,d in enumerate(diagonal) if d>0]
            scale={i:rat(math.sqrt(float(diagonal[i]))) for i in active};assert all(v>0 for v in scale.values())
            # Weighted Young inequality balances the disparate moment units.
            loss=sum(abs(H[i][l])*(diagonal[i]*scale[l]/scale[i]+diagonal[l]*scale[i]/scale[l])/2 for i in active for l in active)
            bounds.append(lower_text(1-loss))
        count=record['stop'];print('PRESSURE IONIZATION AUDIT',count,'/',plan['cells'],flush=True)
    assert count==5735 and max_H<1e-12 and max_value<1e-12 and max_G<1e-12 and max_spectrum<1e-9
    result=read(OUT/'result.json');assert abs(min(spectra)-result['minimum_combined_component_margin'])<1e-9
    minimum=min(map(F,bounds));positive=minimum>0
    p.save('frozen-matrix-bounds.json',dict(classification='Proven',passed=True,positive_lower_bound=positive,
        minimum_lower_bound=lower_text(minimum),lower_bounds=bounds,
        proof='The exact constrained Gram diagonal is evaluated from frozen positive species and frozen coordinate coefficients. For any positive scales s_i, sqrt(G_ii G_jj)<=0.5*(G_ii*s_j/s_i+G_jj*s_i/s_j). Summing this bound against absolute entries of the exact frozen Coulomb-plus-MDH Hessian bounds its normalized spectral norm. Add identity for ideal ions; the exact positive ideal-electron-plus-exchange correction cannot lower the result.',
        polynomial='MDH Hessians are evaluated by independent exact monomial products, including repeated indices, at the declared frozen moment coordinates and scaling coefficients.',
        scope='Exact finite coefficients, reconstructed species, radii and moments. Native continuous or physical accuracy, omitted species and excitation terms are not included. A nonpositive lower bound would mean this sufficient bound is inconclusive, not physical instability.'))
    p.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=count,
        independent_polynomial_value_relative_error=max_value,independent_Hessian_scaled_error=max_H,
        independent_Gram_scaled_error=max_G,independent_spectrum_error=max_spectrum,
        exact_frozen_component_lower_bound_positive=positive,excitation_included=False,physical_EOS_certified=False))
    paths=[a for a in OUT.rglob('*') if a.is_file()]+[g.ROOT/'verification/pressure_ionization_curvature.py',g.ROOT/'verification/verify_pressure_ionization.py']
    p.save('manifest.json',dict(classification='Counterexample candidate',sha256={a.relative_to(g.ROOT).as_posix():g.c.sha(a) for a in paths},
        runtime={str(p.LIB):g.c.sha(p.LIB),str(p.CACHE/'pi.so'):g.c.sha(p.CACHE/'pi.so')},original_library_sha256=plan['original_library_sha256']))
    print('PASS PRESSURE IONIZATION AUDIT',count,'frozen lower bound',lower_text(minimum),'positive',positive,flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in manifest['runtime'].items(): assert g.c.sha(p.c.s.Path(path))==digest,path
    assert g.c.sha(p.c.s.LIB)==manifest['original_library_sha256']
    print('PASS PRESSURE IONIZATION',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

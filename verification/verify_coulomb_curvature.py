"""Independent saved-state audit; exact bounds concern the frozen matrices only."""
from decimal import Decimal, localcontext, ROUND_FLOOR
from fractions import Fraction as F
import json, sys
import numpy as np
import coulomb_curvature as c

g=c.g; OUT=c.OUT


def read(path): return json.loads(path.read_text())
def rat(x): return F(float(x))


def lower_text(x):
    with localcontext() as ctx:
        ctx.prec=30;ctx.rounding=ROUND_FLOOR
        result=str(Decimal(x.numerator)/Decimal(x.denominator))
    assert F(result)<=x
    return result


def exact_diagonal(n, molecules):
    """Diagonal of the exact element-constrained Gram, in saved nu units."""
    diagonal=[F(0)]*3
    for element,row in enumerate(n):
        entries=[(rat(v),1,(j,int(j>0),j*j)) for j,v in enumerate(row) if v>0]
        if element==0:
            entries += [(rat(v),2,(j,j,j)) for j,v in enumerate(molecules) if v>0]
        if not entries: continue
        den=sum(v*a*a for v,a,b in entries)
        for k in range(3):
            first=sum(v*a*b[k] for v,a,b in entries)
            diagonal[k]+=sum(v*b[k]**2 for v,a,b in entries)-first**2/den
    assert all(x>=0 for x in diagonal)
    return diagonal


def run():
    assert not (OUT/'audit.json').exists()
    plan=read(OUT/'plan.json');result=read(OUT/'result.json')
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(c.s.LIB)==plan['library_sha256']
    assert g.c.sha(c.CACHE/'coulomb.so')==read(OUT/'bridge.json')['bridge_sha256']
    state=dict(np.load(g.ROOT/plan['state']));table=dict(np.load(g.OUT/'initial-adiabats-17.npz'))
    eos=c.EOS();constants=eos.constants.copy();records=[];parts=[];species=[];stop=0
    for start in range(0,plan['cells'],plan['block_cells']):
        row=read(OUT/f'block-{start}.json');path=OUT/f'block-{start}.npz'
        assert row['start']==stop and row['stop']==min(start+128,plan['cells'])
        assert row['plan_sha256']==g.c.sha(OUT/'plan.json') and row['state_sha256']==g.c.sha(g.ROOT/plan['state'])
        assert row['output_sha256']==g.c.sha(path)
        assert [r['cell'] for r in row['rows']]==list(range(row['start'],row['stop']))
        assert all(c.gates(r,plan) for r in row['rows'])
        sr=read(c.s.OUT/f'block-{start}.json');sp=c.s.OUT/f'block-{start}.npz'
        assert sr['output_sha256']==g.c.sha(sp)
        parts.append(dict(np.load(path)));species.append(dict(np.load(sp)));records+=row['rows'];stop=row['stop']
    assert stop==plan['cells']==5735
    data={k:np.concatenate([p[k] for p in parts]) for k in parts[0]}
    species={k:np.concatenate([p[k] for p in species]) for k in species[0]}
    assert np.array_equal(data['eos'],table['reference']) and np.array_equal(data['eos'],species['eos'])
    for i in plan['control_cells']:
        control=dict(np.load(OUT/f'control-{i}.npz'))
        for k in ['values','margins','Gram']: assert np.array_equal(data[k][i],control[k])
    max_gram=0.;max_spectrum=0.;bounds=[]
    for i,a in enumerate(data['values']):
        x=state['X'][i];ym=(x/g.c.A)@eos.mapping;cx=float(ym@eos.weights);eps=ym/cx
        # This reconstruction defines frozen populations. Reversing a rounded
        # rho/cx division is not a proof of the original native density.
        scale=float(data['eos'][i,0]*cx*constants[0])
        n=species['number_fractions'][i];mol=eps[0]*species['molecular_H_fractions'][i]/2
        G=np.zeros((3,3),dtype=np.longdouble)
        for j,row in enumerate(n):
            v=row.astype(np.longdouble);z=np.arange(29,dtype=np.longdouble);atoms=np.ones(29,dtype=np.longdouble)
            if j==0: v=np.r_[v,mol];z=np.r_[z,0.,1.];atoms=np.r_[atoms,2.,2.]
            if not v.any(): continue
            B=np.array([z,(z>0).astype(np.longdouble),z*z]);u=B@(v*atoms)
            G+=(B*v)@B.T-np.outer(u,u)/np.sum(v*atoms*atoms)
        G*=np.longdouble(scale);old=data['Gram'][i]
        gram_error=float(abs(G-old).max()/max(1.,abs(old).max(),float(np.sum(n)*scale)))
        max_gram=max(max_gram,gram_error)
        H=np.array([[a[4],a[5],a[6]],[a[5],a[7],a[8]],[a[6],a[8],a[9]]])/a[15]
        for k in range(2):
            v=np.linalg.eigvals(np.asarray(G,float)@H)
            assert abs(v.imag).max()<1e-10
            margin=1+min(0.,float(v.real.min()))
            max_spectrum=max(max_spectrum,abs(margin-data['margins'][i,k]))
            H[0,0]+=a[13]/(a[10]*a[12])
        diag=exact_diagonal(n,mol)
        h=[[rat(a[k]) for k in row] for row in [[4,5,6],[5,7,8],[6,8,9]]]
        loss=sum(abs(h[j][k])*(diag[j]+diag[k])/2 for j in range(3) for k in range(3))*rat(scale)/rat(a[15])
        bounds.append(lower_text(1-loss))
        if i%512==0: print('COULOMB INDEPENDENT AUDIT',i,'/',len(data['values']),flush=True)
    assert max_gram<1e-12 and max_spectrum<1e-10,(max_gram,max_spectrum)
    assert np.array_equal(data['margins'].min(0),result['minimum_component_margins'])
    minimum=min(map(F,bounds));assert minimum>0
    c.save('exact-matrix-bounds.json',dict(classification='Proven',lower_bounds=bounds,minimum_lower_bound=lower_text(minimum),
        constants=constants.tolist(),
        proof='For the exact reconstructed populations, G is a positive Gram matrix. The spectral norm of sqrt(G) H sqrt(G) is at most sum_ij |H_ij| sqrt(G_ii G_jj), at most sum_ij |H_ij|(G_ii+G_jj)/2. Exact rational diagonals, Hessian entries and scale give a conservative lower bound on I plus the Coulomb correction. The saved ideal-electron compressibility is positive and cannot reduce this bound.',
        scope='Exact frozen binary64 coefficients and reconstructed populations, including only their saved positive support. No enclosure of native evaluation, omitted populations, continuous EOS domain, other nonideal terms or physical-model error.'))
    assert np.all(data['values'][:,[10,12,13,15]]>0)
    c.save('audit.json',dict(classification='Counterexample candidate',passed=True,cells=stop,
        all_21_outputs_bitwise_equal_prior_table=True,serial_controls_bitwise=True,
        independent_unscaled_Gram_error=max_gram,independent_spectrum_error=max_spectrum,
        minimum_frozen_matrix_lower_bound=lower_text(minimum),
        full_EOS_Hessian_certified=False,physical_EOS_certified=False))
    paths=[p for p in OUT.rglob('*') if p.is_file()]+[g.ROOT/'verification/coulomb_curvature.py',g.ROOT/'verification/verify_coulomb_curvature.py']
    c.save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},library_sha256=plan['library_sha256']))
    print('PASS COULOMB INDEPENDENT AUDIT',stop,'minimum exact frozen-matrix bound',lower_text(minimum),flush=True)


def verify():
    manifest=read(OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(c.s.LIB)==manifest['library_sha256']
    print('PASS COULOMB',len(manifest['sha256']),'SHA bindings',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()

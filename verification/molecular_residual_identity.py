"""Exact nuclear-constrained reaction bases for positive H/H+/H2/H2+ support."""
import json,sys
import sympy as sp
import direct_eos_gr as g

OUT=g.OUT/'molecular-residual-identities'


def run():
    assert not OUT.exists();OUT.mkdir()
    A=[1,1,2,2];checked=0
    for mask in range(1,16):
        active=[i for i in range(4) if mask&(1<<i)];aa=sp.Matrix([[A[i] for i in active]])
        for reference in range(len(active)):
            columns=[]
            for i in range(len(active)):
                if i==reference:continue
                column=sp.zeros(len(active),1);column[i]=1
                column[reference]=-sp.Rational(A[active[i]],A[active[reference]])
                assert aa*column==sp.zeros(1,1);columns.append(column)
            basis=sp.Matrix.hstack(*columns) if columns else sp.zeros(1,0)
            assert basis.rank()==len(active)-1;checked+=1
    n=sp.symbols('n0:4',positive=True);mu=sp.symbols('g0:4',real=True);lam=sp.symbols('lam',real=True)
    den=sum(n[i]*A[i]**2 for i in range(4));num=sum(n[i]*A[i]*mu[i] for i in range(4))
    norm=sum(n[i]*(mu[i]-A[i]*lam)**2 for i in range(4))
    assert sp.simplify(sp.diff(norm,lam).subs(lam,num/den))==0
    assert sp.expand(sp.diff(norm,lam,2)-2*den)==0
    assert sp.simplify(norm-norm.subs(lam,num/den)-den*(lam-num/den)**2)==0
    dv0,dv1=2*mu[0]-mu[2],mu[2]-mu[3]
    assert sp.expand(dv0+dv1-(2*mu[0]-mu[3]))==0
    result=dict(classification='Proven',passed=True,positive_support_reference_bases_checked=checked,
        nuclear_counts=A,
        basis='For every nonempty positive support and any reference b, e_i-(A_i/A_b)*e_b is nuclear conserving. These vectors have rank |support|-1 and span every internal H redistribution.',
        projection='The minimum of sum nu_i*(g_i-A_i*lambda)^2 occurs at lambda=sum(nu_i*A_i*g_i)/sum(nu_i*A_i^2). Any chosen positive reference gives an upper bound on this minimum.',
        molecular_constants='The H2+ native correction is relative to H2; adding the H2 formation correction gives the H2+ correction relative to two neutral H atoms, with the electron chemical potential already included by charge pullback.',
        scope='Exact linear algebra and projection identities. They do not certify primitive EOS evaluation, current-state numerical derivatives, omitted species, a continuous root neighborhood, or physical stellar stability.')
    (OUT/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    paths=[OUT/'result.json',g.ROOT/'verification/molecular_residual_identity.py']
    (OUT/'manifest.json').write_text(json.dumps(dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths}),indent=2)+'\n')
    print('PASS MOLECULAR REACTION BASES AND FISHER PROJECTION',checked,flush=True)


def verify():
    manifest=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in manifest['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS MOLECULAR RESIDUAL IDENTITIES',len(manifest['sha256']),'artifact SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

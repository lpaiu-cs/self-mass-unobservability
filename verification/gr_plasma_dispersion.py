"""Positive moment series for the full leading-order transverse/longitudinal kernels."""
import json, sys
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad_vec
from scipy.special import expit
import gr_fermi_plasma as plasma

g=plasma.g;OUT=g.OUT/'gr-plasma-dispersion'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();plasma.verify()
    r=json.loads((plasma.OUT/'result.json').read_text())
    assert r['quadrature_passed'] and r['native_density_passed']
    paths=[g.ROOT/'verification/gr_plasma_dispersion.py',g.ROOT/'verification/gr_fermi_plasma.py',plasma.OUT/'manifest.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='98e19da',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source='https://arxiv.org/abs/hep-ph/9302213',equations=[51,52,54,55],
        model='Full leading-order alpha polarization kernels of the declared ideal electron/positron distribution, with finite electron mass. Do not replace the distribution by a characteristic velocity.',
        orders=[16,24],transverse_q=[0,.01,.1,1,10,1000],longitudinal_fractions=[0,.01,.5,.99,.999999],
        actual_controls=[0,1175,2972,3043,4352,5734],relative_absolute_tolerance=1e-10,
        quadrature='All 5735 states: two tolerances for positive velocity moments through order 25 and a longitudinal tail-bound moment. Compare order 16/24 roots. At six actual states, independently integrate original logarithmic kernels in momentum coordinate at 50 decimals at every registered nonzero wave number.',
        acceptance='Finite normalized moment/refinement/residual checks use 1e-10. Series truncation bounds are exact-real theorems; their evaluated values are diagnostics until interval bounds on moments and roundoff are supplied.',
        scope='Conditional full leading-order dispersion roots. Not collision damping, physical all-order plasma EOS, a complete longitudinal thermodynamic state count, non-LTE transport or actual GR evolution.'))


def symbolic():
    z,w=sp.symbols('z w',positive=True)
    A=sp.atanh(sp.sqrt(w*z))/sp.sqrt(w*z)
    Bt=1/z-(1-z)*A/z
    Bl=(2*A-1-(1-z)/(1-w*z))/z
    ct=lambda j:w**j/(2*sp.Integer(j)+1)-w**(j+1)/(2*sp.Integer(j)+3)
    assert sp.series(Bt-sum(ct(j)*z**j for j in range(5)),z,0,5).removeO()==0
    assert sp.series(Bl-sum((2*j+1)*ct(j)*z**j for j in range(5)),z,0,5).removeO()==0
    for N in [1,4,16]:
        assert sp.simplify(sum(ct(j) for j in range(N+1))-(1-w**(N+1)/(2*N+3)))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        measure='dM=(p^2/E)*(f_e+f_p)dp; w=(p/E)^2; I_p=integral(1-w/3)dM>0. Define c_j=integral[w^j/(2j+1)-w^(j+1)/(2j+3)]dM/I_p, hence c_0=1 and c_j>=0.',
        exact_kernels='T(z)=sum c_j*z^j=Pi_t/w_p^2; L(z)=sum(2j+1)*c_j*z^j=Pi_l/(z*w_p^2), z=k^2/w^2 in [0,1]. T(1)=I_m/I_p and L(1)=I_k/I_p.',
        transverse_tail='0<=T-T_N<=z^(N+1)*integral w^(N+1)dM/[(2N+3)*I_p]. At z=1 this is an equality by telescoping. No geometric 1/(1-w) loss is needed.',
        longitudinal_tail='0<=L-L_N<=z^(N+1)*integral w^(N+1)*[1+2*w/((2N+5)*(1-w))]dM/I_p. This follows by telescoping to w^(N+1)+2*sum_{j>=N+2}w^j/(2j+1) and bounding its remaining denominators below by 2N+5.',
        roots='Write h=u^2-q^2 for transverse modes: h=T(q^2/(q^2+h)), 1<=h<=T(1). For longitudinal modes write h=u^2: h=L(q^2/h), max(1,q^2)<=h<=L(1), 0<=q<sqrt(L(1)). Each residual is strictly increasing with derivative >=1; endpoints bracket exactly one root. At q=0 both have u=1. At q=sqrt(L(1)) the longitudinal timelike branch reaches u=q.',
        residual_error='A bound on the absolute residual, including kernel truncation, integral evaluation and floating-point error, bounds the error in h by the same amount. Then |delta u|<=|delta h|/(2*min(u_true,u_hat)). The stored finite diagnostics omit an interval roundoff certificate.',
        scope='The theorems hold for these real leading-order kernels and nonnegative measures with finite stated moments. The longitudinal boundary is that of the timelike branch used in the source; it does not certify all damped/spacelike plasma excitations.'))


def moments(eta,beta,scale,N,tol):
    eta,beta,scale=map(np.asarray,[eta,beta,scale])
    def fn(t):
        gamma=1+beta*t;x=np.sqrt(beta*t*(2+beta*t));w=x*x/gamma**2
        weight=beta*x*(expit(eta-t)+expit(-eta-t-2/beta))/scale
        # Store raw positive velocity moments, plus one positive longitudinal tail bound.
        values=[weight];power=np.ones_like(w)
        for j in range(1,N+2):power=power*w;values.append(weight*power)
        values.append(weight*power*(1+2*x*x/(2*N+5)))
        return np.array(values)
    return quad_vec(fn,0,np.inf,epsabs=tol,epsrel=tol,limit=600)


def coefficients(raw,N):
    ip=raw[0]-raw[1]/3
    return np.array([(raw[j]/(2*j+1)-raw[j+1]/(2*j+3))/ip for j in range(N+1)])


def polynomial(coeff,z):
    ans=np.zeros_like(z)
    for x in coeff[::-1]:ans=ans*z+x
    return ans


def roots(coeff,q,kind):
    q=np.asarray(q);q2=q*q
    hi=np.sum(coeff,axis=0)+np.zeros_like(q)
    lo=np.ones_like(hi) if kind=='T' else np.maximum(1,q2)+np.zeros_like(hi)
    # Polynomial longitudinal cutoff may differ from the exact one by its omitted tail.
    assert np.all(hi>=lo)
    for _ in range(60):
        h=(hi+lo)/2;z=q2/(q2+h) if kind=='T' else q2/h
        residual=h-polynomial(coeff,z);positive=residual>=0
        hi=np.where(positive,h,hi);lo=np.where(positive,lo,h)
    h=(hi+lo)/2;z=q2/(q2+h) if kind=='T' else q2/h
    return h,z,float(np.max(abs(h-polynomial(coeff,z))))


def independent(eta,beta,scale,q,h,kind):
    mp.mp.dps=50;eta,beta,scale,q,h=map(lambda x:mp.mpf(str(x)),[eta,beta,scale,q,h])
    z=q*q/(q*q+h) if kind=='T' else q*q/h
    def fn(x,base):
        gamma=mp.sqrt(1+x*x);v=x/gamma;w=v*v;t=x*x/((gamma+1)*beta)
        weight=x*x/gamma*(1/(1+mp.exp(t-eta))+1/(1+mp.exp(t+eta+2/beta)))/scale
        if base or not x or not z:return weight*(1-w/3)
        r=v*mp.sqrt(z);A=mp.atanh(r)/r
        kernel=(1-(1-z)*A)/z if kind=='T' else (2*A-1-(1-z)/(1-w*z))/z
        return weight*kernel
    points=[mp.mpf(0)]+[mp.sqrt(beta*t*(2+beta*t)) for t in [1,max(2,float(eta)),max(2,float(eta))+8,max(2,float(eta))+40]]+[mp.inf]
    den=mp.quad(lambda x:fn(x,True),points);val=mp.quad(lambda x:fn(x,False),points)/den
    return str(val),float(h-val)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    symbolic();data=dict(np.load(plasma.OUT/'stellar-fermi-plasma.npz'))
    eta=data['eta'];beta=data['beta'];scale=data['dimensionless_density_scale'];N=max(plan['orders'])
    low,err1=moments(eta,beta,scale,N,1e-9);raw,err2=moments(eta,beta,scale,N,2e-12)
    score=float(np.max(abs(raw-low)/np.maximum(abs(raw),1e-10)))
    ip=raw[0]-raw[1]/3;coeff=coefficients(raw,N);Lcoeff=coeff*(2*np.arange(N+1)+1)[:,None]
    exactMt=data['normalized_moments'][2]/data['normalized_moments'][1]
    exactK=data['normalized_moments'][4]/data['normalized_moments'][1]
    tailT=raw[N+1]/((2*N+3)*ip);tailL=raw[N+2]/ip
    all_h={};summaries=[];controls=[];rows=[]
    for kind,qs in [('T',plan['transverse_q']),('L',plan['longitudinal_fractions'])]:
        for index,qvalue in enumerate(qs):
            q=np.full(len(eta),qvalue) if kind=='T' else qvalue*np.sqrt(exactK)
            co=coeff if kind=='T' else Lcoeff;h,z,residual=roots(co,q,kind)
            h16,_,res16=roots(co[:plan['orders'][0]+1],q,kind)
            difference=float(np.max(abs(h-h16)/np.maximum(abs(h),1)))
            all_h[f'{kind}-{index}-h']=h;all_h[f'{kind}-{index}-q']=q
            # Positive-kernel derivative gives the declared group velocity; finite check only.
            derivative=polynomial(co[1:]*np.arange(1,N+1)[:,None],z)
            u=np.sqrt(q*q+h) if kind=='T' else np.sqrt(h)
            group=(q/u)*(1+derivative/u**2)/(1+derivative*z/u**2) if kind=='T' else (q/u)*derivative/(h+derivative*z)
            tail=tailT if kind=='T' else tailL
            rows.append(dict(kind=kind,q_or_fraction=qvalue,order_difference=difference,polynomial_residual=residual,
                evaluated_tail_upper=float(np.max(tail*z**(N+1))),group_velocity_range=[float(group.min()),float(group.max())]))
            for i in plan['actual_controls']:
                if not qvalue:continue
                ref,res=independent(eta[i],beta[i],scale[i],q[i],h[i],kind)
                controls.append(dict(cell=i,kind=kind,q=float(q[i]),h=float(h[i]),logarithmic_integral_50decimal=ref,
                    residual=res,passed=abs(res)<plan['relative_absolute_tolerance']))
            if kind=='T':
                v=data['v_star'];approx=np.ones_like(q)
                for _ in range(30):
                    r=v*q/np.sqrt(q*q+approx);w=r*r
                    approx=np.where(r<.02,1+w/5+3*w*w/35+w**3/21+w**4/33,
                        3/(2*np.maximum(w,1e-100))*(1-(1-w)*plasma.atanh_over_v(np.maximum(r,1e-100))))
                summaries.append(dict(q=qvalue,max_characteristic_velocity_h_difference=float(np.max(abs(approx-h)))))
    np.savez_compressed(OUT/'stellar-dispersion.npz',coefficients=coeff,raw_moments=raw,
        transverse_tail_bound=tailT,longitudinal_tail_bound=tailL,**all_h)
    gates=bool(score<plan['relative_absolute_tolerance'] and np.all(coeff>=0) and all(x['passed'] for x in controls)
        and all(x['order_difference']<plan['relative_absolute_tolerance'] and x['polynomial_residual']<plan['relative_absolute_tolerance'] for x in rows))
    save('result.json',dict(classification='Counterexample candidate',completed=True,finite_gates_passed=gates,cells=len(eta),
        quadrature_difference=score,quadrature_error_estimates=[float(err1),float(err2)],rows=rows,controls=controls,
        maximum_independent_residual=max(abs(x['residual']) for x in controls),
        maximum_transverse_series_tail_diagnostic=float(tailT.max()),maximum_longitudinal_series_tail_diagnostic=float(tailL.max()),
        transverse_endpoint_identity_residual=float(np.max(abs(exactMt-np.sum(coeff,axis=0)-tailT))),
        longitudinal_endpoint_difference=float(np.max(abs(exactK-np.sum(Lcoeff,axis=0)))),
        characteristic_velocity_comparison=summaries,moment_and_roundoff_intervals_certified=False,
        physical_all_order_plasma_EOS_certified=False,native_EOS_replaced=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS full leading-order dispersion provenance; finite checks are not physical certificates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

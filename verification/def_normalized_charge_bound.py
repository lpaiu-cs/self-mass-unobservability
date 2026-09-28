"""Outward-interval derivative of the declared, stored leading charge model.

The exact inputs are the decimal strings of saved numbers, including native
EOS derivatives. Their physical and coefficient-generation errors are outside
this finite-dimensional certificate. No new trajectory or EOS call is made.
"""
from pathlib import Path
import json
import time
import mpmath as mp
import numpy as np
import def_normalized_charge as q

OUT=q.OUT/'derivative-bound'


def solve(d,c,rhs):
    a=list(d);b=list(rhs)
    for i in range(1,len(a)):
        assert a[i-1].a>0
        f=c[i-1]/a[i-1];a[i]-=f*c[i-1];b[i]-=f*b[i-1]
    assert a[-1].a>0
    b[-1]/=a[-1]
    for i in range(len(a)-2,-1,-1):b[i]=(b[i]-c[i]*b[i+1])/a[i]
    return b


def compute():
    mp.iv.dps=70;iv=lambda x:mp.iv.mpf(str(x))
    zero=iv(0);one=iv(1);two=iv(2)
    # Analytic positive-pivot control for the same interval solver.
    control=solve([iv(2),iv(2)],[iv(-1)],[one,zero])
    assert control[0].a<=iv(2)/3<=control[0].b and control[1].a<=one/3<=control[1].b
    data=np.load(q.OUT/'zero.npz');initial=np.load(q.s.OUT/'initial.npz')
    A=lambda xs:[iv(x) for x in xs]
    E,P,rho,a,b,bf=A(data['E']),A(data['P']),A(data['rho']),A(data['a']),A(data['b']),A(data['bf'])
    r,rf,V=A(initial['r']),A(initial['rf']),A(initial['volume']);n=len(r);R=rf[-1]
    # Use the same stored geometric fraction and area conventions as the kernel.
    star,_=q.s.initialize(q.s.old.imported.CachedOnly(),-4,q.ld('.001'),q.ld(1))
    fraction=A(star.fraction);gv=[iv(q.e.GRAV)*v for v in V]
    distance=A(star.distance);area=A(star.scalar_area)
    beta=iv(-4);pi=iv(np.pi);mu=iv(data['mf'][-1])/R
    er=[E[i]+rho[i]*iv(data['raw'][i,9]) for i in range(n)]
    et=[rho[i]*iv(data['raw'][i,10]) for i in range(n)]
    pr=[P[i]*iv(data['raw'][i,5]) for i in range(n)]
    pt=[P[i]*iv(data['raw'][i,6]) for i in range(n)]
    W=[gv[i]*(E[i]+P[i])/(b[i]*r[i]) for i in range(n)]
    C=[1-W[i]*(1-fraction[i]) for i in range(n)]
    D=[1+W[i]*fraction[i] for i in range(n)]
    H=[zero]*n;lam=[zero]*(n+1);lam[-1]=one
    for i in range(n-1,-1,-1):
        assert C[i].a>0 and D[i].a>0
        H[i]=lam[i+1]/D[i];lam[i]=H[i]*C[i]
    c=[area[i+1]/(4*pi*R)*(H[i]+H[i+1])/2*bf[i+1]/(r[i+1]-r[i]) for i in range(n-1)]
    L=-mp.iv.ln(1-2*mu)/(2*mu);internal=distance[-1]/(H[-1]*r[-1]*bf[-1])
    ce=1/(internal+L)
    potential=[gv[i]/R*H[i]*beta*(-E[i]+3*P[i]) for i in range(n)]
    diagonal=[([zero]+c)[i]+(c+[ce])[i]-potential[i] for i in range(n)]
    w=solve(diagonal,[-v for v in c],potential)
    flux=-sum(potential[i]*(1+w[i]) for i in range(n));response=flux/mu
    gradient=[];mass_term=[];source_term=[]
    for j in range(n):
        dm=[zero]*(n+1);de=[zero]*n;dr=[zero]*n
        for i in range(n):
            k=a[i]**2/r[i];f=fraction[i];t=one if i==j else zero
            dm[i+1]=((1-gv[i]*er[i]*k*(1-f))*dm[i]+gv[i]*et[i]*t)/(1+gv[i]*er[i]*k*f)
            de[i]=(dm[i+1]-dm[i])/gv[i]
            dr[i]=-k*((1-f)*dm[i]+f*dm[i+1])
        dp=[pr[i]*dr[i]+(pt[i] if i==j else zero) for i in range(n)]
        db=[-2*((1-fraction[i])*dm[i]+fraction[i]*dm[i+1])/r[i] for i in range(n)]
        dbf=[zero]+[-2*dm[i+1]/rf[i+1] for i in range(n)]
        dW=[W[i]*((de[i]+dp[i])/(E[i]+P[i])-db[i]/b[i]) for i in range(n)]
        dC=[-(1-fraction[i])*dW[i] for i in range(n)];dD=[fraction[i]*dW[i] for i in range(n)]
        dH=[zero]*n;dlam=[zero]*(n+1)
        for i in range(n-1,-1,-1):
            dH[i]=dlam[i+1]/D[i]-H[i]*dD[i]/D[i]
            dlam[i]=dH[i]*C[i]+H[i]*dC[i]
        dc=[c[i]*((dH[i]+dH[i+1])/(H[i]+H[i+1])+dbf[i+1]/bf[i+1]) for i in range(n-1)]
        dmu=dm[-1]/R
        dL=(1/(1-2*mu)-L)/mu*dmu
        dce=-ce**2*(-internal*(dH[-1]/H[-1]+dbf[-1]/bf[-1])+dL)
        dpv=[gv[i]/R*beta*(dH[i]*(-E[i]+3*P[i])+H[i]*(-de[i]+3*dp[i])) for i in range(n)]
        dd=[([zero]+dc)[i]+(dc+[dce])[i]-dpv[i] for i in range(n)]
        rhs=[dpv[i]-dd[i]*w[i] for i in range(n)]
        for i in range(n-1):rhs[i]+=dc[i]*w[i+1];rhs[i+1]+=dc[i]*w[i]
        dw=solve(diagonal,[-v for v in c],rhs)
        df=-sum(dpv[i]*(1+w[i])+potential[i]*dw[i] for i in range(n))
        source_term.append(df/mu);mass_term.append(-flux*dmu/mu**2)
        gradient.append(source_term[-1]+mass_term[-1])
    interval=lambda x:[str(x.a),str(x.b)]
    norm=sum(abs(v) for v in gradient)
    numerical=np.load(q.OUT/'gradient.npz')['gradient']
    discrepancies=[max(abs(float(g.a)-v),abs(float(g.b)-v)) for g,v in zip(gradient,numerical)]
    assert max(discrepancies)<1e-10*float(norm.a)
    # The norm ceiling is a rounded rational bound verified by interval arithmetic.
    ceiling=iv('0.00001144');assert norm.b<ceiling.a
    benchmark=json.loads((q.s.OUT/'companion-benchmark.json').read_text())
    du=iv(benchmark['leading_drive']['maximum_delta_u'])
    target=iv('1e-9');threshold=target/(iv('.001')**2*du*ceiling)
    return dict(classification='Proven',passed=True,exact_model='Stored decimal coefficient model of leading(star, zero, theta), differentiated at theta=0 with fixed B, composition and cell radii. The scalar matrix, lapse and mass normalization all vary.',
        response=interval(response),gradient=[interval(v) for v in gradient],
        gradient_l1=interval(norm),rational_gradient_l1_upper_bound='0.00001144',
        maximum_complex_step_absolute_difference=max(discrepancies),
        source_contribution_sum=interval(sum(source_term)),mass_normalization_contribution_sum=interval(sum(mass_term)),
        necessary_temperature_linf=interval(target/ceiling),
        necessary_normalized_fluid_gain=interval(threshold),
        positive_pivots_and_analytic_control=True,
        scope='Conditional linear derivative bound only. No physical/native derivative error, finite-temperature remainder, finite-phi0 remainder, moving radius/composition, orbital fluid resolvent or observable/no-go certificate.')


if __name__=='__main__':
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    files=[Path(__file__),Path(q.__file__),q.OUT/'zero.npz',q.OUT/'gradient.npz',q.OUT/'manifest.json',q.s.OUT/'initial.npz',q.s.OUT/'companion-benchmark.json']
    q.e.write(OUT/'plan.json',dict(classification='Proven',arithmetic='70-decimal outward mpmath intervals with explicit positive-pivot elimination and analytic implicit differentiation; stored number strings are exact declared inputs.',bindings={p.relative_to(q.s.ROOT).as_posix():q.e.digest(p) for p in files},budget_seconds=60))
    result=compute();result['seconds']=time.monotonic()-start;q.e.write(OUT/'result.json',result)
    q.e.write(OUT/'manifest.json',dict(sha256={p.relative_to(q.s.ROOT).as_posix():q.e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))
    print(json.dumps({k:v for k,v in result.items() if k!='gradient'}))

"""A connected-branch primitive uniqueness theorem and exact frozen-jet audit."""
from fractions import Fraction as F
import gzip,json,sys
import numpy as np
import sympy as sp
import gr_heat_primitive_newton as previous

old=previous.original;heat=old.heat;g=old.g;OUT=g.OUT/'gr-heat-primitive-global';sha=g.c.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    previous.verify();assert not OUT.exists();OUT.mkdir()
    files=[g.ROOT/'verification/gr_heat_primitive_global.py',previous.OUT/'manifest.json',previous.OUT/'result.json',
        old.OUT/'symbolic.json',heat.old.char.OUT/'coefficients.npz',heat.OUT/'initial-rates.npz',heat.OUT/'manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='1f4892ed',cells=5735,
        bindings={p.relative_to(g.ROOT).as_posix():sha(p) for p in files},
        target='Prove uniqueness of the fixed-(D,E,J,Q,X) point primitive on a connected admissible thermal branch. Determine the exact minimum of the frozen normalized inverse determinant on the entire closed velocity interval[-1,1], extending the earlier |v|<=1/2 diagnostic without asserting frozen coefficients equal a physical EOS neighbourhood.',
        assumptions='Natural units c=1. Fixed composition X and specified comoving heat Q. Positive D, rho,T,w and b=T*epsilon_T/w. At each density the specified temperature domain is an interval on which a C1 EOS has epsilon_T>0 throughout; the energy equation then defines at most one T(v) in that domain. The admissible velocities with T(v)>0 form the explicitly required connected interval I subset(-1,1). The pointwise determinant condition must hold at every possible root on that branch.',
        quantities='b=T*epsilon_T/w; r=rho*P_rho/w; d=T*P_T/w; e=1-rho*epsilon_rho/w; j=Q/w; g=r+d*e/b. No substitution e=d is made for the stored native derivatives.',
        gates='Record positive,zero,negative margins; never discard failed cells or replace native derivatives by a Maxwell identity. The old native27-case recovery result and all failures preceding its Newton repair remain frozen.',
        boundary='A theorem with explicit branch-wide EOS premises plus an exact audit of frozen jets. Neither certifies actual native EOS neighbourhoods, root existence without a bracket, unresolved cell averages, the heat closure, atmosphere or a finite GR trajectory.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(g.ROOT/rel)==digest,rel
    return p


def symbolic():
    v,b,r,d,e,j,w,E,J,Q=sp.symbols('v b r d e j w E J Q',real=True);W2=1/(1-v*v);sound=r+d*e/b
    # At fixed conserved D,E,J,Q, eliminate rho and epsilon first.
    rho_log_prime=-v*W2
    T_log_prime=(-(J+Q)/w-(1-e)*rho_log_prime)/b
    pressure_prime=w*(r*rho_log_prime+d*T_log_prime)
    p_at_root=sp.factor(pressure_prime.subs(J,w*W2*(v+2*j)-Q))
    assert sp.factor(p_at_root+w*W2*(sound*v+2*j*d/b))==0
    fprime=w*W2*(1+2*j*v)+v*p_at_root
    determinant=1-sound*v*v+2*j*v*(1-d/b)
    assert sp.factor(fprime-w*W2*determinant)==0
    assert sp.factor(sp.diff(determinant,v,2)+2*sound)==0
    # Recover the point-inverse determinant already derived independently.
    raw=b-v*v*(b*r+d*e)+2*j*v*(b-d)
    assert sp.factor(raw/b-determinant)==0
    # The square-root relation uses the positive W branch, without subtracting rest energies.
    s=sp.symbols('s',positive=True)
    kinetic=sp.factor((1-s)-v*v/(1+s))
    assert sp.factor(kinetic.subs(v*v,1-s*s))==0
    # Exact local-jet negative control: causal equilibrium sound speed and DEC alone are insufficient.
    control={b:sp.Rational(1,100),r:sp.Rational(1,100),d:sp.Rational(3,100),e:sp.Rational(3,100),j:sp.Rational(3,10)}
    assert sound.subs(control)==sp.Rational(1,10)
    vstar=sp.sqrt(46)-6;assert 0<vstar<1
    assert sp.simplify(determinant.subs(control).subs(v,vstar))==0
    assert determinant.subs(control).subs(v,sp.Rational(4,5))==-sp.Rational(3,125)
    eps=sp.Integer(1);P=sp.Rational(1,10);ww=eps+P;QQ=sp.Rational(3,10)*ww
    radical=sp.sqrt(ww*ww-4*QQ*QQ);energy=(eps-P+radical)/2;radial=(-eps+P+radical)/2
    assert energy==sp.Rational(89,100) and radial==-sp.Rational(1,100) and energy>abs(radial) and energy>P
    rho,T=sp.symbols('rho T',positive=True)
    free=sp.Rational(9,10)+(rho-1)/10-(T-1)/10-sp.Rational(189,2000)*(rho-1)**2+sp.Rational(33,1000)*(rho-1)*(T-1)-sp.Rational(11,2000)*(T-1)**2
    epsilon=rho*(free-T*sp.diff(free,T));pressure=rho*rho*sp.diff(free,rho);enthalpy=epsilon+pressure;point={rho:1,T:1}
    assert epsilon.subs(point)==eps and pressure.subs(point)==P
    realized=[T*sp.diff(epsilon,T)/enthalpy,rho*sp.diff(pressure,rho)/enthalpy,T*sp.diff(pressure,T)/enthalpy,1-rho*sp.diff(epsilon,rho)/enthalpy]
    assert [sp.simplify(x.subs(point)) for x in realized]==[control[x] for x in [b,r,d,e]]
    save('symbolic.json',dict(classification='Proven',passed=True,
        reduction='rho(v)=D*sqrt(1-v^2), epsilon(v)=E-v*(J+Q). On a connected interval where epsilon(rho(v),T)=epsilon(v) has a positive C1 thermal solution, define F(v)=v*(E+P(rho(v),T(v)))-(J-Q). Its zeros are exactly the fixed-Q primitive solutions.',
        derivative='At a zero, E+P=w*W^2*(1+2*j*v), J+Q=w*W^2*(v+2*j), dP/dv=-w*W^2*(g*v+2*j*d/b). Therefore F_prime=w*W^2*[1-g*v^2+2*j*v*(1-d/b)]. The bracket equals the prior three-variable Jacobian determinant divided by its positive rho*w^2*W^5*b prefactor.',
        uniqueness='If the bracket is positive at every zero on the connected admissible interval, every zero crosses from negative to positive. Two such zeros would require an intervening zero that does not cross that way, contradicting the positive derivative there. Thus there is at most one zero. Continuity on a closed bracket inside that interval with F(left)<0<F(right) supplies existence by IVT and hence a unique solution. No disconnected thermal branches are silently joined.',
        frozen_extremum='For frozen g>=0 the normalized determinant is concave in v. Its exact minimum on[-1,1] is the smaller endpoint, 1-g-2*abs(j)*abs(1-d/b). A positive value certifies every subluminal velocity for that frozen polynomial only. In a physical branch the coefficients vary with rho,T and the premise must be bounded over that branch.',
        rest_energy='For E_star=E-D*C, the thermal target can be evaluated as epsilon-rho*C=E_star+D*C*v^2/(1+sqrt(1-v^2))-v*(J+Q). This identity avoids the direct E-rho*C subtraction but does not restore lost information in rounded conserved data.',
        negative_control='b=r=1/100,d=e=3/100,j=3/10 gives g=1/10 and normalized determinant1-v^2/10-6*v/5, singular at v=sqrt(46)-6. With epsilon=1,P=1/10,Q=33/100, the Landau energy89/100 exceeds both radial pressure magnitude1/100 and transverse pressure1/10. These EOS jets are realized exactly at rho=T=1 by the saved polynomial specific Helmholtz free energy, with epsilon=rho*(f-T*f_T), P=rho^2*f_rho. Thus positive heat capacity, positive isothermal stiffness, Maxwell e=d, subluminal equilibrium sound and the algebraic dominant-energy condition do not suffice for invertibility at every velocity. This is a consistent local thermodynamic jet counterexample, not a calibrated stellar EOS or a claim that the complete dissipative heat system is causal there.',
        negative_control_specific_free_energy=str(free),
        unresolved='The actual EOS derivative/error bounds over an admissible connected branch, a certified sign bracket for evolved conserved data, composition/reaction changes, cell-moment closure and the evolved heat variable remain necessary.'))


def run():
    p=bindings();symbolic();fields=np.load(heat.old.char.OUT/'coefficients.npz')['coefficients'];heats=np.load(heat.OUT/'initial-rates.npz')['j']
    assert len(fields)==len(heats)==p['cells'];counts=dict(positive=0,zero=0,negative=0);minimum=None;limiting=None
    with gzip.open(OUT/'certificates.jsonl.gz','xt',compresslevel=1) as stream:
        for i,(row,flux) in enumerate(zip(fields,heats,strict=True)):
            B,R,DD,EE,H=map(lambda x:F(float(x)),row);JJ=F(float(flux));assert B>0
            sound=R+DD*EE/B;assert 0<=sound<1
            left=1-sound-2*JJ*(1-DD/B);right=1-sound+2*JJ*(1-DD/B);margin=min(left,right)
            assert margin==1-sound-2*abs(JJ)*abs(1-DD/B)
            verdict='positive' if margin>0 else 'negative' if margin<0 else 'zero';counts[verdict]+=1
            if minimum is None or margin<minimum:minimum=margin;limiting=i
            stream.write(json.dumps(dict(cell=i,sound_squared=str(sound),normalized_endpoint_bounds=[str(left),str(right)],minimum_normalized_determinant=str(margin),verdict=verdict))+'\n')
    prior=json.loads((previous.OUT/'result.json').read_text());assert prior['all_passed'] and len(prior['rows'])==27
    save('result.json',dict(classification='Proven',audit_complete=True,cells=p['cells'],counts=counts,minimum_normalized_determinant=str(minimum),limiting_cell=limiting,
        every_frozen_velocity_positive=counts['positive']==p['cells'],old_native_control_cases=27,
        actual_EOS_branch_certified=False,global_actual_native_inverse_certified=False,unresolved_cell_average_closure=False,finite_GR_evolution=False))
    save('manifest.json',dict(sha256={f.relative_to(g.ROOT).as_posix():sha(f) for f in OUT.iterdir() if f.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['audit_complete'] and sum(r['counts'].values())==p['cells']
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    counts=dict(positive=0,zero=0,negative=0);minimum=None
    with gzip.open(OUT/'certificates.jsonl.gz','rt') as stream:
        for i,line in enumerate(stream):
            row=json.loads(line);assert row['cell']==i;lo=min(map(F,row['normalized_endpoint_bounds']));assert lo==F(row['minimum_normalized_determinant'])
            verdict='positive' if lo>0 else 'negative' if lo<0 else 'zero';assert row['verdict']==verdict;counts[verdict]+=1
            minimum=lo if minimum is None else min(minimum,lo)
    assert counts==r['counts'] and minimum==F(r['minimum_normalized_determinant'])
    print('PASS connected thermal-branch uniqueness theorem and full-velocity frozen-jet audit:',counts,flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

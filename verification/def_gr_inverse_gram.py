"""Square-root coupled resolvent and bounded original heat-pole moments.

R=sigma*(L^-1 D M Q).T*(L^-1 D M Q), with D(K+sigma M)D=L L.T.
The same factor transforms the source. No cancellation-prone T.T*A^-1*T.
"""
from pathlib import Path
import argparse
import time
import json
import signal
import resource
import numpy as np
from scipy.linalg import cholesky_banded,eig,solve
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve_triangular
import def_gr_multishift as multi

task=multi.task;space=multi.space;common=multi.common;fields=multi.fields
OUT=common.OUT/'inverse-gram';write=task.write


def basis(model,size):
    _,_,Q,scales,checks=multi.basis(model,size);sigma=128.**2
    K=(model.K+model.K.T)/2;M=(model.M+model.M.T)/2
    A=(K+sigma*M).tocsc();s=1/np.sqrt(A.diagonal());D=diags(s);scaled=(D@A@D).tocsc()
    coo=scaled.tocoo();width=int(max(abs(coo.row-coo.col)));n=model.size;band=np.zeros((width+1,n))
    for j in range(width+1):band[j,:n-j]=scaled.diagonal(-j)
    cb=cholesky_banded(band,lower=True)
    L=diags([cb[j,:n-j] for j in range(width+1)],-np.arange(width+1),shape=(n,n)).tocsr()
    # Recompute the same lifted load using the symmetric stiffness; F is kept.
    load=(model.load+(K-model.K)@model.H).tocsc();T=M@Q
    B=spsolve_triangular(L,s[:,None]*T,lower=True)
    R=sigma*(B.T@B);skew=float(np.max(abs(R-R.T))/max(np.max(abs(R)),1e-100));assert skew<1e-12
    mass=float(np.max(abs(Q.T@(M@Q)-np.eye(size))));assert mass<1e-10
    heat=model.heat;tc=model.original.radiation.geometry.tc;lam=heat.rates*tc;moments=np.zeros((len(heat.edges),8))
    for p in range(8):moments[heat.face_ids,p]=np.sum(heat.amplitude*tc/lam**p,axis=1)
    F=load@moments[:,:6];Z=spsolve_triangular(L,s[:,None]*F,lower=True)
    projected=sigma*(B.T@Z)
    adjoint=s[:,None]*spsolve_triangular(L.T.tocsr(),B,lower=False)
    remainder=abs(load)@abs(moments[:,6:8])
    # A small algebraic residual is checked, without mistaking it for a forward
    # error certificate. Gram symmetry alone likewise cannot accept a trajectory.
    error=float(np.max(abs(A@adjoint-T)/(abs(A)@abs(adjoint)+abs(T)+1e-100)))
    assert error<1e-9,error
    return R,projected,Q,scales,(adjoint,remainder,float(lam.min())),dict(**checks,resolvent_symmetry=skew,
        full_mass_orthogonality=mass,linear_residual=error,symmetric_assembly_relative_change=float(np.max(abs((K-model.K).data),initial=0)/np.max(abs(model.K.data))))


def response(model,R,F,packed):
    Q,scales,tail=packed;balanced=R*scales[None,:]/scales[:,None];theta,V=eig(balanced)
    assert np.all(theta.imag==0) and theta.real.min()>0;theta=theta.real
    k=128.**2*(1-theta)/theta;assert k.min()>-1.
    coefficients=solve(V,F/scales[:,None]).real/theta[:,None];physical=(Q@(scales[:,None]*V)).real
    eigen_error=float(np.max(abs(balanced@V-V*theta))/max(np.max(abs(balanced)@abs(V)),1e-100));assert eigen_error<1e-10
    adjoint,remainder,lam_min=tail
    factor=128.**2*(solve(V,np.diag(1/scales)).real/theta[:,None])
    modal_adjoint=adjoint@factor.T
    absolute_moments=abs(modal_adjoint).T@remainder
    ratio=abs(k)/lam_min**2;assert ratio.max()<1e-12
    tails=k[:,None]**2/(1-ratio[:,None])*absolute_moments
    values=np.zeros((65,len(k)));velocities=values.copy()
    for j,t in enumerate(common.modes.TIMES[1:],1):
        for i,ki in enumerate(k):
            S,C,T,co=common.modes.functions(ki,t);a,b,c,d,e,f=coefficients[i]
            c-=ki*e;d-=ki*f
            values[j,i]=a*T-b*C+c*S-d*co
            velocities[j,i]=a*C-b*S+c*co+ki*d*S
    r=model.original.native;N,a=model.original.radiation.geometry.metric(r)
    speed=a/N*r*task.fem.base.task.h.gr.C;w=model.original.weights
    qb=abs(physical)@(np.sinh(1)*tails[:,0]+np.cosh(1)*tails[:,1])
    vb=abs(physical)@(np.cosh(1)*tails[:,0]+np.sinh(1)*abs(k)*tails[:,1])
    cv=abs(speed)*(abs(model.nativeV[0])@vb);cf=abs(model.nativeV[1])@qb
    bounds=dict(velocity_mass_RMS_m_s=float(np.sqrt(w@(cv*cv))),scalar_mass_RMS=float(np.sqrt(w@(cf*cf))),
        **{key:float(np.sqrt(w[mask]@cv[mask]**2/w[mask].sum())) for key,mask in zip(task.FIELDS[2:],model.original.masks)})
    assert lam_min*common.modes.TIMES[1]>700
    return physical@values.T,physical@velocities.T,dict(minimum_stiffness=float(k.min()),maximum_stiffness=float(k.max()),
        resolvent_eigenvalue_minimum=float(theta.min()),balanced_eigen_residual=eigen_error,
        heat_pole_ratio_maximum=float(ratio.max()),pole_expansion_absolute_bounds=bounds,
        omitted_exponential_coefficient_factor=float(np.exp(-700)))


def series(model,R,F,Q,scales,tail,label):
    class Packed(tuple):
        @property
        def shape(self):return self[0].shape
    common.OUT=OUT;common.response=response
    row=common.series(model,R,F,Packed((Q,scales,tail)),label)
    for key,bound in row['pole_expansion_absolute_bounds'].items():assert bound<1e-12*max(abs(r[key]) for r in row['history'])
    return row


def control():
    import sympy as s
    lam,k=s.symbols('lam k',positive=True)
    assert s.simplify(1/(lam**2+k)-(1/lam**2-k/lam**4)-k*k/(lam**4*(lam**2+k)))==0
    matrix=np.array([[9.,1e-11],[1e-11,25.]])+128.**2*np.eye(2);L=np.linalg.cholesky(matrix)
    B=np.linalg.solve(L,np.eye(2));R=128.**2*B.T@B
    force=np.array([1.,2e-12]);f=128.**2*B.T@np.linalg.solve(L,force)
    errors=[]
    for z in [2+3j,6+100j]:
        actual=np.linalg.solve(z*z*R+128.**2*(np.eye(2)-R),f)
        expected=np.linalg.solve(matrix+(z*z-128.**2)*np.eye(2),force)
        errors.append(float(np.max(abs(actual-expected)/abs(expected))))
    assert max(errors)<1e-9
    return dict(classification='Proven',gram_identity='A^-1=L^-T*L^-1; apply the same factor to both trial and source.',
        pole_remainder_identity_checked=True,numerical_classification='Counterexample candidate',weak_field_transfer_errors=errors,
        scope='Analytic identities; evaluated tail bounds and full GR acceptance remain numerical and conditional, not interval certification.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Preserve inverse reciprocity through a Cholesky Gram and source transform, while retaining field separation, input-norm balance and the fixed broadband basis.',
        method='Same degree4 and total128/256/512 columns; symmetric roundoff-level weak assembly. Cholesky Gram of the full coupled shifted matrix. Original heat poles use a two-term geometric expansion with an explicit absolute remainder below1e-12 of each readout, not a fitted input.',
        gates=dict(propagation_relative=.02,propagation_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,inverse_symmetry=1e-12,heat_balance=2e-13,pole_tail_relative=1e-12),
        decision='First fixed64 pilot, then degree4 only if symmetry, source-tail and budget gates pass. Stop on any actual component failure. Only if accepted run fixed degrees2,1 and original conditional contrasts; no basis or resolution expansion.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=120,total_cap_seconds=300,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(multi.__file__),Path(fields.__file__),Path(common.__file__),Path(space.__file__),Path(task.__file__),OUT.parent/'reciprocity-refinement.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');setup=time.monotonic()-start
    start=time.monotonic();R,F,Q,S,tail,meta=basis(model,64);build=time.monotonic()-start;row=series(model,R,F,Q,S,tail,'pilot-64')
    forecast=1.4*(setup+64*build+8*row['seconds'])
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,build64_seconds=build,
        response64_seconds=row['seconds'],first_case_forecast_seconds=forecast,
        assumption='Conservative quadratic build and linear response from64 to512;40 percent allowance. Original first-case120s and total300s cap.',**meta))
    signal.alarm(0);print('GRAM FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<120
    signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');R,F,Q,S,tail,meta=basis(model,512)
    rows=[series(model,R[:n,:n],F[:n],Q[:,:n],S[:n],(tail[0][:,:n],tail[1],tail[2]),f'p4-{n}') for n in [128,256,512]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False,**meta)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('GRAM RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()

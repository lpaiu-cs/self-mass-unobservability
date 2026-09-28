"""Bounded same-state EOS and common optical-level checks, no stellar evolution."""
from types import FunctionType
import argparse
import ctypes
import json
import signal
import time
import numpy as np
import sympy as sp
import def_photon_shared_atomic as a


def run(reuse=False):
    out=a.OUT;assert not (out/'result.json').exists();signal.alarm(60);begin=time.monotonic()
    assert np.finfo(np.longdouble).nmant>=63,'Run the optical pair check in the declared WSL extended-precision runtime'
    prior=a.p.a.previous;gas_module=prior.old.previous.matter.old
    gas=object.__new__(gas_module.GasEOS);init=gas_module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=a.CACHE/'gas.so'),closure=init.__closure__)(gas)
    gas.inventory_lib=gas.gas_lib;native=gas.gas_lib.ionization_inventory
    def call(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=call;lib=gas.gas_lib
    def arr(name,n,dtype=ctypes.c_double,module='mod_excitation_block'):
        return np.ctypeslib.as_array((dtype*n).in_dll(lib,'__'+module+'_MOD_'+name)).copy()
    d,_=prior.base.inputs();r=float(d['lnd'][0]);X=d['X'][0]
    t=float(np.log(np.load(a.p.OUT/'eos-state.npz')['T'][0]));T=float(np.exp(t));snaps=[];checks=[]
    def snapshot(dr,dt):
        snap=prior.inventory_reader.InventoryEOS.snapshot(gas,r+dr,t+dt,X)
        checks.append(prior.inventory_reader.check(snap,X,r+dr,gas))
        count=int(arr('extrace_count',1,ctypes.c_int)[0])
        snap.update(ids=arr('extrace_ids',636,ctypes.c_int).reshape(318,2)[:count],
            value=arr('extrace_value',318*6).reshape(318,6)[:count],scale=arr('extrace_scale',318)[:count],
            x=arr('x',5),ground_logw=arr('shared_ground_logw',318,module='mod_excitation'))
        assert all(np.all(np.isfinite(v)) for v in snap.values())
        snaps.append(snap);return snap
    if reuse:
        saved=np.load(out/'eos-state.npz');base={k:saved[k][0] for k in saved.files};v=base['eos']
        previous=json.loads((out/'result-first-double-precision.json').read_text())
        derivatives=previous['derivatives'];checks=previous['checks'];history=previous['history_bitwise']
        inverse_error=previous['pressure_inverse_log_density_error']
    else:
        base=snapshot(0,0);v=base['eos'];derivatives=[]
        for h in [2e-4,1e-4]:
            rm=snapshot(-h,0)['eos'];rp=snapshot(h,0)['eos']
            tm=snapshot(0,-h)['eos'];tp=snapshot(0,h)['eos']
            fd=np.array([(np.log(rp[1])-np.log(rm[1]))/(2*h),
                (np.log(tp[1])-np.log(tm[1]))/(2*h),(rp[2]-rm[2])/(2*h),(tp[2]-tm[2])/(2*h)])
            expected=v[[5,6,9,10]];err=abs(fd-expected)/np.maximum(abs(expected),1.)
            # F=u-Ts; differentiating the equilibrated F tests both EOS hooks.
            fr=((rp[2]-T*rp[3])-(rm[2]-T*rm[3]))/(2*h)
            ft=((tp[2]-T*np.exp(h)*tp[3])-(tm[2]-T*np.exp(-h)*tm[3]))/(2*h)
            ferr=abs(np.array([fr,ft])-np.array([v[1]/v[0],-T*v[3]]))/np.maximum(np.abs([v[1]/v[0],T*v[3]]),1.)
            derivatives.append(dict(h=h,relative_errors=err.tolist(),free_energy_relative=ferr.tolist()))
        repeated=snapshot(0,0);history=all(np.array_equal(base[k],repeated[k]) for k in base)
        inverse=gas(1,float(np.log(v[1])),t,X);inverse_error=float(abs(np.log(inverse[0])-r))
        base_again=snapshot(0,0);assert all(np.array_equal(base[k],base_again[k]) for k in base)
        np.savez_compressed(out/'eos-state.npz',**{k:np.asarray([s[k] for s in snaps]) for k in base})
    catalog=json.loads((out/'catalog.json').read_text());erg,ryd,c2,clight,kb=catalog['constants']
    terms=lib.shared_terms;array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
    terms.argtypes=[ctypes.c_int,ctypes.c_double,array,array];terms.restype=None
    rows=[];optical={}
    x=np.ascontiguousarray(base['x'])
    for k,row in enumerate(catalog['rows']):
        i=row['ion_index'];mask=np.array(row['included']);g=np.array(row['weight'],dtype=np.longdouble)[mask]
        bind=np.array(row['binding_cm_inverse'],dtype=np.longdouble)[mask]
        # One energy representation, extended precision for near-degenerate pairs.
        e=np.longdouble(row['ground_cm_inverse'])-bind;ct=np.longdouble(c2)/np.longdouble(T)
        def evaluate(tval,xval=x):
            output=np.zeros(38*700);terms(i,tval,np.ascontiguousarray(xval),output)
            return output.reshape(700,38)[:len(g)-1].copy()
        q=evaluate(t);h=1e-4;qm=evaluate(t-h);qp=evaluate(t+h)
        partial1=np.max(abs((qp[:,0]-qm[:,0])/(2*h)-q[:,1])/np.maximum(abs(q[:,1]),1e-30))
        partial2=np.max(abs((qp[:,0]-2*q[:,0]+qm[:,0])/h**2-q[:,2])/np.maximum(abs(q[:,2]),1e-30))
        # Optical and EOS sums independently meet through the emitted ground factor.
        qld=q[:,0].astype(np.longdouble)
        ratios=qld*np.exp(-ct*np.longdouble(row['ground_cm_inverse'])-np.longdouble(base['ground_logw'][i-1]))/g[0]
        fraction=np.r_[np.longdouble(1),ratios];fraction/=fraction.sum()
        match=np.flatnonzero(np.all(base['ids']==[row['element_index'],row['charge']],axis=1));assert len(match)==1
        j=int(match[0]);L=base['value'][j,3]/base['scale'][j]
        partition_error=float(abs(np.log1p(ratios.sum())-L)/max(abs(L),1e-30))
        inventory=base['number_fractions'][row['element_index']-1,row['charge']]
        populations=inventory*fraction
        # Relative survival weights include the exact *native* ground correction.
        logw=np.r_[np.longdouble(base['ground_logw'][i-1]),np.log(qld/g[1:])-ct*bind[1:]]
        errors=[];old_errors=[];pair_rows=[]
        for lo in range(len(g)):
            for hi in range(lo+1,len(g)):
                delta=ct*(bind[lo]-bind[hi])
                if delta<=1e-12:continue  # Degenerate pairs have no positive-frequency photon.
                # Declared common-U survival: p_ij=min(w_i,w_j). Both conditional
                # probabilities stay <=1 even when occupations are not ordered.
                ratio=np.exp(logw[hi]-logw[lo]);up_factor=min(np.longdouble(1),ratio)
                down_factor=min(np.longdouble(1),1/ratio)
                up=fraction[lo]*g[hi]/g[lo]*up_factor;down=fraction[hi]*down_factor
                absorption=up-down;emission=down*np.expm1(delta)
                assert 0<up_factor<=1 and 0<down_factor<=1 and absorption>0
                errors.append(abs(absorption-emission)/max(abs(absorption),abs(emission),1e-300))
                ordinary=fraction[lo]*g[hi]/g[lo]-fraction[hi]
                ordinary_emission=fraction[hi]*np.expm1(delta)
                old_errors.append(abs(ordinary-ordinary_emission)/max(abs(ordinary),abs(ordinary_emission),1e-300))
                pair_rows.append([lo,hi,up_factor,down_factor,absorption,down,ratio])
        row_result=dict(Z=row['Z'],charge=row['charge'],levels=len(g),excited_fraction=float(1-fraction[0]),
            partition_relative=partition_error,inventory_relative=float(abs(populations.sum()/inventory-1)),
            positive_frequency_pairs=len(errors),kirchhoff_relative=float(max(errors)),
            unmodified_Einstein_relative=float(max(old_errors)),first_partial_relative=float(partial1),second_partial_relative=float(partial2),
            max_upward_occupation_ratio=float(max(v[6] for v in pair_rows)),conditional_probabilities_bounded=True)
        rows.append(row_result);optical.update({f'{key}_{k}':value for key,value in dict(fraction=fraction,populations=populations,logw=logw,weight=g,excitation_cm_inverse=e,level_terms=q,pairs=np.asarray(pair_rows,dtype=np.longdouble)).items()})
    np.savez_compressed(out/'optical-levels.npz',**optical)
    old=np.load(a.p.OUT/'eos-state.npz')['eos'][0]
    changes={name:dict(old=float(old[i]),candidate=float(v[i]),relative=float(v[i]/old[i]-1)) for name,i in [('pressure',1),('energy',2),('entropy',3),('electron_density_over_NA',13),('heat_capacity_times_T',10)]}
    passed=bool(history and inverse_error<1e-10 and all(c['inventory_error']<1e-10 and c['charge_error']<1e-10 for c in checks)
        and all(max(z['relative_errors']+z['free_energy_relative'])<1e-5 for z in derivatives)
        and all(z['partition_relative']<1e-12 and z['inventory_relative']<1e-12 and z['kirchhoff_relative']<1e-10
            and max(z['first_partial_relative'],z['second_partial_relative'])<1e-5 for z in rows))
    result=dict(classification='Counterexample candidate',passed=passed,seconds=time.monotonic()-begin,native_EOS_calls=0 if reuse else 12,
        reused_native_EOS_calls=12 if reuse else 0,longdouble_precision_bits=int(np.finfo(np.longdouble).nmant),
        bindings={str(f):a.p.a.digest(f) for f in [out/'eos-state.npz',out/'build.json',out/'catalog.json',a.LIB,a.CACHE/'gas.so']},
        checks=checks,derivatives=derivatives,history_bitwise=history,pressure_inverse_log_density_error=inverse_error,
        changes=changes,atomic_stages=rows,new_stellar_steps=0,coupled_runs=0,physical_opacity_certified=False,
        boundary='A common finite PL/MHD state model and same-state thermodynamic/optical provider. No cross-section amplitude, missing-stage or dissolved-continuum certification, and no replacement of the saved stellar history.')
    result['joint_survival_model']='p_ij=min(w_i,w_j), realized by a common uniform survival variable; a declared assumption, not a derivation of physical level correlations.'
    a.p.a.write(out/'result.json',result);print('SHARED ATOMIC',passed,'seconds',result['seconds'],
        'derivative',max(max(x['relative_errors']) for x in derivatives),'partition',max(x['partition_relative'] for x in rows),
        'Kirchhoff',max(x['kirchhoff_relative'] for x in rows),'levels',sum(x['levels'] for x in rows),flush=True)
    n,g,wi,wj,joint,delta=sp.symbols('n g wi wj joint delta',positive=True)
    upper=n*g*wj/wi*sp.exp(-delta)
    assert sp.simplify(n*g*joint/wi-upper*joint/wj-upper*joint/wj*(sp.exp(delta)-1))==0
    a.p.a.write(out/'symbolic.json',dict(classification='Proven',passed=True,
        identity='With populations proportional to g_i w_i exp(-E_i/kT), any common joint survival p_ij gives conditional factors p_ij/w_i and p_ij/w_j and Kirchhoff balance. For p_ij=min(w_i,w_j) both factors are bounded by one. Degenerate pairs have no positive-frequency photon.',
        limitation='Conditional algebra; thermodynamic PL/MHD survival and the common-uniform-variable joint model are Conjectural.'))
    signal.alarm(0);assert passed,'same-state shared atomic gates failed'


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--reuse-state',action='store_true')
    run(parser.parse_args().reuse_state)

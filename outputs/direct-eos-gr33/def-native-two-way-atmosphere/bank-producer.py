"""Recover arbitrary-spectrum H coefficients from saved native constraint roots."""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import def_native_hydrogen_exchange as ex
import def_native_causal_photons as photons

OUT=ex.OUT.parent/'def-native-two-way-atmosphere'
write=ex.write
sha=photons.sha


def timeout(*_):raise TimeoutError('Registered phase112 wall-time cap')


def metadata(native,raw,y,lt):
    lib=native.ion.gas.gas_lib
    def arr(name,n,dtype=ctypes.c_double):return np.ctypeslib.as_array((dtype*n).in_dll(lib,'__mod_excitation_block_MOD_'+name)).copy()
    count=int(arr('extrace_count',1,ctypes.c_int)[0]);ids=arr('extrace_ids',636,ctypes.c_int).reshape(318,2)[:count]
    row=np.flatnonzero(np.all(ids==[1,0],axis=1));assert len(row)==1
    L=float((arr('extrace_value',318*6).reshape(318,6)[:count,3]/arr('extrace_scale',318)[:count])[row[0]])
    terms=np.zeros(30);native.levels(1,float(lt),arr('x',5),terms);tail=terms.reshape(3,10)[0];terms=tail-np.r_[tail[1:],0.]
    assert np.all(terms>=0) and L>=0
    fraction=np.zeros(10);fraction[0]=np.exp(-L)
    if L>0:fraction[1:]=-np.expm1(-L)*terms[1:]/sum(terms[1:])
    return dict(raw=raw,fraction=fraction,affinity=float(native.ion.fields[0]),y=y,lt=lt)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='2029784c4',
        claim='Replace the prescribed atmospheric Planck bath by the same time-dependent angular photons that evolve in the deep native interior, with both directions exchanged at the actual radius.',
        decision='Evolve the moving native interface using actual absorption/emission and paired radiation energy/momentum. Only accepted time/space controls may later supply the final scalar source.',
        reuse='Retain Phase111 deep EOS and its unsplit material-photon integrator, Phase108 physical fluid initialization/recovery/fluxes, and the918 Phase107 constrained atmosphere states. Restore saved native roots once to recover missing spectral level metadata; do not redo the constraint search.',
        geometry='Deep last volume ends at the actual atmosphere inner face; preserve its actual native quadrature point. Remove overlap. Propagate fine atmospheric photons with the same spherical Killing-frequency measure and one shared face flux.',
        method='Keep deep transport/chemistry unsplit. Transparent atmospheric inward/outward transport is eliminated in angular direction order; local atmospheric collisions and conservative moving gas use SSP substeps. Judge this partition on actual coupled controls, not by old separate passes.',
        initial='Keep the registered interior Planck mean and grey flux; extend surface outgoing occupations into initially empty exterior along the declared free-streaming geometry. Incoming exterior is vacuum.',
        gates=dict(native=.002,energy=1e-8,baryon=1e-10,species=1e-9,positive_photons=True,time_trace=.02,space_trace=.02),
        budget=dict(cached_native_replays=918,new_native_controls_call_cap=100,total_native_call_cap=1100,spectral_seconds=60,
                    coupling_pilot_seconds=45,production_seconds=600,CPU_threads=1,memory_GB=3),
        stop='No unmeasured dispatch of all production paths, no automatic added grid/horizon/EOS support or weaker gate. Measure the coupled pilot and require the forecast to fit the declared cap. Preserve failures; final charge remains open.',
        bindings={str(p):sha(p) for p in [Path(__file__),ex.OUT/'repaired-bank.npz',ex.OUT/'bank-native-states.npz',ex.OUT/'repair-native-states.npz',
             Path(ex.__file__),ex.ATOMIC,ex.LEVELS,OUT.parent/'def-native-unsplit-photons/result.json',OUT.parent/'def-native-reactive-interface/physical/result.json']}))


def bank():
    assert not (OUT/'spectrum-bank.npz').exists();start=time.monotonic();signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    native=ex.Native(cap=1100);d=np.load(ex.OUT/'repaired-bank.npz');saved=[np.load(ex.OUT/f) for f in ['bank-native-states.npz','repair-native-states.npz']]
    cache={k:np.concatenate([v[k] for v in saved]) for k in saved[0].files};hy=cache['number_fractions'][:,0,0]/cache['number_fractions'][:,0,:2].sum(1)
    frac=np.zeros((*d['raw'].shape[:3],10));affinity=np.zeros(d['raw'].shape[:3]);done=np.zeros_like(affinity,bool);maximum=0.
    try:
        for iy,y in enumerate(d['ys']):
            for ix,x in enumerate(d['x']):
                for it,lt in enumerate(d['lt']):
                    ids=np.flatnonzero((abs(cache['lrho']-native.lr-x)<1e-10)&(abs(cache['logT']-lt)<1e-12)&(abs(hy/y-1)<1e-7));assert len(ids),('Cached root missing',iy,ix,it)
                    j=ids[-1];state=native.ion.snapshot(float(native.lr+x),float(lt),cache['fields'][j].copy())
                    error=float(np.max(abs(state['eos'][[0,1,2]]/d['raw'][iy,ix,it,[0,1,2]]-1)));maximum=max(maximum,error);assert error<1e-9,('Cached replay',iy,ix,it,error)
                    meta=metadata(native,state['eos'],float(y),float(lt));frac[iy,ix,it]=meta['fraction'];affinity[iy,ix,it]=meta['affinity'];done[iy,ix,it]=True
        np.savez_compressed(OUT/'spectrum-bank.npz',x=d['x'],lt=d['lt'],ys=d['ys'],fraction=frac,affinity=affinity,binding=native.binding,y0=d['y0'],nH=d['nH'],rho0=d['rho0'])
        result=dict(classification='Counterexample candidate',passed=True,restored_states=int(done.sum()),native_calls=native.ion.calls,
            maximum_constitutive_replay_relative=maximum,seconds=time.monotonic()-start,source_sha256=sha(__file__))
        write(OUT/'spectrum-bank.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        np.savez_compressed(OUT/'spectrum-bank-partial.npz',fraction=frac,affinity=affinity,done=done)
        write(OUT/'spectrum-bank-failure.json',dict(error=repr(exc),completed=int(done.sum()),native_calls=native.ion.calls,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


class Spectrum:
    def __init__(self):
        self.d=d=np.load(OUT/'spectrum-bank.npz');self.binding=d['binding'];self.native=ex.Native(cap=100)
        self.rho0=float(d['rho0']);self.ys=d['ys'];self.nH=float(d['nH']);self.f=[];self.rev=[];self.aff=[]
        delta=self.binding[0]-self.binding
        for k,y in enumerate(self.ys):
            frac=d['fraction'][k];T=np.exp(d['lt'])[None,:,None]
            # Known Boltzmann factors are evaluated directly. If a retained
            # native level is zero at a tabulated endpoint, interpolate its
            # untransformed fraction locally instead of inventing a prefactor.
            pref=frac*np.exp(np.minimum(delta[None,None,:]/(ex.K*T),700))
            self.f.append((RegularGridInterpolator((d['x'],d['lt']),pref),RegularGridInterpolator((d['x'],d['lt']),frac),RegularGridInterpolator((d['x'],d['lt']),(frac>0).astype(float))))
            reverse=np.log(y/(1-y))+d['affinity'][k,:, :,None]-self.binding[None,None,:]/(ex.K*T)
            self.rev.append(RegularGridInterpolator((d['x'],d['lt']),reverse))

    def levels(self,rho,lt,y):
        x=np.log(rho/self.rho0);query=np.column_stack([x,lt]);w=(y-self.ys[0])/np.diff(self.ys)[0];T=np.exp(lt)
        assert np.all((w>=-1e-8)&(w<=1+1e-8))
        fractions=[];reverse=[]
        for (pref,raw,mask),rr in zip(self.f,self.rev):
            f=pref(query)*np.exp(-(self.binding[0]-self.binding)[None,:]/(ex.K*T[:,None]));f=np.where(mask(query)>=1-1e-14,f,raw(query))
            f/=f.sum(1)[:,None];fractions.append(f);reverse.append(f*np.exp(rr(query)))
        f=fractions[0]*(1-w[:,None])+fractions[1]*w[:,None]
        r=reverse[0]*(1-w[:,None])+reverse[1]*w[:,None]
        return f*y[:,None],r*(1-y[:,None])

    def cross(self,energy,level):
        good=energy>=self.binding[level-1];out=np.zeros_like(energy)
        nu=np.ascontiguousarray((self.native.oldground/level**2+energy[good]-self.binding[level-1])/float(np.float32(6.6256e-27)))
        bf=np.zeros_like(nu);ff=np.zeros_like(nu)
        if len(nu):self.native.cross(len(nu),1,level,1.e4,nu,bf,ff)
        out[good]=bf;return out

    def coefficients(self,rho,lt,y,energy):
        f,r=self.levels(rho,lt,y);ab=np.zeros_like(energy);em=np.zeros_like(energy);T=np.exp(lt)
        for k in range(10):
            sigma=self.cross(energy,k+1);ab+=f[:,k,None]*sigma
            em+=r[:,k,None]*sigma*np.exp(-np.maximum(energy-self.binding[k],0)/(ex.K*T[:,None]))
        return ab,em


if __name__=='__main__':globals()[sys.argv[1]]()

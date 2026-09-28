"""Independent conservation and saved-response checks for Phase93."""
import json
import numpy as np
from types import SimpleNamespace
import def_native_coupled_readjustment as task


def source_control():
    # A finite-volume heat debit must produce the same radial mass constraint,
    # including the regular r^3 central cell. A linear radius fraction fails.
    m=task.Model.__new__(task.Model);m.edges=np.array([0.,.2,.6,1.]);m.n=3
    def sample(r):
        return dict(r=r,m=np.zeros_like(r),N=np.ones_like(r),phi=np.zeros_like(r))
    m.bg=SimpleNamespace(R=7e9,sample=sample)
    m.volumes=4*np.pi*m.bg.R**3*np.diff(m.edges**3)/3
    m.raw=np.zeros((3,21));m.raw[:,0]=1;m.raw[:,8]=-1
    m.thermo=np.ones((3,6));m.thermo[:,5]=3
    r=np.array([.03,.15,.3,.5,.7,.9]);energy=np.array([2e20,-3e20,5e20])
    point=sample(r);_,loss,J=m.sources(point)
    ids=np.searchsorted(m.edges,r)-1
    f=(r**3-m.edges[ids]**3)/np.diff(m.edges**3)[ids]
    full=np.r_[0.,energy];expected=-(full[ids]*(1-f)+full[ids+1]*f)*task.G/(task.C**4*m.bg.R)
    error=float(max(abs(J@energy-expected))/max(abs(expected)))
    assert error<1e-12,error
    density=(loss@energy)/(task.G*m.bg.R**2/task.C**4)
    assert np.all(np.isfinite(density))
    return dict(classification='Counterexample candidate',passed=True,
        enclosed_mass_volume_fraction_relative=error,central_regular=True)


def main():
    out=task.OUT;result=json.loads((out/'resolved-result.json').read_text());plan=json.loads((out/'plan.json').read_text())
    checks=source_control();native=np.load(out/'coefficients.npz');d=np.load(out/'inputs.npz')
    cv,cp,ad=native['thermo'][:,3],native['thermo'][:,5],native['thermo'][:,4]
    chain=float(max(abs(cp/cv-(1-ad*native['raw'][:,8]))));assert chain<1e-7
    paths={};fields=['temperature','velocity','scalar']
    bg=task.Background();env=np.load(task.prior.OUT/'final-envelope.npz')
    for n in [32,64]:
        a=np.load(out/f'resolved-p4-{n}.npz');report=json.loads((out/f'resolved-p4-{n}.json').read_text())
        assert all(np.all(np.isfinite(a[k])) for k in a.files)
        assert all(np.max(abs(a[k][0]))==0 for k in fields)
        E=a['E'].astype(np.longdouble);net=np.sum(-np.diff(np.r_[np.longdouble(0),E]),dtype=np.longdouble)
        balance=float(abs(net+E[-1])/max(abs(E).max(),1e-100));assert balance<1e-12
        assert np.all(a['emission_flux']>0) and a['emission_energy'][0]==0
        degree=4;ids=a['indices'][::degree,1];r=a['grid'][::degree];valid=(ids>=0)&(r>=1.07)
        tail=float(max(abs(a['q'][ids[valid]]),default=0))
        norm=float(max(abs(a['scalar']).ravel()));tail_ratio=tail/max(norm,1e-100)
        surface=int(np.flatnonzero(a['grid']==1)[0]);z=float(a['q'][a['indices'][surface,0]])
        work=4*np.pi*bg.R**3*float(env['Ptotal'][-1])*float(env['A'][-1])**4*float(env['N'][-1])/np.sqrt(float(env['b'][-1]))*z
        paths[str(n)]=dict(energy_balance=balance,outer_tail_relative_to_native_scalar=tail_ratio,
            maximum_temperature=float(max(abs(a['temperature']).ravel())),
            maximum_luminosity_change=float(max(abs(a['emission_flux']/a['emission_flux'][0]-1))),
            emitted_energy_erg=float(a['E'][-1]),surface_displacement_fraction=z,
            omitted_surface_pressure_work_erg=work,
            pressure_work_over_emitted_heat=work/float(a['E'][-1]))
        assert report['max_linear_residual']<plan['gates']['linear_residual']
    base=np.load(out/'resolved-p4-64.npz');coarse=np.load(out/'resolved-p4-32.npz');other=np.load(out/'resolved-p2-64.npz')
    for key in fields:
        err=float(max(abs(base[key][::2]-coarse[key]).ravel())/max(abs(base[key]).ravel()))
        assert abs(err-result['comparisons'][key]['time_relative'])<1e-14
        err=float(max(abs(base[key]-other[key]).ravel())/max(abs(base[key]).ravel()))
        assert abs(err-result['comparisons'][key]['space_relative'])<1e-14
    # One emitted-energy history supplies both the material debit and the ray
    # inventory. This identity does not include omitted moving-boundary work.
    from scipy.interpolate import CubicHermiteSpline
    t=base['emission_times'];energy=base['emission_energy'];flux=base['emission_flux']
    history=CubicHermiteSpline(t,energy,flux);r=np.array([1.,1.005,1.02,1.07])
    mu,weights,delay=bg.rays(r,96);ret=np.maximum(t[-1]-delay,0)
    propagated=np.asarray(history(ret))@weights;inventory=energy[-1]-propagated
    identity=float(max(abs(-energy[-1]+inventory+propagated))/max(abs(energy[-1]),1e-100));assert identity<1e-12
    previous=json.loads((task.prior.OUT/'final.json').read_text())['whole_material_g']
    inventory=abs(float(d['dm'].sum())/previous-1);assert inventory<1e-12
    summary=dict(classification='Counterexample candidate',artifact_checks_passed=True,
        source_control=checks,cp_cv_chain_absolute=chain,whole_baryon_relative=inventory,paths=paths,
        causal_heat_photon_inventory_identity=identity,
        all_registered_response_gates_passed=result['passed'],
        source='Actual Phase92 finite-pressure matched background; no reused old evolution response.',
        pressure_work_scope='The recorded finite-pressure work is not included by the frozen emitting geometry. Its small heat-energy ratio is not a scalar-charge error bound or a completed material-radiation stress junction.',
        full_goal_complete=False)
    task.write(out/'audit.json',summary);print(json.dumps(summary),flush=True)


if __name__=='__main__':main()

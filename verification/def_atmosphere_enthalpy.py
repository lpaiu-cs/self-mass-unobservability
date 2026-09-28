"""First-law-preserving interpolation of accepted native atmospheric states."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.integrate import solve_ivp
import def_material_surface_enclosure as surface

h=surface.h
OUT=surface.atmosphere.parent.OUT/'enthalpy-coordinate'


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    prefix=surface.atmosphere.parent.OUT/'extended-cold/progress.json'
    source=json.loads(prefix.read_text());rows=source['rows'];raw=np.array([r['raw'] for r in rows])
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),prefix,surface.OUT/'result.json',surface.atmosphere.parent.BACKGROUND/'background-0.001.npz',surface.atmosphere.parent.BACKGROUND/'lapse.npz']},
        cause='Independent Pchip pressure/energy interpolation shifted the 946 K radius by 3.17 m despite tight ODE tolerance. Preserve it; impose dP=rho*dh along the interpolated adiabat.',
        method='Cubic Hermite P(eta), eta=(u+ground_offset+P/rho)/Ebase, node dP/deta=rho*Ebase. Set rho=P_eta/Ebase and u=h-P/rho; use only the accepted native prefix. No cold EOS extrapolation.',
        prefix='The extended native command ended with code 1 and no final result/failure record. Reuse only its 12 saved accepted roots, not a completed 100 K run.',
        gates=dict(node_pressure_relative=1e-12,node_density_relative=1e-12,positive_density=True,first_integral_absolute=1e-16,Just_radius_difference_m=.001),
        tolerances=[1e-10,1e-12],budget=dict(hard_timeout_seconds=30,native_calls=0,evolution_steps=0)))
    data,_=h.inputs();body=h.Structure(.001);bg=np.load(surface.atmosphere.parent.BACKGROUND/'background-0.001.npz');lap=np.load(surface.atmosphere.parent.BACKGROUND/'lapse.npz');y0=bg['faces'][0]
    R=y0[0]*body.R;M=y0[1]*body.B;mu=M/R;phib=.001*(1+body.mu*y0[3]);vb=.001*body.mu*y0[4]/body.R;nub=float(lap['nu_faces'][0]);cx=data['CX'][0];C=h.gr.C
    # Form the small thermal enthalpy before adding the large rest energy.
    excess=raw[:,2]+cx*raw[:,11]+raw[:,1]/raw[:,0];Ebase=excess[0];eta=excess/Ebase
    assert np.all(np.diff(eta)<0)
    h0=cx-cx*raw[0,11]*1e-4/C**2;dh=Ebase*1e-4/C**2
    curve=CubicHermiteSpline(eta[::-1],raw[::-1,1],raw[::-1,0]*Ebase,extrapolate=False);derivative=curve.derivative()
    temperature=PchipInterpolator(eta[::-1],np.array([r['lnT'] for r in rows])[::-1],extrapolate=False)
    nodeP=float(np.max(abs(curve(eta)/raw[:,1]-1)));nodeRho=float(np.max(abs(derivative(eta)/(Ebase*raw[:,0])-1)))
    dense=np.linspace(eta[-1],1,4001);assert np.all(curve(dense)>0) and np.all(derivative(dense)>0)
    # Each density polynomial is quadratic; include every internal extremum.
    for j in range(len(eta)-1):
        co=derivative.c[:,j];width=curve.x[j+1]-curve.x[j]
        t=-co[1]/(2*co[0]) if co[0] else -1
        if 0<t<width:assert (co[0]*t+co[1])*t+co[2]>0
    answers=[]
    def rhs(e,y):
        p=h.gr.G*float(curve(e))*.1/C**4;rho=h.gr.G*float(derivative(e))/Ebase*1000/C**2;H=h0+e*dh;en=rho*H-p
        r=R*(1+y[0]);m=R*(mu+y[1]);phi=phib+.001*y[2];v=vb+.001*y[3]/R;A=np.exp(-2*phi*phi);b=1-2*m/r;N=np.exp(nub+y[4]);pe=A**4*p;ee=A**4*en
        nr=m/(r*r*b)+4*np.pi*r*pe/b+r*v*v/2;dr=-dh/H/(nr-4*phi*v)
        vr=4*np.pi/b*(-4*phi*(ee-3*pe)+r*v*(ee-pe))-2*(r-m)/(r*r*b)*v
        return dr*np.array([1/R,(4*np.pi*r*r*ee+r*r*b*v*v/2)/R,v/.001,vr*R/.001,nr,4*np.pi*r*r*A**3*rho/(np.sqrt(b)*body.B)])
    for tol in [1e-10,1e-12]:
        sol=solve_ivp(rhs,(1,float(eta[-1])),np.zeros(6),method='DOP853',rtol=tol,atol=[tol*1e-5,1e-27,1e-22,1e-21,1e-22,1e-27],max_step=.002,dense_output=True)
        assert sol.success;end=sol.y[:,-1];grid=np.linspace(1,eta[-1],401);states=sol.sol(grid)
        phi=phib+.001*states[2];thermal=dh*(grid-1)/(h0+dh)
        defect=np.log1p(thermal)+states[4]-2*((phi-phib)*(phi+phib))
        answers.append(dict(tolerance=tol,radius_m=float(R*(1+end[0])),added_baryon_fraction=float(end[5]),max_first_integral_defect=float(abs(defect).max()),RHS_evaluations=sol.nfev))
        np.savez_compressed(OUT/f'atmosphere-{tol}.npz',eta=grid,state=states,pressure_cgs=curve(grid),rho_B_cgs=derivative(grid)/Ebase,lnT=temperature(grid),first_integral_defect=defect)
    mp=surface.mp;mp.mp.dps=75
    hb=mp.mpf(float(h0))+mp.mpf(float(dh));he=mp.mpf(float(h0))+mp.mpf(float(dh))*mp.mpf(float(eta[-1]))
    expected=surface.just_surface(mp,mp.mpf(float(mu)),mp.mpf(float(R*vb)),mp.mpf(float(phib)),mp.log(hb/he),mp.mpf(float(R)))[0]
    error=abs(answers[-1]['radius_m']-float(expected))
    record=dict(classification='Counterexample candidate',passed=bool(nodeP<1e-12 and nodeRho<1e-12 and max(a['max_first_integral_defect'] for a in answers)<1e-16 and error<.001),
        reused_native_points=len(rows),reused_extended_prefix_roots=source['new_roots'],last_temperature_K=float(np.exp(rows[-1]['lnT'])),last_pressure_dyn_cm2=float(raw[-1,1]),
        node_pressure_relative_error=nodeP,node_density_relative_error=nodeRho,positive_interpolated_density=True,answers=answers,
        exact_Just_radius_difference_m=error,radius_tolerance_comparison_m=abs(answers[-1]['radius_m']-answers[0]['radius_m']),seconds=time.monotonic()-start,
        native_calls=0,continuous_native_EOS_error_certified=False,physical_cold_tail_certified=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',record);print(json.dumps(record),flush=True);assert record['passed']


if __name__=='__main__':main()

"""Fixed support-interface mesh repair with the accepted full GR propagator."""
from pathlib import Path
import argparse
import json
import signal
import time
import resource
import numpy as np
import sympy as sp
import def_gr_full_weeks as go

PRIOR=go.OUT;OUT=PRIOR.parent/'def-gr-interface-patch';write=go.write
base=go.task.fem.base.task.coupled;OriginalBackground=base.Background
lu=go.splu;go.splu=lambda A:lu(A,permc_spec='NATURAL')
go.OUT=OUT;go.BETA=1024


def patch_cells():
    with np.load(PRIOR/'beta1024/p4-2048.npz') as d:
        mask=d['masks'][1];r=d['native_radius'][mask];cells=d['cells']
        ids=np.searchsorted(cells,r,side='right')-1
        # Freeze the complete original readout window, not a fitted selection
        # of points having the largest error. One neighbouring cell each side.
        return cells,np.arange(ids.min()-1,ids.max()+2)


class Background(OriginalBackground):
    def __init__(self,radiation,outer=2):
        super().__init__(radiation,outer)
        cells,ids=patch_cells();assert np.array_equal(self.grid[:self.surface_index+1],cells[cells<=1])
        extra=(cells[ids,None]+np.diff(cells)[ids,None]*np.arange(1,4)[None,:]/4).ravel()
        self.grid=np.sort(np.r_[self.grid,extra]);assert np.all(np.diff(self.grid)>0)
        self.surface_index=int(np.flatnonzero(self.grid==1)[0])
        self.nodes=self.sample(self.grid);self.mid=self.sample((self.grid[:-1]+self.grid[1:])/2)
        self.node_projection=radiation.projection(self.grid);self.mid_projection=radiation.projection(self.mid['r'])


def install():
    base.Background=Background


def control():
    # Splitting changes test-space resolution, never the input jump measure.
    x=sp.symbols('x');cuts=[sp.Rational(i,4) for i in range(5)]
    for degree in [1,2,4]:
        for n in range(degree+1):
            test=x**n
            integral=sum(sp.integrate(sp.diff(test,x),(x,a,b)) for a,b in zip(cuts[:-1],cuts[1:]))
            assert sp.simplify(integral-(test.subs(x,1)-test.subs(x,0)))==0
    return dict(classification='Proven',passed=True,
        identity='Subcell weak integration of constant g retains integral(test_prime*g)=g*(test(right)-test(left)); internal faces telescope. Original heat support and energy debit are unchanged.',
        scope='Weak constant-source identity only; actual variable-coefficient GR paths must pass unchanged gates.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90);start=time.monotonic()
    cells,ids=patch_cells()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Repair the remaining support-interface spatial wave resolution while retaining the jointly accepted full coupled time propagator.',
        evidence='Frozen p1/p2/p4 differences at the new boundary decrease12.0093% to6.09016%;99.999% of saved weighted squared difference lies near the last active source face. Saved acoustic crossing is1-3 original cells over the horizon. Direct lift subtraction is7e-13 of field difference. Existing source faces already match trial-cell boundaries.',
        reassessment='This is a separately planned local geometry candidate after the original fixed-grid branch stopped. Four subdivisions are chosen once: with the observed roughly first-order spatial reduction,6.09%/4 predicts1.52% below2% (a forecast, not a bound). Do not refine globally or iterate the split count after failure.',
        method='Split only cells containing the entire original new-interface readout mask plus one neighbour each side into4. Same source/background knots,heat poles,initial state,horizon,65 native readouts,degree1/2/4,quadrature6,beta1024,sigma12,contour4096 andcoefficients512/1024/2048. New mechanical trial nodes do not request EOS states or new heat cells.',
        patch=dict(first_cell=int(ids.min()),last_cell=int(ids.max()),original_cells=len(ids),new_cells=4*len(ids),
            radius_interval_R=[float(cells[ids.min()]),float(cells[ids.max()+1])],subdivisions=4),
        gates=dict(propagation_relative=.02,propagation_order=1.5,contour=.0002,spatial_relative=.02,spatial_decrease=True,
            coefficient=.02,outer=.002,quadrature=.002,abscissa=.0002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Pilot16 p4 resolvents. First p4 must forecast<450s and pass all four time gates. Fixed p2/p1 only after measured budget review. Stop on spatial failure without another patch. Four original coefficient/outer/quadrature/sigma contrasts remain conditional, with a measured total budget review before launch.',
        budget=dict(pilot_seconds=90,first_case_seconds=450,spatial_seconds=300,total_compute_seconds=1800,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(go.__file__),Path(go.space.__file__),
            PRIOR/'result.json',PRIOR/'spatial-location.json',PRIOR/'lift-check.json',PRIOR/'beta1024/p4-2048.npz',
            PRIOR.parent/'def-gr-spatial-repair/inspection.json',go.task.BANK/'fine-bank.npz',go.task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());install();setup_start=time.monotonic();p=go.Problem();setup=time.monotonic()-setup_start
    ids=np.linspace(1,go.COUNT//2,16,dtype=int);tick=time.monotonic();answers=np.array([p.transform(go.contour(k,12)[0]) for k in ids]);seconds=time.monotonic()-tick
    np.savez_compressed(OUT/'pilot.npz',ids=ids,answers=answers)
    old=json.loads((PRIOR/'pilot-budget.json').read_text());forecast=1.4*(setup+seconds/16*2048+old['inversion_forecast_seconds']+20)
    original=np.load(PRIOR/'beta1024/p4-2048.npz')
    for key,value in [('native_radius',p.model.original.native),('weights',p.model.original.weights),('masks',p.model.original.masks)]:
        assert np.array_equal(original[key],value),key
    flux,energy=p.model.heat.faces(1.)
    assert np.array_equal(original['heat_energy'],energy) and np.array_equal(original['heat_flux'],flux)
    assert np.all(np.isin(original['cells'],p.model.cells))
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,solve16_seconds=seconds,
        forecast_seconds=forecast,dofs=p.model.size,original_dofs=47417,linear_residual=p.error,seconds=time.monotonic()-start,
        unchanged_native_readouts=True,unchanged_heat_history=True,
        assumption='Actual patched p4 setup and16 complex solves scaled to2048 plus saved native inversion cost,20s output and40 percent margin. Smaller spaces and contrasts must be budgeted separately.'))
    signal.alarm(0);print('PATCH PILOT',forecast,'DOFS',p.model.size,flush=True)


def verify_plan():
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert go.task.digest(Path(p))==h,p


def run():
    verify_plan();assert not (OUT/'p4-result.json').exists()
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<450
    install();signal.alarm(450);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)))
    go.solve(reuse=True);signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()

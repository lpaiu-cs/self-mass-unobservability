"""Independent fixed-table scalar root at a failed conservative state."""
import json
import numpy as np
from scipy.optimize import brentq
import def_native_energy_temperature as task


def main():
    assert not (task.OUT/'independent-root.json').exists()
    d=np.load(task.OUT/'primitive-failure.npz');e=task.EOS(task.OUT/'runtime-columns.npz')
    errors=np.where(d['U'][0]>=e.floor,abs(d['error'])/np.maximum(abs(d['tau']),1e-100),0);i=int(np.argmax(errors));D,S,K=d['U'][:,i];tau=float(d['tau'][i]);lt=float(d['sigma'][i])
    def residual(t):
        p=0.
        for _ in range(4):
            v=S/(e.cx*D+tau+p);root=np.sqrt(1-v*v);rho=D*root;W=1/root;wm=v*v/(root*(1+root))
            pp,u,gamma,T,k=e(np.array([rho]),np.array([t]));p=float(pp[0])
        energy=e.cx*D*wm+(rho*u[0]+p)*W*W-p
        return float((energy-tau)/tau)
    lo,hi=e.limits(np.array([d['rho'][i]]));lower=float(lo[0])+1e-10;upper=float(hi[0])-1e-10
    root=brentq(residual,lower,upper,xtol=5e-15,rtol=1e-15)
    data=dict(classification='Counterexample candidate',cell=i,saved_logT=lt,saved_error=float(errors[i]),fixed_table_saved_error=residual(lt),independent_root_logT=root,
        independent_root_error=residual(root),bracket=[lower,upper],bracket_signs=[residual(lower),residual(upper)],
        nearby=[(float(z),residual(float(z))) for z in lt+np.array([-.001,-.0001,0,.0001,.001])])
    task.old.write(task.OUT/'independent-root.json',data);print(json.dumps(data),flush=True)


if __name__=='__main__':main()

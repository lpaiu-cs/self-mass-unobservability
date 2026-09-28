"""Expose, never accept, the preserved cold native failure; test saved precision."""
import ctypes
import json
import numpy as np
import mpmath as mp
import def_native_atmosphere_cold as cold
import gr_molecular_precision as precision

h=cold.h


def main():
    target=cold.OUT/'precision-probe.json';assert not target.exists()
    failed=json.loads((cold.OUT/'failure.json').read_text());data,_=h.inputs();X=data['X'][0]
    lr=failed['requested_logrho'];lt=failed['requested_lnT'];eos=h.molecular.model.EOS();captured=[];native=eos.call
    def capture(*args):
        native(*args);captured.append(args[-2].copy())
    eos.call=capture
    try:eos(2,lr,lt,X)
    except AssertionError:pass
    mp.mp.dps=65;extended=precision.EOS();ym=(X/precision.g.c.A)@extended.mapping;cx=float(ym@extended.weights);eps=ym/cx
    seed=np.zeros(24);seed[2]=1;seed/=seed@extended.weights
    rows=[]
    for dx in [0.,-1.,-2.]:
        extended.raw(0,mp.mpf(-20),mp.mpf(float(np.log(1e6))),seed,1.)
        try:
            a=extended.raw(2,mp.mpf(lr+dx)+mp.log(cx),mp.mpf(lt+dx/3),eps,cx)
            rows.append(dict(logrho=lr+dx,lnT=lt+dx/3,raw=[mp.nstr(v,40) for v in a],passed=True))
        except Exception as exc:rows.append(dict(logrho=lr+dx,lnT=lt+dx/3,error=repr(exc),passed=False))
    record=dict(classification='Counterexample candidate',failed_raw=[float(v) if np.isfinite(v) else str(v) for v in captured[-1]],nonfinite_columns=np.where(~np.isfinite(captured[-1]))[0].tolist(),extended=rows,
        physical_EOS_certified=False,scope='Same-input arithmetic diagnostic, no accepted entropy root or model substitution.')
    h.write(target,record);print(json.dumps(record),flush=True)


if __name__=='__main__':main()

import inspect,textwrap,numpy as np,json
import apply_retained_native_acoustic as a
a.initialize_photon();m=a.motion.State();m.select(64);m.folder=a.OUT/'collision-diagnostic';m.folder.mkdir(exist_ok=True)
s=textwrap.dedent(inspect.getsource(a.motion.State.point))
# The current public wrapper owns exact-zero inputs; use its enclosed nonzero function.
f=a.motion.State.point.__closure__
fn=next(c.cell_contents for c in f if callable(c.cell_contents))
s=textwrap.dedent(inspect.getsource(a.motion.State.point.__globals__['motion'].State.point)) if False else (a.OUT/'expanded-collision-source.py').read_text()
s=s.replace('if factor==1.:full=dict(t=t,photon=residual,bound=br,escape=er)',"if factor==1.:full=dict(t=t,photon=residual,bound=br,escape=er,cell_rounding=16*np.finfo(LD).eps*np.sum((abs(a['emit'])+abs(b['emit'])+(abs(a['loss'])+abs(b['loss']))*abs(I))*weight,axis=(1,2)),cell_source=np.sum(abs(residual)*weight,axis=(1,2)),nonzero=np.any(z!=0,axis=0),delta=z*AMP,coefficient_change=np.array([b[key]-a[key] for key in ['rho','beta','u','p']]))")
ns=dict(fn.__globals__);exec(compile(s,'acoustic153-small-source-diagnostic','exec'),ns);ns['point'](m,8,False)
z=np.load(m.folder/'point-8.npz');idx=np.argsort(z['cell_rounding'])[-10:][::-1]
print(json.dumps(dict(top=[dict(cell=int(i),rounding=float(z['cell_rounding'][i]),source=float(z['cell_source'][i]),nonzero=bool(z['nonzero'][i]),state_delta=z['delta'][:,i].tolist(),coefficient_change=z['coefficient_change'][:,i].tolist()) for i in idx],nonzero_deep=z['nonzero'][:19].tolist(),ratio_on_nonzero=float(np.sum(z['cell_rounding'][z['nonzero']])/np.sum(z['cell_source'])))))
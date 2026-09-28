"""Counterexample candidate: high-precision evaluation of native mass flux.

Reuse the primitive, reconstruction and HLL formulas. Keep the current stored
background/EOS jets, and promote arithmetic on the material direction. This is
an independent true native evaluation, not acceptance against a fitted matrix.
"""
from types import FunctionType
import inspect,textwrap
import numpy as np
import mpmath as mp


def scalar(v):
    if isinstance(v,mp.mpf):return v
    if hasattr(v,'as_integer_ratio'):
        p,q=v.as_integer_ratio();return mp.mpf(p)/q
    return mp.mpf(v)


def hp(v):
    a=np.asarray(v)
    if not a.ndim:return scalar(a[()])
    return np.frompyfunc(scalar,1,1)(a.astype(object))


def replace(source,changes):
    for old,new in changes:
        assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
    return source


def compile_function(fn,changes,extra):
    source=replace(textwrap.dedent(inspect.getsource(fn)),changes)
    ns=dict(fn.__globals__,**extra);exec(compile(source,__file__,'exec'),ns)
    return ns[fn.__name__]


def build(joint,owner):
    engine=joint.previous.engine;face=engine.face
    primitive=engine.tangent.__globals__['primitive']
    source=face.source.replace('np.zeros(m.n)','np.zeros(m.n,dtype=object)').replace(',float)',',object)')
    source=replace(source,[
        ("q=row['Q'];active=", "q=hp(row['Q']);active="),
        ("p,u,ut,uy,pt,py,*_=b.eos.gas(row['theta'],row['eta']);eta=(1+row['eta'])*dy[:nb]",
         "p,u,ut,uy,pt,py,*_=[hp(v) for v in b.eos.gas(row['theta'],row['eta'])];eta=hp(1+row['eta'])*dy[:nb]"),
        ('W=np.sqrt(W2)','W=hsqrt(W2)')])
    source=source.replace('m.a[','hp(m.a)[').replace('m.V[','hp(m.V)[').replace('model.cx','hp(model.cx)').replace('LD(C)','C')
    ns=dict(primitive.__globals__,LD=object,C=scalar(engine.C),hp=hp,hsqrt=np.frompyfunc(mp.sqrt,1,1))
    exec(compile(source,__file__,'exec'),ns);primitive_impl=ns['primitive']
    def primitive_high(m,k,z,field,bank):
        # The same bank also stores millions of photon coefficients, which
        # this material inverse never reads.
        keys=['beta','rho','p','u','dr_p','dt_p','dy_p','dr_u','dt_u','dy_u']
        return primitive_impl(m,k,hp(z),hp(field),{key:hp(bank[key]) for key in keys})
    reconstruction=compile_function(face.reconstruction,[
        ('ds=np.zeros_like(x)','ds=np.zeros_like(x,dtype=object)'),
        ('ds=np.zeros_like(y)','ds=np.zeros_like(y,dtype=object)')],{})
    # The geometric wrapper owns the actual conserved function; reuse its
    # original flat-metric part and then the unchanged lapse contribution.
    geometric=engine.flux_direction.__globals__['conserved']
    flat=geometric.__globals__['fixed_conserved']
    conserved=compile_function(flat,[
        ('a=np.asarray(f.base.af,LD);rho,v,lt,y=V','a=hp(f.base.af);rho,v,lt,y=hp(V)'),
        ('p,u,gamma=map(lambda x:np.asarray(x,LD),thermo[:3])','p,u,gamma=map(hp,thermo[:3])'),
        ("np.einsum('ajn,an->jn',jets(m,k,label,V,h),coords)","np.sum(hp(jets(m,k,label,V,h))*coords[:,None,:],axis=0)"),
        ('root=np.sqrt(1-v*v)','root=hsqrt(1-v*v)'),
        ('np.maximum(rho,f.eos.floor)','np.maximum(rho,hp(f.eos.floor))'),
        ('rho>f.eos.floor','rho>hp(f.eos.floor)'),
        ("np.maximum(H,LD('1e-100'))","np.maximum(H,hp(LD('1e-100')))"),
        ('cs=np.asarray(thermo[-1],LD)','cs=hp(thermo[-1])'),
        ('out=np.zeros_like(cs)','out=np.zeros_like(cs,dtype=object)'),
        ('return U[:3],F[:3],dU,dF,speed,dspeed','return hp(U[:3]),hp(F[:3]),dU,dF,speed,dspeed')],
        dict(hp=hp,hsqrt=np.frompyfunc(mp.sqrt,1,1)))
    geometric=FunctionType(geometric.__code__,dict(geometric.__globals__,fixed_conserved=conserved))
    flux=compile_function(engine.flux_direction,[
        ('zero=np.zeros(L.shape[1],dtype=LD)','zero=np.zeros(L.shape[1],dtype=object)'),
        ('*m.model.m.af*m.model.m.area','*hp(m.model.m.af)*hp(m.model.m.area)')],dict(conserved=geometric,hp=hp))
    source=replace(engine.source,[
        ("dV=np.array([V[0]*p['dr'][nb:],p['dv'][nb:],p['dt'][nb:],V[3]*p['dy'][nb:]],dtype=LD)",
         "dV=np.array([hp(V[0])*p['dr'][nb:],p['dv'][nb:],p['dt'][nb:],hp(V[3])*p['dy'][nb:]],dtype=object)"),
        ("djoin=np.array([join[0]*p['dr'][nb-1],p['dv'][nb-1],p['dt'][nb-1],join[3]*p['dy'][nb-1]],dtype=LD)",
         "djoin=np.array([hp(join[0])*p['dr'][nb-1],p['dv'][nb-1],p['dt'][nb-1],hp(join[3])*p['dy'][nb-1]],dtype=object)"),
        ('    owner=float(np.max(', '    factor=hp(factor)\n    owner=float(np.max('),
        ('F=np.c_[df,flux]','F=hp(np.c_[df,flux])'),
        ('flux[:3]+=(2*uf+ef)*baseline[:3]*factor[:3,None]',
         'flux[:3]+=hp((2*uf+ef)*baseline[:3])*factor[:3,None]')])
    ns=dict(engine.tangent.__globals__,primitive=primitive_high,reconstruction=reconstruction,flux_direction=flux,hp=hp)
    exec(compile(source,__file__,'exec'),ns)
    return ns['tangent']


def native_B(m,t,g,tangent,probe=1.,return_flux=False):
    material=m.material;material.full_field=m.geometry(t)[0]
    z=m.conserved(g);k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k])
    flux=np.zeros(m.n+1,dtype=object)
    for j,weight in [(k,1-w),(k+1,w)]:
        if weight:
            F,_,_=tangent(material,j,z,probe);flux+=scalar(weight)*hp(F[0])
    rate=-np.diff(flux)/hp(m.bu)
    return (rate,flux) if return_flux else rate


def cast(values):
    return np.frompyfunc(lambda v:np.longdouble(str(v)),1,1)(values).astype(np.longdouble)


def baryon_product(J,g):
    matrix=J[2::4].tocsr();values=hp(g.ravel());data=hp(matrix.data)
    return np.array([sum(data[a:b]*values[matrix.indices[a:b]],mp.mpf(0))
        for a,b in zip(matrix.indptr[:-1],matrix.indptr[1:])],object)

def initialize():
    Model=bind(base.initialize,OUT=OUT)();calls=[]
    template=prior.right.right_solver(calls);original=template.__globals__['gmres']
    inner=FunctionType(original.__code__,dict(original.__globals__,joint=SimpleNamespace(gmres=prior.scipy_gmres)),argdefs=original.__defaults__,closure=original.__closure__)
    source=inspect.getsource(joint.solve)
    for a,b in [('restart=20,maxiter=5','restart=80,maxiter=20'),('range(4)','range(12)'),('assert k<3,','assert k<11,'),
                ('if relative<1e-14 and max(moments)<1e-13:','if relative<1e-14 and max(moments)<1e-13 and gas_gate(gas_relative(m,residual,sol)):')]:
        assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(joint.solve.__globals__,gmres=inner,gas_relative=gas_relative,gas_gate=gas_gate);exec(compile(source,__file__,'exec'),ns)
    def solve(m,op,P,rhs,guess):
        start=len(calls)
        try:return ns['solve'](m,op,P,rhs,guess)
        except BaseException as exc:
            tb=exc.__traceback__
            while tb and not {'sol','residual','rhs'}.issubset(tb.tb_frame.f_locals):tb=tb.tb_next
            if tb:
                v=tb.tb_frame.f_locals;np.savez_compressed(OUT/'failed-linear.npz',rhs=rhs,guess=guess,solution=v['sol'],residual=v['residual'])
                write(OUT/'failed-linear.json',dict(error=repr(exc),gas_relative=gas_relative(m,v['residual'],v['sol'])))
            raise
        finally:write(OUT/'right-calls.json',dict(classification='Counterexample candidate',last_solve_begin=start,calls=calls))
    stage=Model.run.__globals__['stages'];source=(OUT/'expanded-full-stages.py').read_text();source=base.base.replace(base.base.replace(source,'    for newton in range(3):','    for newton in range(12):'),'if newton==2:','if newton==11:')
    old='if relative<1e-12 and max(moments)<1e-13:break'
    source=base.base.replace(source,old,'if relative<1e-12 and max(moments)<1e-13 and gas_gate(gas_relative(m,defect,sol)):break')
    source=base.base.replace(source,'audit.append(dict(relative=relative,moments=moments.astype(float).tolist()))','audit.append(dict(relative=relative,moments=moments.astype(float).tolist(),material_relative=gas_relative(m,defect,sol)))')
    ns_stage=dict(stage.__globals__,solve=solve,gas_relative=gas_relative,gas_gate=gas_gate);exec(compile(source,__file__,'exec'),ns_stage)
    Model.run=bind(Model.run,stages=ns_stage['stages']);return Model

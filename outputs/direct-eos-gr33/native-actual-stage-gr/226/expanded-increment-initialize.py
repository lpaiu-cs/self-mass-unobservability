def initialize():
    global Model
    prior.prior.run.prior.OUT=OUT;prior.prior.run.prior.initialize();Parent=owner.Model
    class Increment(Parent):
        def __init__(self,n):
            super().__init__(n);self.primary=self.driver;self.returned=prior.ReturnOnly(self.primary)
            self.select(self.returned);self.anchor=dict(np.load(prior.prior.saved(n)));self.selected={}
            self.anchor_checks=[];self.branch_checks=[];self.checked_times=set()
        def select(self,driver):self.driver=driver;self.redshift_driver=driver;self.material.driver=driver
        def anchor_state(self,t):
            k=int(np.argmin(abs(self.anchor['joint_stage_times']-t)))
            assert abs(self.anchor['joint_stage_times'][k]-t)<1e-18,('Not an actual saved stage',t)
            q=self.anchor['joint_stage_conserved_scaled'][k]
            return k,restored_gas(self,q)
        def native(self,t,g,probe=1.,details=False,tangent=None,metric=True):
            if tangent is not None:
                self.select(self.returned)
                return Parent.native(self,t,g,probe,details,tangent,metric)
            assert metric
            k,base=self.anchor_state(t);fn,reset,settle=paired_tangent();changes=[]
            try:
                for iteration in range(8):
                    if iteration:reset('high')
                    self.select(self.primary);high=Parent.native(self,t,base,probe,tangent=fn)
                    reset('low');self.select(self.returned);low=Parent.native(self,t,g,probe,True,fn)
                    changed=settle();changes.append(changed)
                    if not changed:break
                else:raise AssertionError(('Unsettled paired branches',changes))
                if t not in self.checked_times:
                    saved=self.anchor['joint_native_rates_scaled'][k]/self.units
                    error=(np.sum(abs(high-saved)*self.units,axis=0)/np.maximum(np.sum(abs(saved)*self.units,axis=0),LD('1e-290'))).astype(float)
                    self.anchor_checks.append(dict(time=float(t),native_relative=error.tolist()))
                    assert max(error)<1e-12,('Saved-stage anchor precision',error.tolist())
                    self.checked_times.add(t)
                self.selected[t]=(fn,lambda capture:reset('replay'))
                self.branch_checks.append(dict(time=float(t),iterations=len(changes),tie_switches=changes))
                return low if details else low[0]
            finally:self.select(self.returned)
    source=(OUT/'expanded-affine-jacobian.py').read_text()
    old="base=self.native(t,g);tangent,reset=selected_tangent()\n    selected=self.native(t,g,tangent=tangent)\n    assert np.array_equal(selected,base),'Selected branches must reproduce the original native direction'"
    assert source.count(old)==1
    source=source.replace(old,'base=self.native(t,g);tangent,reset=self.selected[t]')
    ns=dict(Parent.jacobian.__globals__);exec(compile(source,__file__,'exec'),ns);Increment.jacobian=ns['jacobian']
    (OUT/'expanded-increment-jacobian.py').write_text(source);Model=Increment

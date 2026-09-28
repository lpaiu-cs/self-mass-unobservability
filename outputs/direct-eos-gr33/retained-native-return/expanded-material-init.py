def __init__(self,reference,steps=128):
    base_init(self,reference);self.steps=steps;self.finite_resolution=0.;self.finite_calls=0;m=self.model
    common='momentum=self.area_gas*(self.face_p+ps+self.cx*self.face_rho*vf*vf)'
    centered='momentum=self.area_gas*(ps+self.cx*self.face_rho*vf*vf)'
    join='mass[-1],momentum[-1],energy[-1],neutral[-1]=self.mflux'
    subtract=join+'\n    momentum[-1]-=self.base_momentum_flux[-1]'
    # Change the common owner and its flux-export version together.
    physical=old.old.physical.branch.base.centered.method
    physical=replace(physical,common,centered);physical=replace(physical,join,subtract)
    physical=replace(physical,'np.diff(momentum-self.base_momentum_flux)','np.diff(momentum)')
    scope=dict(m.material_rhs.__func__.__globals__);exec(compile(physical,__file__,'exec'),scope)
    m.material_rhs=MethodType(scope['material_rhs'],m)
    deep=replace(self.expanded[1],common,centered);deep=replace(deep,join,subtract)
    deep=replace(deep,'np.diff(self.base_momentum_flux)+b.volume*self.f0["initial_support"]+geometry','b.volume*self.f0["initial_support"]+geometry')
    scope=dict(self.deep.__func__.__globals__);exec(compile(deep,__file__,'exec'),scope);self.deep=MethodType(scope['material_rhs'],m)
    self.centered_sources=(physical,deep)

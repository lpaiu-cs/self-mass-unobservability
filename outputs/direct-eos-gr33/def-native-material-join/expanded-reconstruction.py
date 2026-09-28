def reconstruct(self,V,t):
    m=self.base;delta=V[:3]-self.background_cell
    if self.boundary=='copy':ghost_delta=self.join_state[:3]-self.background_left[:,0]
    else:ghost_delta=np.array([0.,np.interp(t,m.hist_t,m.hist_v),0.])
    ext=np.column_stack([ghost_delta,delta,np.zeros(3)]);slope=np.zeros_like(ext)
    slope[:,1:-1]=thermo.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
    L=self.background_left+ext[:,:-1]+slope[:,:-1]/2;R=self.background_right+ext[:,1:]-slope[:,1:]/2
    if self.boundary=='copy':ghost=self.background_left[:,0]+ghost_delta
    else:ghost=np.array(m.left);ghost[1]=np.interp(t,m.hist_t,m.hist_v)
    ext=np.column_stack([ghost,V[:3],np.zeros(3)]);slope=np.zeros_like(ext)
    slope[:,1:-1]=thermo.task.minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
    left=ext[:,:-1]+slope[:,:-1]/2;right=ext[:,1:]-slope[:,1:]/2;left[:,:2]=L[:,:2];right[:,:2]=R[:,:2]
    assert min(left[0])>=-1e-13 and min(right[0])>=-1e-13
    left[0]=np.maximum(left[0],0);right[0]=np.maximum(right[0],0)
    y=V[3];ext=np.r_[self.incoming_y(y),y,self.eos.y0];slope=np.zeros_like(ext)
    slope[1:-1]=thermo.task.minmod(ext[1:-1]-ext[:-2],ext[2:]-ext[1:-1])
    return np.vstack([left,ext[:-1]+slope[:-1]/2]),np.vstack([right,ext[1:]-slope[1:]/2])

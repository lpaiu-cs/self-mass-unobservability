def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
    q=z['joint_stage_conserved_scaled'][2*step-1]
    initial=restored_gas(m,q);initial[~m.material.active(t)]=0
    v=m.pack(x,initial);src=[];rates=[];states=[];identity=[]
    for j,(xx,g,now,c,s) in enumerate(zip(photons,gas,times,cs,ss)):
        J,b=m.jacobian(now,g);affine=b-(J@g.ravel()).reshape(m.n,4)
        src.append(m.pack(s[0]/(m.scale*AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+affine))
        ph,q,*_=m.collision(c,xx,g,True);native=m.native(now,g)
        rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+ph+s[0]/(m.scale*AMP),q+native));states.append(m.pack(xx,g))
        old=z['joint_native_rates_scaled'][2*step+j];identity.append(float(np.max(abs(native*m.units-old))))
    rhs=(v+h*(A@np.array(src))).ravel();sol=np.array(states).ravel();defect=(np.array(states)-v-h*(A@np.array(rates))).ravel()
    rel=float(np.linalg.norm(defect)/np.linalg.norm(rhs));physical=(base.joint.physical_norm(m,defect)/base.joint.scales(m,rhs,sol)).astype(float).tolist()
    row=dict(step=step,relative=rel,physical=physical,native_identity_absolute=identity,passed=bool(rel<1e-12 and max(physical)<1e-13 and max(identity)==0))
    write(OUT/f'original-equation-{step}.json',row);assert row['passed'],row
    return row

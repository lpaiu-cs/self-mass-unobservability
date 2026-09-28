import signal
import def_photon_loss_bounds as m
original=m.FunctionType
def bounded(code,env,argdefs=None):
    fn=original(code,env,argdefs=argdefs)
    def invoke(row):
        signal.alarm(20)
        try:return fn(row)
        finally:signal.alarm(0)
    return invoke
m.FunctionType=bounded
m.fetch()

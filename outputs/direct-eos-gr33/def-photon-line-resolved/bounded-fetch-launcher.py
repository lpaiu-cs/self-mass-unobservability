import signal
import def_photon_line_resolved as m
original=m.previous.FunctionType
def bounded(code,env,argdefs=None):
    fn=original(code,env,argdefs=argdefs)
    def invoke(row):
        signal.alarm(20)
        try:return fn(row)
        finally:signal.alarm(0)
    return invoke
m.previous.FunctionType=bounded
m.fetch()

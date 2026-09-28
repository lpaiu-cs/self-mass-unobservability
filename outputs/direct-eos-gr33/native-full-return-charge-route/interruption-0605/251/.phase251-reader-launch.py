"""Serialize legacy shared-code initialization, retaining parallel field solves."""
from pathlib import Path
from types import FunctionType
import fcntl,runpy,sys
sys.argv=sys.argv[1:];sys.path.insert(0,'verification')
import read_complete_radau_history as base
original=base.endpoint.initialize
OUT=base.endpoint.OUT
def initialize():
    with Path('.native-reader-initialization.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        FunctionType(original.__code__,dict(original.__globals__,OUT=OUT),argdefs=original.__defaults__)()
base.endpoint.initialize=initialize
runpy.run_path(sys.argv[0],run_name='__main__')

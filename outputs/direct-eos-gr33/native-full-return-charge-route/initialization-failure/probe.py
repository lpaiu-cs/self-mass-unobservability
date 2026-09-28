from pathlib import Path
import sys
sys.path.insert(0,'verification');sys.argv=['probe','check_charge_field1288']
import read_full_return_charge as a
def trace(frame,event,arg):
 if event=='line' and frame.f_code.co_filename.endswith('solve_native_incident_lift.py') and frame.f_lineno==59:
  p=frame.f_globals['base'].drive.native.prior.OUT/'expanded-photon-run.py'
  print('GENERATED_PATH',str(p),flush=True)
 return trace
sys.settrace(trace)
try:a.bind(a.base.endpoint.initialize,OUT=a.charge.OUT/'field1288')()
finally:sys.settrace(None)
print('INITIALIZER_PASS')

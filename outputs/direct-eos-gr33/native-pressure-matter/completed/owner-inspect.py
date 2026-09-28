import inspect,json,numpy as np
import return_native_pressure_matter as run
run.initialize();m=run.c.Material(128,64)
for cls in type(m).__mro__:
 if 'rhs' in cls.__dict__:
  print(cls.__module__,cls.__name__,inspect.getsource(cls.__dict__['rhs']))
print('bank',run.c.feedback.old.OUT)
print('fields',max(np.max(abs(v)) for v in m.fields(m.t[1])[2:]))
print('raw',m.raw.__func__.__code__.co_filename)
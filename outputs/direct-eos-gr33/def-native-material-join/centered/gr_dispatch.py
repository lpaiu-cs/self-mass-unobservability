"""Reuse saved-source GR response on the corrected material trajectory."""
from types import FunctionType
import def_native_material_centered as run
import verify_native_material_join as prior

FunctionType(prior.gr.__code__,dict(vars(prior),run=run,OUT=run.OUT))()

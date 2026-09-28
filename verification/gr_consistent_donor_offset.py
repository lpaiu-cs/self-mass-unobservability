"""Counterexample candidate: retain the native composition anchor across pivots.

The branch extension must change both flux derivatives and the affine species
response. Recentring the species anchor on each artificial donor discards that
offset. Only the iteration operator changes; native acceptance stays fixed.
"""
from types import FunctionType

import gr_consistent_donor_tangent as prior

core, m, e = prior.core, prior.m, prior.e
OUT = m.BASE/'step42-consistent-offset'


class ConsistentOffset(core.PredictedDonorTangent):
    evaluate = prior.evaluate


context = dict(vars(core), OUT=OUT, __file__=__file__, PredictedDonorTangent=ConsistentOffset)
correction = FunctionType(core.correction.__code__, context)
context['correction'] = correction
finish = FunctionType(core.finish.__code__, context)
context['finish'] = finish
check = FunctionType(core.check.__code__, context)


if __name__ == '__main__':
    check()

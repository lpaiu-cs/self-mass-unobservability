"""Serialize finite mpmath intervals with outward decimal endpoints."""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction


def exact_endpoint(raw):
    sign,mantissa,exponent,bits=map(int,raw)
    assert bits>=0 and (mantissa or exponent==0),'Finite endpoints only'
    value=(-1 if sign else 1)*mantissa
    return Fraction(value*(1<<max(exponent,0)),1<<max(-exponent,0))


def interval_text(value,digits=60):
    endpoints=[]
    for raw,rounding in zip(value._mpi_,[ROUND_FLOOR,ROUND_CEILING]):
        number=exact_endpoint(raw)
        with localcontext() as ctx:
            ctx.prec=digits;ctx.rounding=rounding
            decimal=Decimal(number.numerator)/Decimal(number.denominator)
        represented=Fraction(decimal)
        assert represented<=number if rounding==ROUND_FLOOR else represented>=number
        endpoints.append(str(decimal))
    return '['+', '.join(endpoints)+']'


def control():
    from mpmath import iv
    iv.dps=40;failed_default=0;checked=0
    for denominator in range(3,80):
        for sign in [-1,1]:
            for scale in ['1e-35','1','1e35']:
                value=sign*iv.mpf(scale)/denominator
                original=list(map(exact_endpoint,value._mpi_))
                old=[Fraction(s.strip()) for s in str(value)[1:-1].split(',')]
                failed_default+=int(old[0]>original[0] or old[1]<original[1])
                new=[Fraction(s.strip()) for s in interval_text(value)[1:-1].split(',')]
                assert new[0]<=original[0] and new[1]>=original[1]
                checked+=1
    assert failed_default>0
    return dict(classification='Proven',passed=True,exact_endpoint_checks=checked,
        default_decimal_literal_inward_cases=failed_default,
        scope='Exact rational containment of serialized decimal endpoints. No claim that same-precision directed parsing of the former strings necessarily lost containment; no scientific tolerance was changed.')


if __name__=='__main__': print(control())

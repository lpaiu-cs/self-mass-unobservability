def build():
    fn=g.d.build_at
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save),argdefs=fn.__defaults__)(CACHE/'source',CACHE/'build',NAME,BRIDGE)
    lib=ctypes.CDLL(str(BRIDGE));a=np.zeros(5,dtype=np.int32)
    lib.precision_info.argtypes=[np.ctypeslib.ndpointer(np.int32,flags='C_CONTIGUOUS')];lib.precision_info(a)
    assert list(a)==[33,113,4931,128,53],a
    save('precision.json',dict(classification='Counterexample candidate',decimal_digits=int(a[0]),binary_significand_bits=int(a[1]),
        decimal_exponent_range=int(a[2]),storage_bits=int(a[3]),LAPACK_significand_bits=int(a[4])))

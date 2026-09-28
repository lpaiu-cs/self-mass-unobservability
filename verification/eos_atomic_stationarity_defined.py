"""Read PL weights only for the atomic elements whose branch defines them."""
import shutil, sys
import eos_atomic_stationarity as a

e=a.e;g=a.g
FAILED_OUT=e.OUT
e.OUT=g.OUT/'atomic-stationarity-defined';e.CACHE=g.CACHE/'atomic-stationarity-defined'
e.NAME='free_eos_direct24_atomic_defined'
e.LIB=e.CACHE/'build/src'/('lib'+e.NAME+'.so.1.0.0');a.OUT=e.OUT;OUT=e.OUT


def prepare():
    a.prepare()
    path=e.CACHE/'source/src/ionize.f90';text=path.read_text()
    old='  if(ifpl.eq.1) station_plop=plop';assert text.count(old)==1;text=text.replace(old,'')
    old='        station_active(ielement)=1;station_zero(ielement)=dvzero(ielement)'
    assert text.count(old)==1
    text=text.replace(old,old+'''
        if(ifpl.eq.1) station_plop(ion0+1:ion0+iatomic_number(ielement))= &
             plop(ion0+1:ion0+iatomic_number(ielement))''')
    path.write_text(text);shutil.copy2(path,OUT/'ionize.f90')
    plan=e.read(OUT/'plan.json')
    plan['source_changes']['ionize.f90']['after']=g.c.sha(path)
    plan['bindings'].update({rel:g.c.sha(g.ROOT/rel) for rel in [
        'verification/eos_atomic_stationarity_defined.py','outputs/direct-eos-gr33/atomic-stationarity/failure-manifest.json']})
    plan['read_contract_revision']=dict(
        reason='The original diagnostic encountered a nonfinite PL slot 261 after states 0,1,2. It belongs to an inactive element. A fresh later repeat had no nonfinite slot, consistent with undefined storage rather than a stable physical value. Only the source-defined active atomic PL slices are read; all other trace entries are padding and excluded from comparisons.',
        excluded_molecular_PL_slots=True,original_failure_preserved=True,physics_and_tolerances_unchanged=True)
    e.save('plan.json',plan)


build=a.build;control=a.control;run=a.run
if __name__=='__main__': globals()[sys.argv[1]]()

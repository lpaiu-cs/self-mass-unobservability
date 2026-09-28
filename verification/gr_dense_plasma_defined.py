"""Literal-aware return audit and derivatives of the unchanged screening fit.

Keep the failed first audits and their frozen scripts. Reuse their workflows
through explicit, reversible source substitutions in separately bound files.
"""
from types import FunctionType
import ctypes,inspect,json,shutil,sys
import gr_dense_plasma as d
import gr_dense_plasma_audit as audit
import gr_dense_plasma_screening_repair as screening

g=d.g;OUT=g.OUT/'gr-dense-plasma-defined';CACHE=g.CACHE/'dense-plasma-defined'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    assert not json.loads((screening.OUT/'controls.json').read_text())['passed']
    assert not json.loads((audit.OUT/'return-audit.json').read_text())['source_identity_passed']
    replacements={
        'return-run.py':(inspect.getsource(audit.run),[('-1.323515','-float(np.float32(1.323515))')]),
        'screen-free.py':(inspect.getsource(screening.free_mp),[("q('.001')","mp.mpf(float('.001'))")]),
        'screen-controls.py':(inspect.getsource(screening.controls),[('before[:6],after[:6]','before[[0,1,3,4]],after[[0,1,3,4]]')]),
        'screen-run.py':(inspect.getsource(screening.run),[('b[:6],saved[\'corrected\'][k,:6]',"b[[0,1,3,4]],saved['corrected'][k,[0,1,3,4]]")])}
    changes={}
    for name,(before,pairs) in replacements.items():
        after=before
        for old,new in pairs:assert after.count(old)==1;after=after.replace(old,new)
        reverse=after
        for old,new in reversed(pairs):reverse=reverse.replace(new,old)
        assert reverse==before;(OUT/name).write_text(after);changes[name]=pairs
    for branch,module in [('return',audit),('screening',screening)]:
        folder=OUT/branch;folder.mkdir();plan=json.loads((module.OUT/'plan.json').read_text())
        plan['literal_correction']='Keep the first failed audit. The ideal-ion constant 1.323515 is a default-real literal. TY1 uses double literal 1.d-3. The screening free energy H1 numerator is 1+X^2/5, whose derivative coefficient is 2/5, not the independently binary32-rounded literal .4.'
        if branch=='screening':
            plan['all_unaffected_outputs_bitwise']='F,U,S,CV only; P,PDT,PDR are now derived consistently and may change.'
            before=(screening.OUT/'eos22-screening.f').read_text()
            pairs=[('H1X=.4*X/H1U','H1X=(2.d0/5.d0)*X/H1U'),
                ('H1*(.4/H1U-(.4*X/H1U)**2','H1*((2.d0/5.d0)/H1U-((2.d0/5.d0)*X/H1U)**2')]
            after=before
            for old,new in pairs:assert after.count(old)==1;after=after.replace(old,new)
            reverse=after
            for old,new in reversed(pairs):reverse=reverse.replace(new,old)
            assert reverse==before;(folder/'eos22-screening.f').write_text(after);changes['screening/eos22-screening.f']=pairs
        else:shutil.copy2(audit.OUT/'calibration.json',folder/'calibration.json')
        (folder/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
    paths=[g.ROOT/'verification/gr_dense_plasma_defined.py',screening.OUT/'plan.json',screening.OUT/'controls.json',
        audit.OUT/'plan.json',audit.OUT/'return-audit.json']+[p for p in OUT.rglob('*') if p.is_file()]
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='e9500e5',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},reversible_substitutions=changes,
        policy='Retain both original audit failures. Same 1e-10 independent and group-identity gates, same 1e-5 finite mixture-derivative gate. Classical screening free energy is unchanged; no physical coefficients are fitted. Only derivative algebra and the independent literal reader change.',
        numerical_not_physical_certification=True))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel


class Provider(d.Provider):
    def __init__(self):self.lib=ctypes.CDLL(str(CACHE/'pc.so'))


def namespace(module,branch):
    folder=OUT/branch
    def save_branch(name,value):(folder/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    def bound():
        bindings()
        return FunctionType(module.bindings.__code__,dict(module.bindings.__globals__,OUT=folder))()
    ns=dict(module.run.__globals__,OUT=folder,CACHE=CACHE,save=save_branch,bindings=bound)
    def verify_branch():FunctionType(module.verify.__code__,ns)()
    ns['verify']=verify_branch
    if branch=='screening':
        ns['Provider']=Provider
        exec(compile((OUT/'screen-free.py').read_text(),str(OUT/'screen-free.py'),'exec'),ns)
        ns['independent']=FunctionType(module.independent.__code__,ns)
    return ns


def return_audit():
    ns=namespace(audit,'return');exec(compile((OUT/'return-run.py').read_text(),str(OUT/'return-run.py'),'exec'),ns);ns['run']()


def build():
    ns=namespace(screening,'screening');FunctionType(screening.build.__code__,ns)()


def controls():
    ns=namespace(screening,'screening');exec(compile((OUT/'screen-controls.py').read_text(),str(OUT/'screen-controls.py'),'exec'),ns);ns['controls']()


def run():
    ns=namespace(screening,'screening');exec(compile((OUT/'screen-run.py').read_text(),str(OUT/'screen-run.py'),'exec'),ns);ns['run']()
    namespace(audit,'return')['verify']()
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.rglob('*') if p.is_file()}));verify()


def verify():
    bindings()
    for module,branch in [(audit,'return'),(screening,'screening')]:namespace(module,branch)['verify']()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS literal-aware dense plasma contract and derivative repairs; original failures preserved',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()

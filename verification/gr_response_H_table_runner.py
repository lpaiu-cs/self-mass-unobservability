"""Preserve the Python-version failure; replay hashing with the portable API."""
from types import ModuleType
import json,sys
import gr_response_H_table as original

OUT=original.OUT.parent/'gr-response-H-table-defined'


def source():
    before=(original.ROOT/'verification/gr_response_H_table.py').read_text()
    changes={"OUT=product.OUT.parent/'gr-response-H-table'":"OUT=product.OUT.parent/'gr-response-H-table-defined'",
        "hashlib.file_digest(source,'sha256').hexdigest()":"hashlib.sha256(source.read()).hexdigest()"}
    after=before
    for old,new in changes.items():assert after.count(old)==1;after=after.replace(old,new)
    return after,changes


def module(text):
    name='gr_response_H_table_defined';obj=ModuleType(name);sys.modules[name]=obj;exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    assert not OUT.exists();text,changes=source();obj=module(text);obj.prepare();(OUT/'candidate.py').write_text(text)
    error=original.ROOT/'outputs/gr-response-H-table33-preflight.log';assert "module 'hashlib' has no attribute 'file_digest'" in error.read_text()
    p=json.loads((OUT/'plan.json').read_text());p.update(checkpoint='2028df24',substitutions=changes,
        correction='The WSL Python lacks hashlib.file_digest. The first native table and its lossless gzip were written, then the replay hash check raised AttributeError. Use the existing hashlib.sha256 on the decompressed bytes instead. Each table is only several MB. Native equations, all inputs, budgets, state order, coefficient checks and independent controls are unchanged.')
    files=[original.ROOT/'verification/gr_response_H_table_runner.py',OUT/'candidate.py',original.OUT/'plan.json',error]
    p['bindings'].update({x.relative_to(original.ROOT).as_posix():original.g.sha(x) for x in files});obj.save('plan.json',p)


def run(name):
    text,_=source();assert text==(OUT/'candidate.py').read_text();obj=module(text);obj.bindings();getattr(obj,name)()


def verify_preflight():run('verify_preflight')
def verify():run('verify')


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:run(sys.argv[1])

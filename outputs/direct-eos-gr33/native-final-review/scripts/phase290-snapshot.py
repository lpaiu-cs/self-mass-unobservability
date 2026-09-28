"""Phase 290: build a public snapshot commit of local main on top of origin/main without touching local history.

Included: every blob of main's tree except *.npz, *.npy, *.gz and files larger than 50 MB. A PUBLIC_SNAPSHOT_OMITTED.tsv (path, size,
blob id) is added. Included text blobs up to 5 MB are scanned for common credential patterns; any hit aborts. Prints the commit id; the
push is a separate step (git push origin <commit>:refs/heads/main, a fast-forward of origin/main).
Usage: python phase290-snapshot.py <repo dir>
"""
import os, re, subprocess, sys, tempfile
repo = sys.argv[1]
def git(*args, inp=None, env=None):
    return subprocess.run(['git', '-c', 'gc.auto=0', '-c', 'maintenance.auto=false', *args], cwd=repo, input=inp, capture_output=True, check=True, env=env).stdout
main = git('rev-parse', 'main').decode().strip(); base = git('rev-parse', 'origin/main').decode().strip()
entries = []
for rec in git('ls-tree', '-r', '-l', '-z', 'main').split(b'\0'):
    if not rec: continue
    meta, path = rec.split(b'\t', 1); mode, typ, sha, size = meta.split()
    entries.append((mode.decode(), typ.decode(), sha.decode(), int(size) if size != b'-' else 0, path.decode('utf-8', 'surrogateescape')))
assert all(t in ('blob', 'commit') for _, t, _, _, _ in entries), 'unexpected tree entry type'
print('gitlinks kept as links:', [e[4] for e in entries if e[1] == 'commit'])
EXCL = ('.npz', '.npy', '.gz')
keep = [e for e in entries if e[1] == 'commit' or (not e[4].lower().endswith(EXCL) and e[3] <= 50_000_000)]
omit = [e for e in entries if e not in keep]
print('main', main[:12], 'base', base[:12], 'files', len(entries), 'keep', len(keep), f'{sum(e[3] for e in keep)/1e9:.3f} GB', 'omit', len(omit), f'{sum(e[3] for e in omit)/1e9:.2f} GB')
pat = re.compile(rb'(ghp_[A-Za-z0-9]{36}|github_pat_[A-Za-z0-9_]{40,}|sk-[A-Za-z0-9_-]{32,}|AKIA[0-9A-Z]{16}|-----BEGIN [A-Z ]*PRIVATE KEY-----|xox[baprs]-[0-9A-Za-z-]{10,})')
scan = [e for e in keep if e[1] == 'blob' and e[3] <= 5_000_000]
proc = subprocess.Popen(['git', 'cat-file', '--batch'], cwd=repo, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
hits = []
for mode, typ, sha, size, path in scan:
    proc.stdin.write(sha.encode() + b'\n'); proc.stdin.flush()
    header = proc.stdout.readline(); n = int(header.split()[2]); data = proc.stdout.read(n); proc.stdout.read(1)
    if b'\0' in data[:8000]: continue
    for m in pat.finditer(data): hits.append((path, m.group(0)[:12]))
proc.stdin.close(); proc.wait()
print('scanned', len(scan), 'credential-pattern hits', len(hits), hits[:10])
if hits: sys.exit('aborted: credential-like strings found')
tsv = 'path\tsize_bytes\tgit_blob\n' + ''.join(f'{p}\t{s}\t{h}\n' for _, _, h, s, p in sorted(omit, key=lambda e: e[4]))
tsv_sha = git('hash-object', '-w', '--stdin', inp=tsv.encode('utf-8')).decode().strip()
with tempfile.TemporaryDirectory() as d:
    env = dict(os.environ, GIT_INDEX_FILE=os.path.join(d, 'index'))
    info = ''.join(f'{m} {h}\t{p}\n' for m, _, h, _, p in keep) + f'100644 {tsv_sha}\tPUBLIC_SNAPSHOT_OMITTED.tsv\n'
    git('update-index', '--add', '--index-info', inp=info.encode('utf-8', 'surrogateescape'), env=env)
    tree = git('write-tree', env=env).decode().strip()
msg = (f'Public snapshot of {main} (2026-09-28)\n\nWorking-tree snapshot without binary arrays (*.npz, *.npy, *.gz) and files larger than 50 MB; '
       f'{len(omit)} omitted files are listed in PUBLIC_SNAPSHOT_OMITTED.tsv and bound by SHA-256 in the manifests. See PUBLIC_SNAPSHOT.md.\n')
commit = git('commit-tree', tree, '-p', base, '-m', msg).decode().strip()
print('snapshot commit', commit, 'tree', tree)

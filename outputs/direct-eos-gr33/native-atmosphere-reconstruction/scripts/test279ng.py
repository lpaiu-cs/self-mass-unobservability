import sys, time
src = open('phase279-nongray.py', encoding='utf-8').read().split('# self-checks')[0]
exec(src)
G = atmosphere(15898.27, 5.73962, gray=True); t = time.monotonic()
N = atmosphere(15898.27, 5.73962, gray=False, T_init=G['T'], verbose=True)
print('iterations', N['iterations'], N['final'], 'T0/Teff', N['T'][0]/15898.27, '%.0fs' % (time.monotonic() - t))

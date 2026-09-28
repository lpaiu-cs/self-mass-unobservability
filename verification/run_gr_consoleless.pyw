"""Run the frozen GR continuation without a Windows terminal console."""
import ctypes
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-resume271-consoleless'


def write(name, value):
    with (OUT/name).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2)


def main():
    assert os.name == 'nt'
    console = ctypes.windll.kernel32.GetConsoleWindow()
    assert console == 0, 'Use pythonw.exe; the launcher must have no console.'
    write('windows-launcher.json', dict(pid=os.getpid(),
        started_utc=datetime.now(timezone.utc).isoformat(),
        task_name='SelfMassGR-Consoleless-20260917', executable=sys.executable,
        console_hwnd=console, purpose='Console-free launcher and CREATE_NO_WINDOW WSL child.'))
    command = [str(Path(os.environ['SystemRoot'])/'System32/wsl.exe'),
        '-d', 'Ubuntu-22.04', '--cd', '/mnt/e/lab/self-mass-unobservability', '--',
        '/usr/bin/env', 'OPENBLAS_NUM_THREADS=1',
        'PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification',
        '/usr/bin/taskset', '-c', '0-15', '/usr/bin/python3', '-u',
        'verification/gr_resume271_consoleless.py', 'chain', '--workers', '15']
    with (OUT/'run.log').open('xb') as log:
        child = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL,
            stdout=log, stderr=subprocess.STDOUT, creationflags=subprocess.CREATE_NO_WINDOW)
        write('windows-wsl-child.json', dict(pid=child.pid, parent_pid=os.getpid(),
            creationflags=subprocess.CREATE_NO_WINDOW, command=command))
        code = child.wait()
    write('windows-runner-exit.json', dict(exit_code=code,
        ended_utc=datetime.now(timezone.utc).isoformat()))
    return code


if __name__ == '__main__':
    try:
        exit_code = main()
    except Exception:
        write('windows-runner-error.json', dict(traceback=traceback.format_exc(),
            ended_utc=datetime.now(timezone.utc).isoformat()))
        raise
    sys.exit(exit_code)

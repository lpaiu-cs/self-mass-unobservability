$ErrorActionPreference = 'Stop'
$grRoot = 'E:\lab\self-mass-unobservability'
$grOutput = Join-Path $grRoot 'outputs\direct-eos-gr33\gr-compatible-equilibrium-v2\production-resume271'
$grMarker = Join-Path $grOutput 'windows-launcher.json'
$grEncoding = [System.Text.UTF8Encoding]::new($false)
$grRecord = [ordered]@{
    pid = $PID
    started_utc = [DateTime]::UtcNow.ToString('o')
    task_name = 'SelfMassGR-Resume271-20260917'
    purpose = 'Demand-only Windows launcher holding the foreground WSL computation outside the harness.'
}
$grStream = [System.IO.File]::Open($grMarker, [System.IO.FileMode]::CreateNew, [System.IO.FileAccess]::Write)
try {
    $grBytes = $grEncoding.GetBytes(($grRecord | ConvertTo-Json))
    $grStream.Write($grBytes, 0, $grBytes.Length)
} finally { $grStream.Dispose() }
$grCommand = 'set -C; exec /usr/bin/env OPENBLAS_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification /usr/bin/taskset -c 0-15 /usr/bin/python3 -u verification/gr_resume271.py chain --workers 15 > outputs/direct-eos-gr33/gr-compatible-equilibrium-v2/production-resume271/run.log 2>&1'
& "$env:SystemRoot\System32\wsl.exe" -d Ubuntu-22.04 --cd /mnt/e/lab/self-mass-unobservability -- /bin/bash -c $grCommand
$grExit = $LASTEXITCODE
$grFinished = @{ exit_code = $grExit; ended_utc = [DateTime]::UtcNow.ToString('o') } | ConvertTo-Json
[System.IO.File]::WriteAllText((Join-Path $grOutput 'windows-runner-exit.json'), $grFinished, $grEncoding)
exit $grExit

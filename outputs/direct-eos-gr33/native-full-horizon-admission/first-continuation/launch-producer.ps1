$ErrorActionPreference = 'Stop'
$researchRuntime = '\\wsl.localhost\Ubuntu-22.04\home\lpaiu\work\native-retained-tail-runtime'
$researchWork = Join-Path $researchRuntime 'native-full-horizon185-work'
$launchReceipt = Join-Path $researchWork 'launch.json'
if (Test-Path -LiteralPath $launchReceipt) { throw 'Existing launch receipt: inspect the exact process before any restart.' }
$status = Get-Content -LiteralPath (Join-Path $researchWork 'status.json') -Raw | ConvertFrom-Json
if ($status.state -ne 'prepared') { throw 'Producer is not in the prepared state.' }
$researchLogs = Join-Path $PSScriptRoot 'outputs\direct-eos-gr33\native-full-horizon-live'
New-Item -ItemType Directory -Path $researchLogs -Force | Out-Null
$stdout = Join-Path $researchLogs 'controller.stdout.log'
$stderr = Join-Path $researchLogs 'controller.stderr.log'
if ((Test-Path -LiteralPath $stdout) -or (Test-Path -LiteralPath $stderr)) { throw 'Existing controller logs must be inspected, not overwritten.' }
$arguments = @('-d','Ubuntu-22.04','--cd','/home/lpaiu/work/native-retained-tail-runtime','--exec','env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification','/usr/bin/python3','verification/complete_full_incident_horizon.py','run')
$process = Start-Process -FilePath (Get-Command wsl.exe).Source -ArgumentList $arguments -WindowStyle Hidden -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru
$receipt = [ordered]@{
    windows_pid = $process.Id
    started_utc = [DateTime]::UtcNow.ToString('o')
    arguments = $arguments
    stdout = $stdout
    stderr = $stderr
    source_sha256 = (Get-FileHash -LiteralPath (Join-Path $researchRuntime 'verification\complete_full_incident_horizon.py') -Algorithm SHA256).Hash.ToLowerInvariant()
    plan_sha256 = (Get-FileHash -LiteralPath (Join-Path $researchWork 'plan.json') -Algorithm SHA256).Hash.ToLowerInvariant()
}
$receipt | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $launchReceipt -Encoding utf8NoBOM
Copy-Item -LiteralPath $PSCommandPath -Destination (Join-Path $researchWork 'launch-producer.ps1')
$receipt | ConvertTo-Json -Depth 5

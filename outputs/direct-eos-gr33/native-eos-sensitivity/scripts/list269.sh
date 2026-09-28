# List the phase-269 artifacts and sizes (read-only).
cd /home/lpaiu/work/native-retained-tail-runtime && ls -la .phase269-* | awk '{print $5, $9}'
echo == libs; ls -la /home/lpaiu/work/direct-eos-gr33/native-cold-population-*/phase269-build.json /home/lpaiu/work/direct-eos-gr33/photon-eos-levels-*/phase269-build.json | awk '{print $5, $9}'
echo == ploff runtime; cd /home/lpaiu/work/native-eos269-runtime && ls -la .phase269-* | awk '{print $5, $9}'
ls primary269-eos64-work | grep -E 'json$' | tr '\n' ' '; echo
ls .phase269-eos64-logs | tr '\n' ' '; echo
ls readout269-eos61-work | tr '\n' ' '; echo
echo == attempt1; cd /home/lpaiu/work/native-eos269-runtime-attempt1 && ls -la .phase269-* | awk '{print $5, $9}'

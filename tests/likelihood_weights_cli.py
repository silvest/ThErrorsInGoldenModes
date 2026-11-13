"""Real GoldenModes driver: per-term ROOT output in PT and ordinary MPI runs."""
import argparse
from pathlib import Path
import subprocess
import tempfile

p = argparse.ArgumentParser()
p.add_argument('--executable', required=True)
p.add_argument('--checker', required=True)
p.add_argument('--mpiexec', required=True)
p.add_argument('--mpi-flags', default='')
a = p.parse_args()
exe, checker = str(Path(a.executable).resolve()), str(Path(a.checker).resolve())
help_text = subprocess.run([exe, '--help'], capture_output=True, text=True, check=True).stdout
assert '--write-mcmc-weights' in help_text
flags = [f for f in a.mpi_flags.split(';') if f]
with tempfile.TemporaryDirectory(prefix='gm-weights-') as temporary:
    work = Path(temporary)
    for pt in (False, True):
        prefix = 'pt_' if pt else 'ordinary_'
        command = [a.mpiexec, '-n', '3', *flags, exe,
                   '--indSU3', '--flagBJPSIP', '--minimal', '--Univariate',
                   '--write-mcmc-weights', '--nChains', '2', '--nPreRun', '8', '--nRun', '4',
                   '--nUpdateMax', '4', '--prerun-policy', 'staged', '--settling-sweeps', '3',
                   '--outfile', prefix]
        if pt:
            command += ['--pt-temperatures', '1,2,4', '--pt-full-priors', '--pt-swap-interval', '2']
        result = subprocess.run(command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
        assert result.returncode == 0, result.stdout[-15000:]
        path = work / (prefix + ('cold_chains.root' if pt else 'chains.root'))
        checked = subprocess.run([checker, 'check', str(path), '8'], cwd=work, capture_output=True, text=True, timeout=15)
        assert checked.returncode == 0, checked.stdout + checked.stderr
print('PASS: GoldenModes flag enables aligned likelihood-weight output in ordinary and PT runs')

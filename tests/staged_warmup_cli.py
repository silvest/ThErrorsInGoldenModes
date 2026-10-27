"""CLI validation and a short GoldenModes PT run with discarded settling."""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile
p=argparse.ArgumentParser()
p.add_argument('--executable',required=True)
p.add_argument('--mpiexec',required=True)
p.add_argument('--mpi-flags',default='')
p.add_argument('--smooth',action='store_true')
a=p.parse_args();exe=str(Path(a.executable).resolve())
help_result=subprocess.run([exe,'--help'],text=True,capture_output=True,check=True)
assert '--prerun-policy' in help_result.stdout and '--settling-sweeps' in help_result.stdout and '--smooth-adaptation' in help_result.stdout
for flags in [ ['--prerun-policy','unknown'], ['--prerun-policy'],
               ['--settling-sweeps','-1'], ['--settling-sweeps','2x'], ['--settling-sweeps','1'],
               ['--prerun-policy','staged','--load-proposal','unused'],
               ['--smooth-adaptation','--load-proposal','unused'],
               ['--smooth-adaptation','--nPreRun','12','--nPreRunFactorized','12'],
               ['--prerun-policy','staged','--nPreRun','12','--nPreRunFactorized','12'],
               ['--prerun-policy','staged','--Univariate','--nPreRun','2147483647','--settling-sweeps','1'] ]:
    result=subprocess.run([exe,*flags],text=True,capture_output=True,timeout=10)
    assert result.returncode!=0,(flags,result.stdout,result.stderr)
with tempfile.TemporaryDirectory(prefix='gm-staged-') as directory:
    flags=[f for f in a.mpi_flags.split(';') if f]
    command=[a.mpiexec,'-n','3',*flags,exe,'--indSU3','--flagBJPSIP','--minimal','--Univariate',
             '--pt-full-priors','--pt-temperatures','1,4','--nChains','2','--pt-swap-interval','2',
             '--prerun-policy','staged','--nPreRun','12','--settling-sweeps','7','--nRun','9',
             '--nUpdateMax','5','--pt-diagnostics-interval','4','--pt-trace-interval','2','--outfile','run_']
    if a.smooth:command.append('--smooth-adaptation')
    result=subprocess.run(command,cwd=directory,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=100)
    assert result.returncode==0,result.stdout[-15000:]
    base=Path(directory)
    with (base/'run_warmup_status.csv').open() as f:status=list(csv.DictReader(f))[0]
    assert status['policy']=='staged' and status['stage']=='complete'
    assert status['scale_adaptation']==('smooth' if a.smooth else 'legacy')
    assert float(status['scale_target'])==.25
    assert int(status['factorized_scale_updates'])==(3 if a.smooth else 0)
    assert status['adaptation_sweeps']=='12' and status['settling_sweeps']=='7'
    assert status['stop_reason']=='adaptation_budget_exhausted' and status['posterior_convergence_assessed']=='0'
    with (base/'run_pt_checkpoints.csv').open() as f:ledger=list(csv.DictReader(f))
    assert {x['phase'] for x in ledger}=={'prerun','warmup','production'}
    for phase,n in [('prerun',12),('warmup',7),('production',9)]:
        phase_rows=[x for x in ledger if x['phase']==phase]
        assert phase_rows[-1]['last_sweep']==str(n),(phase,phase_rows)
    assert (base/'run_proposal.txt').is_file() and (base/'run_cold_chains.root').stat().st_size>1000
print('PASS: staged CLI validation, real GoldenModes MPI run, independent phase checkpoints/status/cache')

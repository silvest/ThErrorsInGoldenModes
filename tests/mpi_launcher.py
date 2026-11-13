"""Dry-run PBS launchers; verify real MPI startup rejects a wrong allocation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

p=argparse.ArgumentParser()
p.add_argument('--script',required=True)
p.add_argument('--executable',required=True)
p.add_argument('--mpiexec',required=True)
p.add_argument('--mpi-flags',default='')
a=p.parse_args()
with tempfile.TemporaryDirectory(prefix='gm-launcher-') as directory:
    base=Path(directory);work=base/'working directory';work.mkdir()
    fake=base/'fake-mpiexec'
    fake.write_text('#!/usr/bin/env python3\nimport json,os,sys\nfrom pathlib import Path\n'
                    'Path(os.environ["GM_TEST_CAPTURE"]).write_text(json.dumps({"argv":sys.argv[1:],'
                    '"expected":os.environ["GOLDENMODES_MPI_EXPECTED_RANKS"],'
                    '"threads":os.environ["OMP_NUM_THREADS"],"cwd":os.getcwd()}))\n')
    fake.chmod(0o755)
    capture=base/'capture.json';nodefile=base/'nodefile'
    env={**os.environ,'PBS_NODEFILE':str(nodefile),'PBS_O_WORKDIR':str(work),
         'GOLDENMODES_MPIEXEC':str(fake),'GOLDENMODES_BAT_LIBDIR':str(base),
         'GM_TEST_CAPTURE':str(capture),'ARGS':'--nChains 4 --Univariate'}
    for slots,nodes in [(35,1),(69,1),(35,2)]:
        nodefile.write_text(''.join(f'node-{i%nodes}\n' for i in range(slots)))
        result=subprocess.run(['bash',a.script],env=env,text=True,capture_output=True,timeout=10)
        assert result.returncode==0,result.stderr
        seen=json.loads(capture.read_text())
        placement=['--use-hwthread-cpus','--map-by','hwthread','--bind-to','hwthread','--report-bindings']
        assert seen['argv']==['-machinefile',str(nodefile),'-n',str(slots),*placement,str(work/'mainGM'),'--nChains','4','--Univariate'],seen
        assert seen['expected']==str(slots) and seen['threads']=='1' and Path(seen['cwd']).resolve()==work.resolve(),seen
        assert f'on {nodes} node(s)' in result.stdout
        assert 'MPI placement: '+ ' '.join(placement) in result.stdout
        assert not any('overload' in value or 'oversubscribe' in value for value in seen['argv'])
    nodefile.write_text('')
    assert subprocess.run(['bash',a.script],env=env,capture_output=True,timeout=10).returncode!=0
    # The executable must abort all ranks before starting a mismatched sampler.
    flags=[v for v in a.mpi_flags.split(';') if v]
    env={**os.environ,'GOLDENMODES_MPI_EXPECTED_RANKS':'4'}
    result=subprocess.run([a.mpiexec,'-n','3',*flags,a.executable,'--indSU3'],cwd=base,env=env,
                          text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    assert result.returncode!=0 and 'MPI launch mismatch: expected 4 MPI processes, but started 3' in result.stdout,result.stdout
print('PASS: 35/69-slot PBS launches, multi-node counting and real MPI mismatch abort')

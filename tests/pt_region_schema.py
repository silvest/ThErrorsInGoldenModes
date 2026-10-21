"""Keep sampler region cuts identical to the model's offline JSON definitions."""
import csv, io, json, subprocess, sys
rows=list(csv.DictReader(io.StringIO(subprocess.check_output([sys.argv[1], '--definitions'],text=True))))
parts=json.load(open(sys.argv[2]))['partitions']
assert len(rows)==len(parts)
for row, part in zip(rows,parts):
    cut=part['conditions'][0]
    assert len(part['conditions'])==1 and part['optional']
    assert row['parameter']==cut['parameter']==part['name']
    assert {'c':'cos','s':'sin','i':'identity'}[row['transform']]==cut['transform']
    assert ('<' if row['transform']=='i' else '>=')==cut['op']
    assert float(row['threshold'])==cut['value']
    assert row['true_label']==part['true_label'] and row['false_label']==part['false_label']
print('PASS: all model-owned sampler cuts match the offline JSON')

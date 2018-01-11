#
# Visualise a schema to validate a simple wf
#
from cerberus import Validator
from yaml import load


yml = """
inputs:
    fastqs:
        SequenceReadArchivePaired:
          forward-pattern: '*_R1.fastq.gz'
          backward-pattern: '*_R2.fastq.gz'
"""

pr = [
    'SequenceReadArchivePaired',
    {
        'schema': {
            'forward-pattern': {'type': 'string'},
            'backward-pattern': {'type': 'string'}
        }
    }
]
bm = [
    'BAM',
    {
        'schema': {
            'path': {'type': 'string'}
        }
    }
]

schema = {
    'inputs': {
        'type': 'dict',
        'keyschema': {'type': 'string'},
        'valueschema': {
            'schema': {
                pr[0]: pr[1],
                bm[0]: bm[1]
            }
        }
    }
}

doc = load(yml)

v = Validator()

print()
if not v.validate(doc, schema):
    print(doc)
    print(v.errors)
    print('failure')
else:
    print('success')

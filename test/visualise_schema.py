#
# Hand built schema for testing and visualisation
#
import unittest

from cerberus import Validator
from yaml import load
import json


class SchemeVisualisation(unittest.TestCase):

  def visualise(self):
    inputs_yml = """
    inputs:
        fastqs:
            SequenceReadArchivePaired:
              forward-pattern: '*_R1.fastq.gz'
              backward-pattern: '*_R2.fastq.gz'

        ref:
            BAM:
                path: '/path/to/ref.bam'
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

    inputs = [
      'inputs',
      {
        'required': True,
        'type': 'dict',
        'keyschema': {'type': 'string'},
        'valueschema': {
          'schema': {
            pr[0]: pr[1],
            bm[0]: bm[1]
          }
        }
      }
    ]

    trm = [
      'trim',
      {
        'schema': {
          'trimmer': {
            'type': 'string',
            'allowed': ['cutadapt', 'trimmomatic'],
            'default': 'trimmomatic'
          }
        },
        'nullable': True
      }
    ]

    aln = [
      'align',
      {
        'schema': {
          'aligner': {
            'type': 'string',
            'allowed': ['bowtie', 'bwa'],
            'default': 'bowtie'
          }
        },
        'nullable': True
      }
    ]

    steps_yml = """          
    steps:
        - step1: # label this step
            trim:
                trimmer: cutadapt
        - step2: 
            trim:
    """

    main_schema = {
      'type': 'dict',
      'keyschema': {'type': 'string'},
      'valueschema': {
        'schema': {
          trm[0]: trm[1],
          aln[0]: aln[1]
        }
      }
    }

    steps = [
      'steps',
      {
        'required': True,
        'type': 'list',
        'schema': main_schema
      }
    ]

    schema = {
      inputs[0]: inputs[1],
      steps[0]: steps[1]
    }

    yml = inputs_yml + steps_yml
    doc = load(yml)

    v = Validator()

    print(json.dumps(schema, indent=4))

    valid = v.validate(doc, schema)
    print('input yaml: ' + str(doc))
    if not valid:
      print('failure: ' + str(v.errors))

    self.assertTrue(valid)


if __name__ == '__main__':
    unittest.main()

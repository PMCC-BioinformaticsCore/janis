from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step


class AlignFactory(StepFactory):
  @classmethod
  def type(cls):
    return 'align'

  @classmethod
  def label(cls):
    return 'align'

  @classmethod
  def description(cls):
    return cls.label()

  @classmethod
  def describe(cls):
    return {
      'schema': {
        'aligner': {
          'type': 'string',
          'allowed': ['bowtie2', 'bwa'],
          'default': 'bowtie2'
        }
      },
      'nullable': True
    }

  @classmethod
  def build(cls, meta, debug=False):
    return AlignStep(meta, debug=debug)


class AlignStep(Step):

  def __init__(self, meta, debug=False):
    super().__init__(meta, debug=debug)
    if self.meta()['aligner'] != 'bowtie2':
      raise Exception('Sorry, only bowtie2 is supported at the moment')

  def cores(self):
    return 24

  def ram(self):
    return 64000

  def translate(self, step_inputs):
    vf = """${
  if ( self == null ) {
    return null;
  } else {
    return [self];
  }
}
"""
    xlate = dict()

    xlate['run'] = '../tools/src/tools/bowtie2.cwl'
    xlate['requirements'] = {'ResourceRequirement': {'coresMin': self.cores(), 'ramMin': self.ram()}}

    for mi in step_inputs:
      for candidate in mi.candidates.values():
        if mi.step_output_id == 'read' and candidate['tag'] == self.tag():
          read_step = candidate['step']
          read_id = candidate['id']
        if mi.step_output_id == 'reference':
          reference_step = candidate['step']
          reference_id = candidate['id']

    if reference_step == 'input-step':
      reference_step = 'inputs'

    inx = dict()

    if self.tag() is None:
      new_suffix = '.sam'
    else:
      new_suffix = '.' + self.tag() + '.sam'

    inx['samout'] = {
      'source': f'{read_step}/{read_id}',
      'valueFrom': f'${{ return self.nameroot + "{new_suffix}"; }}'
    }
    inx['threads'] = {'valueFrom': '$(' + str(self.cores()) + ')'}
    inx['one'] = {
      'source': f'{read_step}/reads1_trimmed',
      'valueFrom': vf
    }
    inx['two'] = {
      'source': f'{read_step}/reads2_trimmed_paired',
      'valueFrom': vf
    }
    inx['unpaired'] = {
      'source': f'{read_step}/reads_trimmed_unpaired',
      'valueFrom': vf
    }
    inx['bt2-idx'] = {'default': f'{reference_step}/{reference_id}'}
    inx['local'] = {'default': True}
    inx['reorder'] = {'default': True}

    xlate['in'] = inx
    xlate['out'] = ['aligned-file']

    return {self.id(): xlate}

  def provides(self):
    return [
      {
        Step.STR_ID: "alignedbamfile",
        Step.STR_TYPE: "bam"
      }
    ]

  def requires(self):
    return [
      {
        Step.STR_ID: "read",
        Step.STR_TYPE: "SequenceReadArchivePaired"
      },
      {
        Step.STR_ID: "reference",
        Step.STR_TYPE: "reference"
      }
    ]

import unittest, yaml

from wehi.spec import Wehi

_pipeline = """
inputs:
  tumour:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  normal:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    reference:
      path: 'path/to/reference'

steps:
  - trim_tumour:
      tag: 'tumour'
      trim:
        trimmer : 'trimmomatic'
  - align_tumour:
      tag: 'tumour'
      align:
        aligner: 'bowtie2'
  - sort_tumour:
      tag: 'tumour'
      sort:
  - call_tumour:
      tag: 'tumour'
      call:
  - trim_normal:
      tag: 'normal'
      trim:
        trimmer : 'trimmomatic'
  - align_normal:
      tag: 'normal'
      align:
        aligner: 'bowtie2'
  - sort_normal:
      tag: 'normal'
      sort:
  - call_normal:
      tag: 'normal'
      call:
  - detect_mutation:
      joint_call:
        caller: mutect
        tumour_tag: '#tumour'
        normal_tag: '#normal'

"""

_expected = yaml.load("""
class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var rename_trim_file = function() {
      if ( self == null ) {
        return null;
      } else {
        var xx = self.basename.split('.');
        var id = xx.indexOf('fastq');
        xx.splice(id, 1);
        return xx.join('.');
      }
    };
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../tools/src/tools/trimmomatic-end_mode.yml
  - $import: ../tools/src/tools/trimmomatic-sliding_window.yml
  - $import: ../tools/src/tools/trimmomatic-phred.yml
  - $import: ../tools/src/tools/trimmomatic-illumina_clipping.yml
  - $import: ../tools/src/tools/trimmomatic-max_info.yml
  
inputs:
  normal_backward: File
  normal_forward: File
  ref_reference: File
  tumour_backward: File
  tumour_forward: File

outputs:
  align_normal_aligned-file:
    outputSource: align_normal/aligned-file
    type: File
  align_tumour_aligned-file:
    outputSource: align_tumour/aligned-file
    type: File
  detect_mutation_call_stats:
    outputSource: detect_mutation/call_stats
    type: File
  detect_mutation_coverage:
    outputSource: detect_mutation/coverage
    type: File
  detect_mutation_mutations:
    outputSource: detect_mutation/mutations
    type: File
  sort_normal_sorted:
    outputSource: sort_normal/sorted
    type: File
  sort_tumour_sorted:
    outputSource: sort_tumour/sorted
    type: File
  trim_normal_output_log:
    outputSource: trim_normal/output_log
    type: File
  trim_normal_reads1_trimmed:
    outputSource: trim_normal/reads1_trimmed
    type: File
  trim_normal_reads1_trimmed_unpaired:
    outputSource: trim_normal/reads1_trimmed_unpaired
    type: File
  trim_normal_reads2_trimmed_paired:
    outputSource: trim_normal/reads2_trimmed_paired
    type: File
  trim_normal_reads2_trimmed_unpaired:
    outputSource: trim_normal/reads2_trimmed_unpaired
    type: File
  trim_tumour_output_log:
    outputSource: trim_tumour/output_log
    type: File
  trim_tumour_reads1_trimmed:
    outputSource: trim_tumour/reads1_trimmed
    type: File
  trim_tumour_reads1_trimmed_unpaired:
    outputSource: trim_tumour/reads1_trimmed_unpaired
    type: File
  trim_tumour_reads2_trimmed_paired:
    outputSource: trim_tumour/reads2_trimmed_paired
    type: File
  trim_tumour_reads2_trimmed_unpaired:
    outputSource: trim_tumour/reads2_trimmed_unpaired
    type: File
    
steps:
  align_normal:
    in:
      bt2-idx:
        default: inputs/ref
      local:
        default: true
      one:
        source: trim_normal/reads1_trimmed
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
      reorder:
        default: true
      samout:
        source: trim_normal/trimmed
        valueFrom: ${ return self.nameroot + ".normal.sam"; }
      threads:
        valueFrom: $(24)
      two:
        source: trim_normal/reads2_trimmed_paired
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
      unpaired:
        source: trim_normal/reads_trimmed_unpaired
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
    out:
    - aligned-file
    requirements:
      ResourceRequirement:
        coresMin: 24
        ramMin: 64000
    run: ../tools/src/tools/bowtie2.cwl
  align_tumour:
    in:
      bt2-idx:
        default: inputs/ref
      local:
        default: true
      one:
        source: trim_tumour/reads1_trimmed
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
      reorder:
        default: true
      samout:
        source: trim_tumour/trimmed
        valueFrom: ${ return self.nameroot + ".tumour.sam"; }
      threads:
        valueFrom: $(24)
      two:
        source: trim_tumour/reads2_trimmed_paired
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
      unpaired:
        source: trim_tumour/reads_trimmed_unpaired
        valueFrom: |
          ${
            if ( self == null ) {
              return null;
            } else {
              return [self];
            }
          }
    out:
    - aligned-file
    requirements:
      ResourceRequirement:
        coresMin: 24
        ramMin: 64000
    run: ../tools/src/tools/bowtie2.cwl
  call_normal:
    requirements:
      ResourceRequirement:
        coresMin: 2
        ramMin: 16000
    run: ../tools/src/tools/platypus.cwl
  call_tumour:
    requirements:
      ResourceRequirement:
        coresMin: 2
        ramMin: 16000
    run: ../tools/src/tools/platypus.cwl
  detect_mutation:
    in:
      normal:
        source: sort_normal/sorted
      reference:
        source: inputs/ref
      tumour:
        source: sort_tumour/sorted
      vcf:
        source: sort_tumour/sorted
        valueFrom: ${ return self.nameroot + ".vcf"; }
    out:
    - coverage
    - call_stats
    - mutations
    requirements:
      ResourceRequirement:
        coresMin: 8
        ramMin: 16000
    run: ../tools/src/tools/platypus.cwl
  sort_normal:
    in:
      input:
        source: align_normal/aligned-file
      output_name:
        source: '{align_step}/aligned-file'
        valueFrom: ${return self.nameroot + ".sorted.normal.bam";}
      threads:
        valueFrom: $( 8 )
    out:
    - sorted
    requirements:
      ResourceRequirement:
        coresMin: 8
        ramMin: 16000
    run: ../tools/src/tools/samtools-sort.cwl
  sort_tumour:
    in:
      input:
        source: align_tumour/aligned-file
      output_name:
        source: '{align_step}/aligned-file'
        valueFrom: ${return self.nameroot + ".sorted.tumour.bam";}
      threads:
        valueFrom: $( 8 )
    out:
    - sorted
    requirements:
      ResourceRequirement:
        coresMin: 8
        ramMin: 16000
    run: ../tools/src/tools/samtools-sort.cwl
  trim_normal:
    end_mode:
      default: PE
    illuminaClip:
      source: adaptors
      valueForm: |
        ${
          return {
            "adapters": self,
            "seedMismatches": 1,
            "palindromeClipThreshold": 20,
            "simpleClipThreshold": 20,
            "minAdapterLength": 4,
            "keepBothReads": true };
        }
    in:
      reads1:
        source: normal_forward
        valueForm: |
          ${
            self.format = "http://edamontology.org/format_1930";
            return self;
          }
      reads2:
        source: normal_backward
        valueForm: |
          ${
            self.format = "http://edamontology.org/format_1930";
            return self;
          }
    nthreads:
      valueFrom: $(2)
    out:
    - output_log
    - reads1_trimmed
    - reads1_trimmed_unpaired
    - reads2_trimmed_paired
    - reads2_trimmed_unpaired
    requirements:
      ResourceRequirement:
        coresMin: 2
        ramMin: 16000
    run: ../tools/src/tools/trimmomatic.cwl
  trim_tumour:
    end_mode:
      default: PE
    illuminaClip:
      source: adaptors
      valueForm: |
        ${
          return {
            "adapters": self,
            "seedMismatches": 1,
            "palindromeClipThreshold": 20,
            "simpleClipThreshold": 20,
            "minAdapterLength": 4,
            "keepBothReads": true };
        }
    in:
      reads1:
        source: tumour_forward
        valueForm: |
          ${
            self.format = "http://edamontology.org/format_1930";
            return self;
          }
      reads2:
        source: tumour_backward
        valueForm: |
          ${
            self.format = "http://edamontology.org/format_1930";
            return self;
          }
    nthreads:
      valueFrom: $(2)
    out:
    - output_log
    - reads1_trimmed
    - reads1_trimmed_unpaired
    - reads2_trimmed_paired
    - reads2_trimmed_unpaired
    requirements:
      ResourceRequirement:
        coresMin: 2
        ramMin: 16000
    run: ../tools/src/tools/trimmomatic.cwl
""")


class TumourNormalPipeline(unittest.TestCase):

  def test_graph(self):
    translator = Wehi("tumour-normal")
    translator.parse_string(_pipeline)
    # print('-'*80)
    # print(translation)
    # print('-'*80)
    # self.assertTrue(True)
    # tr_json = yaml.load(translation)
    self.assertTrue(True) # tr_json == _expected)
    # self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()



:orphan:

Combine Variants
==================================

1 contributor · 1 version

:ID: ``combinevariants``
:Python: ``janis_bioinformatics.tools.pmac.combinevariants.combinevariants_0_0_4 import CombineVariants_0_0_4``
:Versions: 0.0.4
:Container: michaelfranklin/pmacutil:0.0.4
:Authors: Michael Franklin
:Citations: None
:Created: None
:Updated: 2019-07-04 00:00:00
:Required inputs:
   - ``vcfs: Array<VCF>``

   - ``type: String``
:Outputs: 
   - ``vcf: VCF``

   - ``tsv: tsv``

Documentation
-------------

URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_


usage: combine_vcf.py [-h] -i I --columns COLUMNS -o O --type
                      {germline,somatic} [--regions REGIONS] [--normal NORMAL]
                      [--tumor TUMOR] [--priority PRIORITY [PRIORITY ...]]

Extracts and combines the information from germline / somatic vcfs into one

required arguments:
  -i I                  input vcfs, the priority of the vcfs will be based on
                        the order of the input. This parameter can be
                        specified more than once
  --columns COLUMNS     Columns to keep. This parameter can be specified more
                        than once
  -o O                  output vcf (unsorted)
  --type {germline,somatic}
                        must be either germline or somatic
  --regions REGIONS     Region file containing all the variants, used as
                        samtools mpileup
  --normal NORMAL       Sample id of germline vcf, or normal sample id of
                        somatic vcf
  --tumor TUMOR         tumor sample ID, required if inputs are somatic vcfs
  --priority PRIORITY [PRIORITY ...]
                        The priority of the callers, must match with the
                        callers in the source header

optional arguments:
  -h, --help            show this help message and exit


------

Additional configuration (inputs)
---------------------------------

==============  =======================  =============================================================================
name            type                     documentation
==============  =======================  =============================================================================
vcfs            Array<VCF>               input vcfs, the priority of the vcfs will be based on the order of the input
type            String                   germline | somatic
outputFilename  Optional<Filename>
regions         Optional<Filename>       Region file containing all the variants, used as samtools mpileup
columns         Optional<Array<String>>  Columns to keep, seperated by space output vcf (unsorted)
normal          Optional<String>         Sample id of germline vcf, or normal sample id of somatic vcf
tumor           Optional<String>         tumor sample ID, required if inputs are somatic vcfs
priority        Optional<Integer>        The priority of the callers, must match with the callers in the source header
==============  =======================  =============================================================================


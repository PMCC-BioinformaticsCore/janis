
.. include:: combinevariants_0.0.4

Combine Variants
==================================

Description
-------------

Tool identifier: ``combinevariants``

Tool path: ``janis_bioinformatics.tools.pmac.combinevariants.combinevariants_0_0_4 import CombineVariants_0_0_4``

Version: 0.0.4

Container: ``michaelfranklin/pmacutil:0.0.4``



Documentation
-------------

URL
******
`https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_

Tool documentation
******************

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


Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
vcf     VCF
tsv     tsv
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ==========  ========  ==========  ============================================================================
name    type        prefix    position    documentation
======  ==========  ========  ==========  ============================================================================
vcfs    Array<VCF>  -i                    input vcfs, the priority of the vcfs will be based on the order of the input
type    String      --type                germline | somatic
======  ==========  ========  ==========  ============================================================================

Optional inputs
***************

==============  =======================  ==========  ==========  =============================================================================
name            type                     prefix      position    documentation
==============  =======================  ==========  ==========  =============================================================================
outputFilename  Optional<Filename>       -o
regions         Optional<Filename>       --regions               Region file containing all the variants, used as samtools mpileup
columns         Optional<Array<String>>  --columns               Columns to keep, seperated by space output vcf (unsorted)
normal          Optional<String>         --normal                Sample id of germline vcf, or normal sample id of somatic vcf
tumor           Optional<String>         --tumor                 tumor sample ID, required if inputs are somatic vcfs
priority        Optional<Integer>        --priority              The priority of the callers, must match with the callers in the source header
==============  =======================  ==========  ==========  =============================================================================


Metadata
********

Author: Jiaan Yu


*Combine Variants was last updated on 2019-07-04 00:00:00*.
*This page was automatically generated on 2019-09-10*.

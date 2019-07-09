
Combine Variants
==================================
Tool identifier: ``combinevariants``

Tool path: ``from janis_bioinformatics.tools.peter maccallum cancer centre import CombineVariants_0_0_4``

Documentation
-------------

Docker
******
``michaelfranklin/pmacutil:0.0.4``

URL
******
*No URL to the documentation was provided*

Docstring
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

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

Author: **Unknown**


*Combine Variants was last updated on **Unknown***.
*This page was automatically generated on 2019-07-09*.

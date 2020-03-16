:orphan:

Split Multiple Alleles
=========================================

*0 contributors · 1 version*

:ID: ``SplitMultiAllele``
:Python: ``janis_bioinformatics.tools.common.splitmultiallele import SplitMultiAllele``
:Versions: v0.5772
:Container: heuermh/vt
:Authors: 
:Citations: None
:Created: None
:Updated: None
:Required inputs:
   - ``vcf: CompressedVCF``

   - ``reference: FastaWithIndexes``
:Outputs: 
   - ``out: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

------

Arguments
----------

========================================  ========  ==========  ===============
value                                     prefix      position  documentation
========================================  ========  ==========  ===============
zcat                                                         0
|                                                            1
sed 's/ID=AD,Number=./ID=AD,Number=R/' <                     2
|                                                            4
vt decompose -s - -o -                                       5
|                                                            6
vt normalize -n -q - -o -                                    7
|                                                            9
sed 's/ID=AD,Number=./ID=AD,Number=1/'                      10
========================================  ========  ==========  ===============

Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
vcf             CompressedVCF                          3
reference       FastaWithIndexes    -r                 8
outputFilename  Optional<Filename>  >                 10
==============  ==================  ========  ==========  ===============


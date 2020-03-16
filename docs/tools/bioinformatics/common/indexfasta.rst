:orphan:

Index Fasta reference
==================================

*0 contributors · 1 version*

:ID: ``IndexFasta``
:Python: ``janis_bioinformatics.tools.common.indexfasta import IndexFasta``
:Versions: 1.0.0
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``reference: Fasta``
:Outputs: 
   - ``bwa: FastaFai``

   - ``samtools: FastaFai``

   - ``dict: FastDict``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

===============================  =========================================
BWA-Index                        ``bwaIndex/v0.7.15``
SamTools: faidx                  ``SamToolsIndex/1.9.0``
GATK4: CreateSequenceDictionary  ``Gatk4CreateSequenceDictionary/4.1.3.0``
===============================  =========================================

------

Additional configuration (inputs)
---------------------------------

=========  ======  ===============
name       type    documentation
=========  ======  ===============
reference  Fasta
=========  ======  ===============

.

:orphan:

GATK4: CreateSequenceDictionary
===============================================================

*1 contributor · 3 versions*

:ID: ``Gatk4CreateSequenceDictionary``
:Python: ``janis_bioinformatics.tools.gatk4.createsequencedictionary.versions import Gatk4CreateSequenceDictionary_4_1_4``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: TBD
:Created: 2020-02-14
:Updated: 2020-02-14
:Required inputs:
   - ``reference: Fasta``
:Outputs: 
   - ``out: FastDict``

Documentation
-------------

URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard->`_

Creates a sequence dictionary for a reference sequence.  This tool creates a sequence dictionary file (with ".dict"
extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools.
The output file contains a header but no SAMRecords, and the header contains only sequence records.

The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

Usage example:

    java -jar picard.jar CreateSequenceDictionary \
        R=reference.fasta \
        O=reference.dict

------

None

Additional configuration (inputs)
---------------------------------

=========  ======  ===========  ==========  =================================================
name       type    prefix       position    documentation
=========  ======  ===========  ==========  =================================================
reference  Fasta   --REFERENCE              (-R) Input reference fasta or fasta.gz  Required.
=========  ======  ===========  ==========  =================================================


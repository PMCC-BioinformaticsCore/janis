:orphan:

SamTools: Index
===============================

*1 contributor · 2 versions*

:ID: ``SamToolsIndex``
:Python: ``janis_bioinformatics.tools.samtools.index.versions import SamToolsIndex_1_9``
:Versions: 1.9.0, 1.7.0
:Container: quay.io/biocontainers/samtools:1.9--h8571acd_11
:Authors: Michael Franklin
:Citations: None
:Created: 2019-12-17
:Updated: 2019-12-17
:Required inputs:
   - ``bam: BAM``
:Outputs: 
   - ``out: IndexedBam``

Documentation
-------------

URL: `http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS <http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS>`_

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

------

Arguments
----------

=======  ========  ==========  =========================
value    prefix      position  documentation
=======  ========  ==========  =========================
-b                          4  Output in the BAM format.
=======  ========  ==========  =========================

Additional configuration (inputs)
---------------------------------

=======  =================  ========  ==========  ===============
name     type               prefix      position  documentation
=======  =================  ========  ==========  ===============
bam      BAM                                  10
threads  Optional<Integer>  -@                10
=======  =================  ========  ==========  ===============


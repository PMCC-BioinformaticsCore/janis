:orphan:

BWA-Index
====================

*1 contributor · 1 version*

:ID: ``bwaIndex``
:Python: ``janis_bioinformatics.tools.bwa.index.versions import BwaIndex_0_7_15``
:Versions: v0.7.15
:Container: biocontainers/bwa:v0.7.15_cv3
:Authors: Michael Franklin
:Citations: The BWA-MEM algorithm has not been published yet.
:Created: 2020-02-14
:Updated: 2020-02-14
:Required inputs:
   - ``reference: Fasta``
:Outputs: 
   - ``out: FastaFai``

Documentation
-------------

URL: `http://bio-bwa.sourceforge.net/bwa.shtml#3 <http://bio-bwa.sourceforge.net/bwa.shtml#3>`_

bwa - Burrows-Wheeler Alignment Tool
Index database sequences in the FASTA format.

Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
         `-a div' do not work not for long genomes.

------

None

Additional configuration (inputs)
---------------------------------

=========  =================  ========  ==========  =======================================================================
name       type               prefix      position  documentation
=========  =================  ========  ==========  =======================================================================
reference  Fasta                                 1
blockSize  Optional<Integer>  -b                    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
algorithm  Optional<String>   -a                    BWT construction algorithm: bwtsw, is or rb2 [auto]
=========  =================  ========  ==========  =======================================================================


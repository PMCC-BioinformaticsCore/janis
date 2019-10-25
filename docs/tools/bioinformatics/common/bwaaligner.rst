:orphan:

Align and sort reads
=================================

0 contributors · 1 version

:ID: ``BwaAligner``
:Python: ``janis_bioinformatics.tools.common.bwaaligner import BwaAligner``
:Versions: 1.0.0
:Authors: 
:Citations: 
:Created: 2018-12-24
:Updated: None
:Required inputs:
   - ``sampleName: String``

   - ``reference: FastaWithDict``

   - ``fastq: FastqGzPair``
:Outputs: 
   - ``out: BamPair``

Documentation
-------------

URL: *No URL to the documentation was provided*

Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam

Embedded Tools
***************

=======================  =================================
Cutadapt                 ``cutadapt/1.18``
Bwa mem + Samtools View  ``BwaMemSamtoolsView/0.7.17|1.9``
GATK4: SortSAM           ``gatk4sortsam/4.1.3.0``
=======================  =================================

------

Additional configuration (inputs)
---------------------------------

=============================  =================  ===============
name                           type               documentation
=============================  =================  ===============
sampleName                     String
reference                      FastaWithDict
fastq                          FastqGzPair
cutadapt_adapter               Optional<String>
cutadapt_adapter_g             Optional<String>
cutadapt_removeMiddle5Adapter  Optional<String>
cutadapt_removeMiddle3Adapter  Optional<String>
cutadapt_qualityCutoff         Optional<Integer>
cutadapt_minReadLength         Optional<Integer>
bwamem_markShorterSplits       Optional<Boolean>
sortsam_sortOrder              Optional<String>
sortsam_createIndex            Optional<Boolean>
sortsam_validationStringency   Optional<String>
sortsam_maxRecordsInRam        Optional<Integer>
sortsam_tmpDir                 Optional<String>
=============================  =================  ===============

.

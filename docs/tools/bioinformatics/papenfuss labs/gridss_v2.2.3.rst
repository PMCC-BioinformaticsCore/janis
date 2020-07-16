:orphan:

Gridss
===============

*1 contributor · 4 versions*

GRIDSS: the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements. 
GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina 
sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide 
break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design, 
any step (such as split read identification) can be replaced by another implementation that also outputs 
using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant 
pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback 
from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.papenfuss.gridss.gridss import Gridss_2_2_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gridss_step",
           Gridss_2_2_3(
               reference=None,
               bams=None,
           )
       )
       wf.output("vcf", source=gridss_step.vcf)
       wf.output("assembly", source=gridss_step.assembly)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for gridss:

.. code-block:: bash

   # user inputs
   janis inputs gridss > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       reference: reference.fasta




5. Run gridss with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       gridss





Information
------------


:ID: ``gridss``
:URL: `https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation <https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation>`_
:Versions: v2.6.2, v2.5.1-dev, v2.4.0, v2.2.3
:Container: gridss/gridss:v2.2.3
:Authors: Michael Franklin
:Citations: Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. Genome Research, 2017 doi: 10.1101/gr.222109.117
:DOI: 10.1101/gr.222109.117
:Created: 2019-06-19
:Updated: 2019-07-03



Outputs
-----------

========  ==========  ===============
name      type        documentation
========  ==========  ===============
vcf       VCF
assembly  IndexedBam
========  ==========  ===============



Additional configuration (inputs)
---------------------------------

=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================
name                       type                prefix                         position    documentation
=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================
reference                  FastaWithIndexes    REFERENCE_SEQUENCE=
bams                       Array<IndexedBam>   INPUT=                                     (I=File Coordinate-sorted input BAM file. Default value: null. This option may be specified 0 or more times.
outputFilename             Optional<Filename>  OUTPUT=                                    (O=) VCF structural variation calls. Required.
assemblyFilename           Optional<Filename>  ASSEMBLY=                                  Breakend assemblies which have undergone split read identification Required.
inputLabel                 Optional<String>    INPUT_LABEL=                               Input label. Variant calling evidence breakdowns are reported for each label. Default labels correspond to INPUT filenames. When specifying labels, labels must be provided for all input files. Default value: null. This option may be specified 0 or more times.
inputMaxFragmentSize       Optional<Integer>   INPUT_MAX_FRAGMENT_SIZE=                   Per input maximum concordant fragment size. Default value: null. This option may be specified 0 or more times.
inputMinFragmentSize       Optional<Integer>   INPUT_MIN_FRAGMENT_SIZE=                   Per input minimum concordant fragment size. Default value: null. This option may be specified 0 or more times.
readPairConcordantPercent  Optional<Float>     READ_PAIR_CONCORDANT_PERCENT=              Percent of read pairs considered concorant (0.0-1.0). If this is unset, the SAM proper pair flag is used to determine whether a read is discordantly aligned. Explicit fragment size specification overrides this setting. Default value: 0.995. This option can be set to 'null' to clear the default value.
blacklist                  Optional<bed>       BLACKLIST=                                 (BL=File) BED blacklist of regions to ignore. Assembly of regions such as high-coverage centromeric repeats is slow, and if such regions are to be filtered in downstream analysis anyway, blacklisting those region will improve runtime performance. For human WGS, the ENCODE DAC blacklist is recommended. Default value: null.
configurationFile          Optional<File>      CONFIGURATION_FILE=                        (C=File) gridss configuration file containing overrides Default value: null.
workerThreads              Optional<Integer>   WORKER_THREADS=                            (THREADS=Integer  Number of worker threads to spawn. Defaults to number of cores available. Note that I/O threads are not included in this worker thread count so CPU usage can be higher than the number of worker thread. Default value: 6. This option can be set to 'null' to clear the default value.
workingDir                 Optional<String>    WORKING_DIR=                               Directory to place intermediate results directories. Default location is the same directory as the associated input or output file. Default value: null.
ignoreDuplicates           Optional<Boolean>   IGNORE_DUPLICATES=                         Ignore reads marked as duplicates. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false}
=========================  ==================  =============================  ==========  ===================================================================================================================================================================================================================================================================================================================================

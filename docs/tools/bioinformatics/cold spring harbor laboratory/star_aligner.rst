:orphan:

STAR Aligner
===========================

*1 contributor · 1 version*

Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
Spliced Transcripts Alignment to a Reference (c) Alexander Dobin, 2009-2019

For more details see:
<https://github.com/alexdobin/STAR>
<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
            


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.star.versions import StarAligner_2_7_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "star_aligner_step",
           StarAligner_2_7_1(

           )
       )
       wf.output("logFinalOut", source=star_aligner_step.logFinalOut)
   wf.output("logOut", source=star_aligner_step.logOut)
   wf.output("logProgressOut", source=star_aligner_step.logProgressOut)
   wf.output("sjOutTab", source=star_aligner_step.sjOutTab)
   wf.output("out", source=star_aligner_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for star_aligner:

.. code-block:: bash

   # user inputs
   janis inputs star_aligner > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run star_aligner with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       star_aligner





Information
------------


:ID: ``star_aligner``
:URL: `https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>`_
:Versions: v2.7.1a
:Container: quay.io/biocontainers/star:2.7.3a--0
:Authors: Jiaan Yu
:Citations: Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics (Oxford, England), 29(1), 15–21. https://doi.org/10.1093/bioinformatics/bts635
:DOI: https://doi.org/10.1093/bioinformatics/bts635
:Created: 2020-04-02
:Updated: 2020-04-02



Outputs
-----------

==============  ======  ===============
name            type    documentation
==============  ======  ===============
logFinalOut     File
logOut          File
logProgressOut  File
sjOutTab        File
out             BAM
==============  ======  ===============



Additional configuration (inputs)
---------------------------------

=================  ========================  ===================  =================================================================  ========================================================================================================================================================================================================================================================================================================================
name               type                      prefix               position                                                           documentation
=================  ========================  ===================  =================================================================  ========================================================================================================================================================================================================================================================================================================================
help               Optional<Boolean>         --help                                                                                  help page
runThreadN         Optional<Integer>         --runThreadN                                                                            int: number of threads to run STAR. Default: 1.
genomeDir          Optional<Directory>       --genomeDir                                                                             string: path to the directory where genome files are stored (for –runMode alignReads) or will be generated (for –runMode generateGenome). Default: ./GenomeDir
readFilesIn        Optional<Array<FastqGz>>  --readFilesIn                                                                           string(s): paths to files that contain input read1 (and, if needed, read2). Default: Read1,Read2.
outFileNamePrefix  Optional<Filename>        --outFileNamePrefix  <janis_core.types.common_data_types.String object at 0x1033f5ac8>  string: output files name prefix (including full or relative path). Can only be defined on the command line.
outSAMtype         Optional<Array<String>>   --outSAMtype                                                                            strings: type of SAM/BAM output. 1st word: "BAM": outputBAMwithoutsorting, "SAM": outputSAMwithoutsorting, "None": no SAM/BAM output. 2nd,3rd: "Unsorted": standard unsorted. "SortedByCoordinate": sorted by coordinate. This option will allocate extra memory for sorting which can be specified by –limitBAMsortRAM.
outSAMunmapped     Optional<String>          --outSAMunmapped                                                                        string(s): output of unmapped reads in the SAM format
outSAMattributes   Optional<String>          --outSAMattributes                                                                      string: a string of desired SAM attributes, in the order desired for the output SAM
readFilesCommand   Optional<String>          --readFilesCommand                                                                      string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
=================  ========================  ===================  =================================================================  ========================================================================================================================================================================================================================================================================================================================

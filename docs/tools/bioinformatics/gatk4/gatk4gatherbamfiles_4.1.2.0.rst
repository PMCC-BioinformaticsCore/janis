:orphan:

GATK4: GatherBamFiles
===========================================

*0 contributors · 3 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.gatherbamfiles.versions import Gatk4GatherBamFiles_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4gatherbamfiles_step",
           Gatk4GatherBamFiles_4_1_2(

           )
       )
       wf.output("out", source=gatk4gatherbamfiles_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GatherBamFiles:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GatherBamFiles > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run Gatk4GatherBamFiles with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GatherBamFiles





Information
------------


:ID: ``Gatk4GatherBamFiles``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================
name                   type                     prefix                   position    documentation
=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================
javaOptions            Optional<Array<String>>
compression_level      Optional<Integer>                                             Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
bams                   Optional<Array<BAM>>     --INPUT                              (-I) Two or more BAM files or text files containing lists of BAM files (one per line). This argument must be specified at least once. Required.
outputFilename         Optional<Filename>       --OUTPUT                             (-O) The output BAM file to write to. Required.
arguments_file         Optional<File>           --arguments_file                     read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
create_index           Optional<Boolean>        --CREATE_INDEX                       Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. Possible values: {true, false}
create_md5_file        Optional<Boolean>        --CREATE_MD5_FILE                    Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. Possible values: {true, false}
ga4gh_client_secrets   Optional<Boolean>        --GA4GH_CLIENT_SECRETS               Default value: client_secrets.json.
help                   Optional<Boolean>        --help                               (-h) display the help message Default value: false. Possible values: {true, false}
max_records_in_ram     Optional<Integer>        --MAX_RECORDS_IN_RAM                 When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.  Default value: 500000.
quiet                  Optional<Boolean>        --QUIET                              Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
reference_sequence     Optional<File>           --REFERENCE_SEQUENCE                 (-R) Reference sequence file. Default value: null.
tmp_dir                Optional<File>           --TMP_DIR                            One or more directories with space available to be used by this program for temporary storage of working files  This argument may be specified 0 or more times. Default value: null.
use_jdk_deflater       Optional<Boolean>        --USE_JDK_DEFLATER                   (-use_jdk_deflater)  Use the JDK Deflater instead of the Intel Deflater for writing compressed output  Default value: false. Possible values: {true, false}
use_jdk_inflater       Optional<Boolean>        --USE_JDK_INFLATER                   (-use_jdk_inflater)  Use the JDK Inflater instead of the Intel Inflater for reading compressed input  Default value: false. Possible values: {true, false}
validation_stringency  Optional<Boolean>        --VALIDATION_STRINGENCY              Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.  Default value: STRICT. Possible values: {STRICT, LENIENT, SILENT}
verbosity              Optional<Boolean>        --VERBOSITY                          Control verbosity of logging. Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                Optional<Boolean>        --version                            display the version number for this tool Default value: false. Possible values: {true, false}
showhidden             Optional<Boolean>        --showHidden                         (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
=====================  =======================  =======================  ==========  ============================================================================================================================================================================================================================================================================================================

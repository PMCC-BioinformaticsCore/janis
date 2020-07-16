:orphan:

GATK4: GatherBQSRReports
=================================================

*0 contributors · 4 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.gatherbqsrreports.versions import Gatk4GatherBQSRReports_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4gatherbqsrreports_step",
           Gatk4GatherBQSRReports_4_1_2(

           )
       )
       wf.output("out", source=gatk4gatherbqsrreports_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GatherBQSRReports:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GatherBQSRReports > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run Gatk4GatherBQSRReports with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GatherBQSRReports





Information
------------


:ID: ``Gatk4GatherBQSRReports``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============



Additional configuration (inputs)
---------------------------------

==========================  =======================  ================================  ==========  ======================================================================================================================================
name                        type                     prefix                            position    documentation
==========================  =======================  ================================  ==========  ======================================================================================================================================
javaOptions                 Optional<Array<String>>
compression_level           Optional<Integer>                                                      Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
reports                     Optional<Array<tsv>>     --input                                       (-I) List of scattered BQSR report files This argument must be specified at least once. Required.
outputFilename              Optional<Filename>       --output                                      (-O) File to output the gathered file to Required.
arguments_file              Optional<File>           --arguments_file                              read one or more arguments files and add them to the command line This argument may be specified 0 or more times. Default value: null.
gatkConfigFile              Optional<String>         --gatk-config-file                            A configuration file to use with the GATK. Default value: null.
gcsMaxRetries               Optional<Integer>        --gcs-max-retries                             (-gcs-retries)  If the GCS bucket channel errors out, how many times it will attempt to re-initiate the connection  Default value: 20.
gcsProjectForRequesterPays  Optional<String>         --gcs-project-for-requester-pays              Project to bill when accessing 'requester pays' buckets. If unset, these buckets cannot be accessed.  Default value: .
help                        Optional<Boolean>        --help                                        (-h) display the help message Default value: false. Possible values: {true, false}
quiet                       Optional<Boolean>        --QUIET                                       Whether to suppress job-summary info on System.err. Default value: false. Possible values: {true, false}
tmpDir                      Optional<Boolean>        --tmp-dir                                     Temp directory to use. Default value: null.
useJdkDeflater              Optional<Boolean>        --use-jdk-deflater                            (-jdk-deflater)  Whether to use the JdkDeflater (as opposed to IntelDeflater)  Default value: false. Possible values: {true, false}
useJdkInflater              Optional<Boolean>        --use-jdk-inflater                            (-jdk-inflater)  Whether to use the JdkInflater (as opposed to IntelInflater)  Default value: false. Possible values: {true, false}
verbosity                   Optional<Boolean>        --verbosity                                   (-verbosity)  Control verbosity of logging.  Default value: INFO. Possible values: {ERROR, WARNING, INFO, DEBUG}
version                     Optional<Boolean>        --version                                     display the version number for this tool Default value: false. Possible values: {true, false}
showhidden                  Optional<Boolean>        --showHidden                                  (-showHidden)  display hidden arguments  Default value: false. Possible values: {true, false}
==========================  =======================  ================================  ==========  ======================================================================================================================================

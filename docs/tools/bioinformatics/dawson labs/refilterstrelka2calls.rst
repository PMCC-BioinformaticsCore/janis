:orphan:

Refilter Strelka2 Variant Calls
=======================================================

*1 contributor · 1 version*

Usage: filterStrelkaCalls.R [options]



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.refilterstrelka2calls.latest import RefilterStrelka2CallsLatest

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "refilterstrelka2calls_step",
           RefilterStrelka2CallsLatest(
               inputFiles=None,
               MQ=None,
               DP=None,
               EVS=None,
               RPRS=None,
               threads=None,
               outputFolder=None,
           )
       )
       wf.output("out", source=refilterstrelka2calls_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for refilterStrelka2Calls:

.. code-block:: bash

   # user inputs
   janis inputs refilterStrelka2Calls > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputFiles:
       - inputFiles_0.vcf.gz
       - inputFiles_1.vcf.gz




5. Run refilterStrelka2Calls with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       refilterStrelka2Calls





Information
------------


:ID: ``refilterStrelka2Calls``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.6
:Container: shollizeck/dawsontoolkit:0.1.6.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-19
:Updated: 2019-10-25



Outputs
-----------

======  ==========  =================
name    type        documentation
======  ==========  =================
out     Array<VCF>  To determine type
======  ==========  =================



Additional configuration (inputs)
---------------------------------

============  ===========================  =============  ==========  ========================================================================
name          type                         prefix         position    documentation
============  ===========================  =============  ==========  ========================================================================
inputFiles    Array<CompressedIndexedVCF>  -i                         comma seperated list of vcfs
MQ            Integer                      --mq                       minimum mapping quality for a variant to be accepted (default: 15)
DP            Integer                      --dp                       minimum depth of coverage for a variant to be accepted (default: 10)
EVS           Integer                      --evs                      minimum phred scaled evidence for a variant to be accepted (default: 20)
RPRS          Integer                      --rprs                     minimum phred scaled evidence for a variant to be accepted (default: 20)
threads       Integer                      -t                         amount of threads to use for parallelization (default: 5)
outputFolder  String                       -o                         Name of the normal sample (default: infered from all sample names)
interval      Optional<String>             -L                         interval to call on (default: everything)
normalName    Optional<String>             -n                         Name of the normal sample (default: infered from all sample names)
sampleNames   Optional<Array<String>>      --sampleNames              Name of the normal sample (default: infered from all sample names)
============  ===========================  =============  ==========  ========================================================================

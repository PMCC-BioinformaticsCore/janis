:orphan:

GATK4: MergeMutectStats
===============================================

*1 contributor · 3 versions*

TBD

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mergemutectstats.versions import Gatk4MergeMutectStats_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mergemutectstats_step",
           Gatk4MergeMutectStats(
               statsFiles=None,
           )
       )
       wf.output("out", source=gatk4mergemutectstats_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4MergeMutectStats:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4MergeMutectStats > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       statsFiles:
       - statsFiles_0.txt
       - statsFiles_1.txt




5. Run Gatk4MergeMutectStats with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4MergeMutectStats





Information
------------


:ID: ``Gatk4MergeMutectStats``
:URL: `TBD <TBD>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09



Outputs
-----------

======  ========  ========================
name    type      documentation
======  ========  ========================
out     TextFile  Merged callability stats
======  ========  ========================



Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  =================
name            type                prefix      position  documentation
==============  ==================  ========  ==========  =================
statsFiles      Array<TextFile>     --stats            0  Callability stats
mergedStatsOut  Optional<Filename>  -O                 1
==============  ==================  ========  ==========  =================

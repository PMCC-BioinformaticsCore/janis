:orphan:

GATK4: LearnReadOrientationModel
=================================================================

*1 contributor · 3 versions*

TBD

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.learnreadorientationmodel.versions import Gatk4LearnReadOrientationModel_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4learnreadorientationmodel_step",
           Gatk4LearnReadOrientationModel(
               f1r2CountsFiles=None,
           )
       )
       wf.output("out", source=gatk4learnreadorientationmodel_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4LearnReadOrientationModel:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4LearnReadOrientationModel > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       f1r2CountsFiles:
       - f1r2CountsFiles_0
       - f1r2CountsFiles_1




5. Run Gatk4LearnReadOrientationModel with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4LearnReadOrientationModel





Information
------------


:ID: ``Gatk4LearnReadOrientationModel``
:URL: `TBD <TBD>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09



Outputs
-----------

======  =================  =============================================================
name    type               documentation
======  =================  =============================================================
out     CompressedTarFile  Model file containing information about fragment orientations
======  =================  =============================================================



Additional configuration (inputs)
---------------------------------

===============  ========================  ===================  ==========  =======================================================
name             type                      prefix                 position  documentation
===============  ========================  ===================  ==========  =======================================================
f1r2CountsFiles  Array<CompressedTarFile>  -I                            0  Counts for the read orientation of fragments
numEmIterations  Optional<Integer>         --num-em-iterations           1  Amount of iterations for the em process before it bails
modelFileOut     Optional<Filename>        -O                            3
===============  ========================  ===================  ==========  =======================================================

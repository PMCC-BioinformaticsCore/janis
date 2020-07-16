:orphan:

Add Sym to DepthOfCoverage
====================================================

*1 contributor · 2 versions*

usage: add_sym_to_DepthOfCoverage.py [-h] -i INPUT -o OUTPUT -bed BED

Performance summary of bam

optional arguments:
  -h, --help  show this help message and exit
  -i INPUT    Gatk3 DepthOfCoverage interval_summary output
  -o OUTPUT   Output file name
  -bed BED    Annotated bed file
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addsymtodepthofcoverage.versions import AddSymToDepthOfCoverage_dev

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addsymtodepthofcoverage_step",
           AddSymToDepthOfCoverage_dev(
               inputFile=None,
               bed=None,
           )
       )
       wf.output("out", source=addsymtodepthofcoverage_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for addSymToDepthOfCoverage:

.. code-block:: bash

   # user inputs
   janis inputs addSymToDepthOfCoverage > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bed: bed.bed
       inputFile: inputFile




5. Run addSymToDepthOfCoverage with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       addSymToDepthOfCoverage





Information
------------


:ID: ``addSymToDepthOfCoverage``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance>`_
:Versions: dev, 0.0.7
:Container: jyu/pmacutil:dev
:Authors: Jiaan Yu
:Citations: None
:Created: None
:Updated: 2020-04-09 00:00:00



Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  =============================================
name            type                prefix    position    documentation
==============  ==================  ========  ==========  =============================================
inputFile       File                -i                    Gatk3 DepthOfCoverage interval_summary output
bed             bed                 -bed                  Annotated bed file
outputFilename  Optional<Filename>  -o                    Output file name
==============  ==================  ========  ==========  =============================================

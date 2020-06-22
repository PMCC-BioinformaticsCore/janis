:orphan:

Sequenza: seqz binning
====================================

*2 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.sequenza.seqz_binning.versions import SequenzaBinning_3_0_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "seqzbinning_step",
           SequenzaBinning_3_0_0(
               seqz=None,
               window=None,
           )
       )
       wf.output("out", source=seqzbinning_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SeqzBinning:

.. code-block:: bash

   # user inputs
   janis inputs SeqzBinning > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       seqz: seqz
       window: 0




5. Run SeqzBinning with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SeqzBinning





Information
------------


:ID: ``SeqzBinning``
:URL: `http://www.cbs.dtu.dk/biotools/sequenza/ <http://www.cbs.dtu.dk/biotools/sequenza/>`_
:Versions: 3.0.0, 2.2.0.9000
:Container: sequenza/sequenza:3.0.0
:Authors: mumbler, evanwehi
:Citations: None
:Created: 2019-12-16
:Updated: 2019-12-16



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============



Additional configuration (inputs)
---------------------------------

===============  ==================  ========  ==========  ===================================================================
name             type                prefix      position  documentation
===============  ==================  ========  ==========  ===================================================================
seqz             File                --seqz             2  A seqz file.
window           Integer             --window           4  Window size used for binning the original seqz file. Default is 50.
output_filename  Optional<Filename>  -o                 6  Output file "-" for STDOUT
===============  ==================  ========  ==========  ===================================================================

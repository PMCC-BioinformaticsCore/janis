:orphan:

Tar (archive)
===================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.tar import Tar

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "tar_step",
           Tar(
               files=None,
           )
       )
       wf.output("out", source=tar_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Tar:

.. code-block:: bash

   # user inputs
   janis inputs Tar > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       files:
       - files_0
       - files_1




5. Run Tar with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Tar





Information
------------


:ID: ``Tar``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  =======  ===============
name    type     documentation
======  =======  ===============
out     TarFile
======  =======  ===============



Additional configuration (inputs)
---------------------------------

==============  =====================  ========  ==========  ===============
name            type                   prefix      position  documentation
==============  =====================  ========  ==========  ===============
files           Array<File>                               2
files2          Optional<Array<File>>                     3
outputFilename  Optional<Filename>                        1
==============  =====================  ========  ==========  ===============

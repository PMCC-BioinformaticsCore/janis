:orphan:

Tar (unarchive)
=======================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.untar import Untar

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "untar_step",
           Untar(
               tarfile=None,
           )
       )
       wf.output("out", source=untar_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for untar:

.. code-block:: bash

   # user inputs
   janis inputs untar > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       tarfile: tarfile




5. Run untar with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       untar





Information
------------


:ID: ``untar``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     Array<File>
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

=======  =======  ========  ==========  ===============
name     type     prefix      position  documentation
=======  =======  ========  ==========  ===============
tarfile  TarFile                     0
=======  =======  ========  ==========  ===============

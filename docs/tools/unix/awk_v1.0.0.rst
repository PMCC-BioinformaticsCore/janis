:orphan:

Awk
=========

*0 contributors · 1 version*

run an awk script


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.awk import Awk

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "awk_step",
           Awk(
               script=None,
               input_files=None,
           )
       )
       wf.output("out", source=awk_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for awk:

.. code-block:: bash

   # user inputs
   janis inputs awk > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       input_files:
       - input_files_0
       - input_files_1
       script: script




5. Run awk with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       awk





Information
------------


:ID: ``awk``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============



Additional configuration (inputs)
---------------------------------

===========  ===========  ========  ==========  ===============
name         type         prefix      position  documentation
===========  ===========  ========  ==========  ===============
script       File         -f                 1
input_files  Array<File>                     2
===========  ===========  ========  ==========  ===============

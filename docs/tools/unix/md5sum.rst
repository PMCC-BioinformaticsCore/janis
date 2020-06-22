:orphan:

MD5 Sum
================

*0 contributors · 1 version*

Compute the MD5 message digest of the given file.


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.md5sum import MD5Sum

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "md5sum_step",
           MD5Sum(
               input_file=None,
           )
       )
       wf.output("out", source=md5sum_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for md5sum:

.. code-block:: bash

   # user inputs
   janis inputs md5sum > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       input_file: input_file




5. Run md5sum with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       md5sum





Information
------------


:ID: ``md5sum``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: 2020-06-09 00:00:00



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============



Additional configuration (inputs)
---------------------------------

==========  ======  ========  ==========  ===============
name        type    prefix      position  documentation
==========  ======  ========  ==========  ===============
input_file  File                       1
==========  ======  ========  ==========  ===============

:orphan:

Java compiler
============================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_

Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.compile import Compile

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "javacompiler_step",
           javacompiler(
               file=None,
           )
       )
       wf.output("out", source=javacompiler_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for javacompiler:

.. code-block:: bash

   # user inputs
   janis inputs javacompiler > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file




5. Run javacompiler with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       javacompiler





Information
------------


:ID: ``javacompiler``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: openjdk:8
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============



Additional configuration (inputs)
---------------------------------

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
file    File                       1
======  ======  ========  ==========  ===============

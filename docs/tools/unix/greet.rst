:orphan:

Greet
=============

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_

Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.greet import Greet

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "greet_step",
           Greet(
               name=None,
           )
       )
       wf.output("out", source=greet_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for greet:

.. code-block:: bash

   # user inputs
   janis inputs greet > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       name: <value>




5. Run greet with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       greet





Information
------------


:ID: ``greet``
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

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
name    String                     1
======  ======  ========  ==========  ===============

:orphan:

Concatenate
=================

*0 contributors · 1 version*

The cat utility reads files sequentially, writing them to the standard output. The file operands are processed in command-line order. If file is a single dash (`-') or absent,cat reads from the standard input. If file is a UNIX domain socket, cat connects to it and then reads it until EOF. This complements the UNIX domain binding capability available in inetd(8).

Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.cat import Cat

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cat_step",
           Cat(

           )
       )
       wf.output("out", source=cat_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for cat:

.. code-block:: bash

   # user inputs
   janis inputs cat > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run cat with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       cat





Information
------------


:ID: ``cat``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: 2019-07-26 00:00:00



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============



Additional configuration (inputs)
---------------------------------

======================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
name                    type                   prefix      position  documentation
======================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
file                    Optional<File>
files                   Optional<Array<File>>                     1
number_output           Optional<Boolean>      -n                    Number the output lines, starting at 1.
number_non_blank        Optional<Boolean>      -b                    Number the non-blank output lines, starting at 1.
disable_output_buffer   Optional<Boolean>      -u                    Disable output buffering.
squeeze                 Optional<Boolean>      -s                    Squeeze multiple adjacent empty lines, causing the output to be single spaced.
display_nonprint_chars  Optional<Boolean>      -e                    Display non-printing characters (see the -v option), and display a dollar sign (`$') at the end of each line.
display_non             Optional<Boolean>      -t                    Display non-printing characters (see the -v option), and display tab characters as `^I'.
number_non_blank        Optional<Boolean>      -v                    Display non-printing characters so they are visible.  Control characters print as `^X' for control-X; the delete character (octal 0177) prints as `^?'.  Non-ASCII characters (with the high bit set) are printed as `M-' (for meta) followed by the character for the low 7 bits.
======================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================

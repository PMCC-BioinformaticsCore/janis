:orphan:

Concatenate
=================

``cat`` · *0 contributors · 1 version*

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
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
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

==============================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
name                            type                   prefix      position  documentation
==============================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
file                            Optional<File>
files                           Optional<Array<File>>                     1
number_output                   Optional<Boolean>      -n                    Number the output lines, starting at 1.
number_non_blank                Optional<Boolean>      -b                    Number the non-blank output lines, starting at 1.
disable_output_buffer           Optional<Boolean>      -u                    Disable output buffering.
squeeze                         Optional<Boolean>      -s                    Squeeze multiple adjacent empty lines, causing the output to be single spaced.
display_nonprint_and_eol_chars  Optional<Boolean>      -e                    Display non-printing characters (see the -v option), and display a dollar sign (`$') at the end of each line.
display_nonprint_and_tab_chars  Optional<Boolean>      -t                    Display non-printing characters (see the -v option), and display tab characters as `^I'.
display_nonprint_chars          Optional<Boolean>      -v                    Display non-printing characters so they are visible.  Control characters print as `^X' for control-X; the delete character (octal 0177) prints as `^?'.  Non-ASCII characters (with the high bit set) are printed as `M-' (for meta) followed by the character for the low 7 bits.
==============================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task cat {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? file
       Array[File]? files
       Boolean? number_output
       Boolean? number_non_blank
       Boolean? disable_output_buffer
       Boolean? squeeze
       Boolean? display_nonprint_and_eol_chars
       Boolean? display_nonprint_and_tab_chars
       Boolean? display_nonprint_chars
     }
     command <<<
       set -e
       cat \
         ~{if (defined(number_output) && select_first([number_output])) then "-n" else ""} \
         ~{if (defined(number_non_blank) && select_first([number_non_blank])) then "-b" else ""} \
         ~{if (defined(disable_output_buffer) && select_first([disable_output_buffer])) then "-u" else ""} \
         ~{if (defined(squeeze) && select_first([squeeze])) then "-s" else ""} \
         ~{if (defined(display_nonprint_and_eol_chars) && select_first([display_nonprint_and_eol_chars])) then "-e" else ""} \
         ~{if (defined(display_nonprint_and_tab_chars) && select_first([display_nonprint_and_tab_chars])) then "-t" else ""} \
         ~{if (defined(display_nonprint_chars) && select_first([display_nonprint_chars])) then "-v" else ""} \
         ~{if (defined(files) && length(select_first([files])) > 0) then "'" + sep("' '", select_first([files])) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Concatenate
   doc: |-
     The cat utility reads files sequentially, writing them to the standard output. The file operands are processed in command-line order. If file is a single dash (`-') or absent,cat reads from the standard input. If file is a UNIX domain socket, cat connects to it and then reads it until EOF. This complements the UNIX domain binding capability available in inetd(8).

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: file
     label: file
     type:
     - File
     - 'null'
   - id: files
     label: files
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       position: 1
   - id: number_output
     label: number_output
     doc: Number the output lines, starting at 1.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -n
   - id: number_non_blank
     label: number_non_blank
     doc: Number the non-blank output lines, starting at 1.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -b
   - id: disable_output_buffer
     label: disable_output_buffer
     doc: Disable output buffering.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -u
   - id: squeeze
     label: squeeze
     doc: Squeeze multiple adjacent empty lines, causing the output to be single spaced.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -s
   - id: display_nonprint_and_eol_chars
     label: display_nonprint_and_eol_chars
     doc: |-
       Display non-printing characters (see the -v option), and display a dollar sign (`$') at the end of each line.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -e
   - id: display_nonprint_and_tab_chars
     label: display_nonprint_and_tab_chars
     doc: |-
       Display non-printing characters (see the -v option), and display tab characters as `^I'.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -t
   - id: display_nonprint_chars
     label: display_nonprint_chars
     doc: |-
       Display non-printing characters so they are visible.  Control characters print as `^X' for control-X; the delete character (octal 0177) prints as `^?'.  Non-ASCII characters (with the high bit set) are printed as `M-' (for meta) followed by the character for the low 7 bits.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -v

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: cat
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: cat



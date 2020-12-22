:orphan:

Echo
===========

``echo`` · *0 contributors · 1 version*

The echo utility writes any specified operands, separated by single blank (` ') characters and followed by a newline (`
') character, to the standard output.

Some shells may provide a builtin echo command which is similar or identical to this utility. Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page.


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.echo import Echo

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "echo_step",
           Echo(
               inp=None,
           )
       )
       wf.output("out", source=echo_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for echo:

.. code-block:: bash

   # user inputs
   janis inputs echo > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp: <value>




5. Run echo with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       echo





Information
------------

:ID: ``echo``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
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

===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
name             type               prefix      position  documentation
===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
inp              String                                1
include_newline  Optional<Boolean>  -n                    Do not print the trailing newline character.  This may also be achieved by appending `\c' to the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) to suppress the newline character.
===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task echo {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String inp
       Boolean? include_newline
     }
     command <<<
       set -e
       echo \
         ~{if (defined(include_newline) && select_first([include_newline])) then "-n" else ""} \
         '~{inp}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
       duration: select_first([runtime_seconds, 60, 86400])
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
   label: Echo
   doc: |-
     The echo utility writes any specified operands, separated by single blank (` ') characters and followed by a newline (`
     ') character, to the standard output.

     Some shells may provide a builtin echo command which is similar or identical to this utility. Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: inp
     label: inp
     type: string
     inputBinding:
       position: 1
   - id: include_newline
     label: include_newline
     doc: |-
       Do not print the trailing newline character.  This may also be achieved by appending `\c' to the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) to suppress the newline character.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -n

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: echo
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 60, 86400].filter(function (inner) { return inner != null })[0])
   id: echo



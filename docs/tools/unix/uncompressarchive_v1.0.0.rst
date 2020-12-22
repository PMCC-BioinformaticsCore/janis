:orphan:

UncompressArchive
=================

``UncompressArchive`` · *0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.uncompressarchive import UncompressArchive

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "uncompressarchive_step",
           UncompressArchive(
               file=None,
           )
       )
       wf.output("out", source=uncompressarchive_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for UncompressArchive:

.. code-block:: bash

   # user inputs
   janis inputs UncompressArchive > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file.gz




5. Run UncompressArchive with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       UncompressArchive





Information
------------

:ID: ``UncompressArchive``
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

==========  =================  ===========  ==========  =======================================================
name        type               prefix         position  documentation
==========  =================  ===========  ==========  =======================================================
file        Gzipped<File>                            1
stdout      Optional<Boolean>  -c                       write on standard output, keep original files unchanged
decompress  Optional<Boolean>  -d                       decompress
force       Optional<Boolean>  -f                       force overwrite of output file and compress links
keep        Optional<Boolean>  -k                       keep (don't delete) input files
list        Optional<Boolean>  -l                       list compressed file contents
noName      Optional<Boolean>  -n                       do not save or restore the original name and time stamp
name        Optional<Boolean>  -N                       save or restore the original name and time stamp
quiet       Optional<Boolean>  -q                       suppress all warnings
recursive   Optional<Boolean>  -r                       operate recursively on directories
suffix      Optional<String>   -s                       use suffix SUF on compressed files
test        Optional<Boolean>  -t                       test compressed file integrity
fast        Optional<Boolean>  -1                       compress faster
best        Optional<Boolean>  -9                       compress better
rsyncable   Optional<Boolean>  --rsyncable              Make rsync-friendly archive
==========  =================  ===========  ==========  =======================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task UncompressArchive {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File file
       Boolean? stdout
       Boolean? decompress
       Boolean? force
       Boolean? keep
       Boolean? list
       Boolean? noName
       Boolean? name
       Boolean? quiet
       Boolean? recursive
       String? suffix
       Boolean? test
       Boolean? fast
       Boolean? best
       Boolean? rsyncable
     }
     command <<<
       set -e
       gunzip \
         ~{if select_first([stdout, true]) then "-c" else ""} \
         ~{if (defined(decompress) && select_first([decompress])) then "-d" else ""} \
         ~{if (defined(force) && select_first([force])) then "-f" else ""} \
         ~{if (defined(keep) && select_first([keep])) then "-k" else ""} \
         ~{if (defined(list) && select_first([list])) then "-l" else ""} \
         ~{if (defined(noName) && select_first([noName])) then "-n" else ""} \
         ~{if (defined(name) && select_first([name])) then "-N" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "-q" else ""} \
         ~{if (defined(recursive) && select_first([recursive])) then "-r" else ""} \
         ~{if defined(suffix) then ("-s '" + suffix + "'") else ""} \
         ~{if (defined(test) && select_first([test])) then "-t" else ""} \
         ~{if (defined(fast) && select_first([fast])) then "-1" else ""} \
         ~{if (defined(best) && select_first([best])) then "-9" else ""} \
         ~{if (defined(rsyncable) && select_first([rsyncable])) then "--rsyncable" else ""} \
         '~{file}'
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
   label: UncompressArchive

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: file
     label: file
     type: File
     inputBinding:
       position: 1
   - id: stdout
     label: stdout
     doc: write on standard output, keep original files unchanged
     type: boolean
     default: true
     inputBinding:
       prefix: -c
   - id: decompress
     label: decompress
     doc: decompress
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -d
   - id: force
     label: force
     doc: force overwrite of output file and compress links
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -f
   - id: keep
     label: keep
     doc: keep (don't delete) input files
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -k
   - id: list
     label: list
     doc: list compressed file contents
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -l
   - id: noName
     label: noName
     doc: do not save or restore the original name and time stamp
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -n
   - id: name
     label: name
     doc: save or restore the original name and time stamp
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -N
   - id: quiet
     label: quiet
     doc: suppress all warnings
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -q
   - id: recursive
     label: recursive
     doc: operate recursively on directories
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -r
   - id: suffix
     label: suffix
     doc: use suffix SUF on compressed files
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -s
   - id: test
     label: test
     doc: test compressed file integrity
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -t
   - id: fast
     label: fast
     doc: compress faster
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-1'
   - id: best
     label: best
     doc: compress better
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-9'
   - id: rsyncable
     label: rsyncable
     doc: Make rsync-friendly archive
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rsyncable

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: gunzip
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: UncompressArchive



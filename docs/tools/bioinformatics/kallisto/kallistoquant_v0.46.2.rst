:orphan:

Kallisto-Quant
==============================

``kallistoQuant`` · *1 contributor · 1 version*

Builds a kallisto index


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.kallisto.quant.versions import KallistoQuant_0_46_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "kallistoquant_step",
           KallistoQuant_0_46_2(
               index=None,
               fastq=None,
           )
       )
       wf.output("out", source=kallistoquant_step.out)
       wf.output("stats", source=kallistoquant_step.stats)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for kallistoQuant:

.. code-block:: bash

   # user inputs
   janis inputs kallistoQuant > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq
       - fastq_1.fastq
       index: index.kidx




5. Run kallistoQuant with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       kallistoQuant





Information
------------

:ID: ``kallistoQuant``
:URL: `https://pachterlab.github.io/kallisto/manual.html <https://pachterlab.github.io/kallisto/manual.html>`_
:Versions: v0.46.2
:Container: quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1
:Authors: Thomas Conway
:Citations: NL Bray, H Pimentel, P Melsted and L Pachter, Near optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, p 525--527 (2016).
:DOI: https://doi.org/10.1038/nbt.3519
:Created: 2020-05-25
:Updated: 2020-05-25


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
stats   File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

===============  ==================  =================  ==========  ========================================================================================
name             type                prefix               position  documentation
===============  ==================  =================  ==========  ========================================================================================
index            KallistoIdx         -i                          2  Filename for the kallisto index to be constructed
fastq            Array<Fastq>                                    4  FASTQ files to process
outdir           Optional<Filename>  -o                          3  directory to put outputs in
bias             Optional<Boolean>   --bias                         Perform sequence based bias correction
fusion           Optional<Boolean>   --fusion                       Search for fusions for Pizzly
single           Optional<Boolean>   --single                       Quantify single-end reads
overhang         Optional<Boolean>   --single-overhang              Include reads where unobserved rest of fragment is predicted to lie outside a transcript
fr_stranded      Optional<Boolean>   --fr-stranded                  Strand specific reads, first read forward
rf_stranded      Optional<Boolean>   --rf-stranded                  Strand specific reads, first read reverse
fragment_length  Optional<Double>    -l                             Estimated average fragment length
fragment_sd      Optional<Double>    -s                             Estimated standard deviation of fragment length
===============  ==================  =================  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task kallistoQuant {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File index
       String? outdir
       Array[File] fastq
       Boolean? bias
       Boolean? fusion
       Boolean? single
       Boolean? overhang
       Boolean? fr_stranded
       Boolean? rf_stranded
       Float? fragment_length
       Float? fragment_sd
     }
     command <<<
       set -e
       kallisto quant \
         ~{if (defined(bias) && select_first([bias])) then "--bias" else ""} \
         ~{if (defined(fusion) && select_first([fusion])) then "--fusion" else ""} \
         ~{if (defined(single) && select_first([single])) then "--single" else ""} \
         ~{if (defined(overhang) && select_first([overhang])) then "--single-overhang" else ""} \
         ~{if (defined(fr_stranded) && select_first([fr_stranded])) then "--fr-stranded" else ""} \
         ~{if (defined(rf_stranded) && select_first([rf_stranded])) then "--rf-stranded" else ""} \
         ~{if defined(fragment_length) then ("-l " + fragment_length) else ''} \
         ~{if defined(fragment_sd) then ("-s " + fragment_sd) else ''} \
         -i '~{index}' \
         -o '~{select_first([outdir, "generated"])}' \
         ~{if length(fastq) > 0 then "'" + sep("' '", fastq) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 2, 4])}G"
       preemptible: 2
     }
     output {
       File out = (select_first([outdir, "generated"]) + "/abundance.tsv")
       File stats = (select_first([outdir, "generated"]) + "/run_info.json")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Kallisto-Quant
   doc: Builds a kallisto index

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1

   inputs:
   - id: index
     label: index
     doc: Filename for the kallisto index to be constructed
     type: File
     inputBinding:
       prefix: -i
       position: 2
   - id: outdir
     label: outdir
     doc: directory to put outputs in
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: -o
       position: 3
   - id: fastq
     label: fastq
     doc: FASTQ files to process
     type:
       type: array
       items: File
     inputBinding:
       position: 4
   - id: bias
     label: bias
     doc: Perform sequence based bias correction
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --bias
   - id: fusion
     label: fusion
     doc: Search for fusions for Pizzly
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fusion
   - id: single
     label: single
     doc: Quantify single-end reads
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --single
   - id: overhang
     label: overhang
     doc: |-
       Include reads where unobserved rest of fragment is predicted to lie outside a transcript
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --single-overhang
   - id: fr_stranded
     label: fr_stranded
     doc: Strand specific reads, first read forward
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fr-stranded
   - id: rf_stranded
     label: rf_stranded
     doc: Strand specific reads, first read reverse
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rf-stranded
   - id: fragment_length
     label: fragment_length
     doc: Estimated average fragment length
     type:
     - double
     - 'null'
     inputBinding:
       prefix: -l
   - id: fragment_sd
     label: fragment_sd
     doc: Estimated standard deviation of fragment length
     type:
     - double
     - 'null'
     inputBinding:
       prefix: -s

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $((inputs.outdir + "/abundance.tsv"))
       outputEval: $((inputs.outdir.basename + "/abundance.tsv"))
       loadContents: false
   - id: stats
     label: stats
     type: File
     outputBinding:
       glob: $((inputs.outdir + "/run_info.json"))
       outputEval: $((inputs.outdir.basename + "/run_info.json"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - kallisto
   - quant
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: kallistoQuant



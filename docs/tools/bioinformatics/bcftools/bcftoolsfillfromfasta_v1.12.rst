:orphan:

BCFTools: fill-from-fasta
=================================================

``bcftoolsFillFromFasta`` · *1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.fill_from_fasta.versions import BcfToolsFillFromFasta_1_12

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolsfillfromfasta_step",
           BcfToolsFillFromFasta_1_12(
               vcf=None,
               column=None,
               fasta=None,
           )
       )
       wf.output("out", source=bcftoolsfillfromfasta_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for bcftoolsFillFromFasta:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolsFillFromFasta > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       column: <value>
       fasta: fasta.fasta
       vcf: null




5. Run bcftoolsFillFromFasta with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolsFillFromFasta

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          bcftoolsFillFromFasta

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``bcftoolsFillFromFasta``
:URL: *No URL to the documentation was provided*
:Versions: v1.12
:Container: quay.io/biocontainers/bcftools:1.12--h45bccc9_1
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-05-19
:Updated: 2021-05-19


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=================  ========================  ===================  ==========  ===============================================
name               type                      prefix                 position  documentation
=================  ========================  ===================  ==========  ===============================================
vcf                Union<VCF, Gzipped<VCF>>                                1  Input vcf
column             String                    --column                      3  REF or INFO tag, e.g. AA for ancestral allele
fasta              Fasta                     --fasta                       3  fasta file
outputFilename     Optional<Filename>                                      6  Output vcf
header_lines       Optional<File>            --header-lines                3  optional file containing header lines to append
include            Optional<String>          --include                     3  annotate only records passing filter expression
exclude            Optional<String>          --exclude                     3  annotate only records failing filter expression
replace_non_ACGTN  Optional<Boolean>         --replace-non-ACGTN           3  replace non-ACGTN characters with N
=================  ========================  ===================  ==========  ===============================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcftoolsFillFromFasta {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File vcf
       String? outputFilename
       String column
       File fasta
       File? header_lines
       String? include
       String? exclude
       Boolean? replace_non_ACGTN
     }

     command <<<
       set -e
       bcftools +fill-from-fasta \
         ~{vcf} \
         -- \
         --column '~{column}' \
         --fasta '~{fasta}' \
         ~{if defined(header_lines) then ("--header-lines '" + header_lines + "'") else ""} \
         ~{if defined(include) then ("--include '" + include + "'") else ""} \
         ~{if defined(exclude) then ("--exclude '" + exclude + "'") else ""} \
         ~{if (defined(replace_non_ACGTN) && select_first([replace_non_ACGTN])) then "--replace-non-ACGTN" else ""} \
         > \
         '~{select_first([outputFilename, "~{basename(vcf)}.fill.vcf"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "quay.io/biocontainers/bcftools:1.12--h45bccc9_1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "~{basename(vcf)}.fill.vcf"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BCFTools: fill-from-fasta'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/bcftools:1.12--h45bccc9_1

   inputs:
   - id: vcf
     label: vcf
     doc: Input vcf
     type: File
     inputBinding:
       position: 1
   - id: outputFilename
     label: outputFilename
     doc: Output vcf
     type:
     - string
     - 'null'
     default: generated.fill.vcf
     inputBinding:
       position: 6
       valueFrom: $(inputs.vcf.basename.replace(/.vcf$/, "").replace(/.vcf.gz$/, "")).fill.vcf
   - id: column
     label: column
     doc: REF or INFO tag, e.g. AA for ancestral allele
     type: string
     inputBinding:
       prefix: --column
       position: 3
   - id: fasta
     label: fasta
     doc: fasta file
     type: File
     inputBinding:
       prefix: --fasta
       position: 3
   - id: header_lines
     label: header_lines
     doc: optional file containing header lines to append
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --header-lines
       position: 3
   - id: include
     label: include
     doc: annotate only records passing filter expression
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --include
       position: 3
   - id: exclude
     label: exclude
     doc: annotate only records failing filter expression
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude
       position: 3
   - id: replace_non_ACGTN
     label: replace_non_ACGTN
     doc: replace non-ACGTN characters with N
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --replace-non-ACGTN
       position: 3

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.vcf.basename.replace(/.vcf$/, "").replace(/.vcf.gz$/, "")).fill.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bcftools
   - +fill-from-fasta
   arguments:
   - position: 2
     valueFrom: --
     shellQuote: false
   - position: 5
     valueFrom: '>'
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcftoolsFillFromFasta



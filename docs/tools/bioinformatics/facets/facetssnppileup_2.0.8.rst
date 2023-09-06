:orphan:

Facets: snp-pileup
====================================

``FacetsSnpPileup`` · *3 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.facets.snp_pileup.versions import FacetsSnpPileup_2_0_8

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "facetssnppileup_step",
           FacetsSnpPileup_2_0_8(
               vcf_file=None,
               normal_bam=None,
               tumor_bam=None,
           )
       )
       wf.output("out", source=facetssnppileup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for FacetsSnpPileup:

.. code-block:: bash

   # user inputs
   janis inputs FacetsSnpPileup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal_bam: normal_bam.bam
       tumor_bam: tumor_bam.bam
       vcf_file: vcf_file




5. Run FacetsSnpPileup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FacetsSnpPileup

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          FacetsSnpPileup

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``FacetsSnpPileup``
:URL: `https://github.com/mskcc/facets-suite <https://github.com/mskcc/facets-suite>`_
:Versions: 2.0.8
:Container: stevekm/facets-suite:2.0.8
:Authors: mumbler, evanwehi, Jiaan Yu
:Citations: Ronglai Shen, Venkatraman E. Seshan; FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing, Nucleic Acids Research, Volume 44, Issue 16, 19 September 2016, Pages e131,
:DOI: https://doi.org/10.1093/nar/gkw520
:Created: 2019-12-16
:Updated: 2021-03-04


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=============  ==================  ===============  ==========  ======================================================
name           type                prefix           position    documentation
=============  ==================  ===============  ==========  ======================================================
vcf_file       File                --vcf-file                   Path to VCF file containing SNP positions
normal_bam     IndexedBam          --normal-bam                 Path to normal sample BAM file
tumor_bam      IndexedBam          --tumor-bam                  Path to tumor sample BAM file
output_prefix  Optional<Filename>  --output-prefix              Path to VCF file containing SNP positions
pseudo_snps    Optional<Integer>   --pseudo-snps                Do pileup at every p:th position [default %(default)s]
max_depth      Optional<Integer>   --max-depth                  Maximum read depth [default %(default)s]
=============  ==================  ===============  ==========  ======================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task FacetsSnpPileup {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File vcf_file
       File normal_bam
       File normal_bam_bai
       File tumor_bam
       File tumor_bam_bai
       String? output_prefix
       Int? pseudo_snps
       Int? max_depth
     }

     command <<<
       set -e
        snp-pileup-wrapper.R \
         --vcf-file '~{vcf_file}' \
         --normal-bam '~{normal_bam}' \
         --tumor-bam '~{tumor_bam}' \
         --output-prefix '~{select_first([output_prefix, "generated"])}' \
         ~{if defined(pseudo_snps) then ("--pseudo-snps " + pseudo_snps) else ''} \
         ~{if defined(max_depth) then ("--max-depth " + max_depth) else ''}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "stevekm/facets-suite:2.0.8"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       File out = (select_first([output_prefix, "generated"]) + ".snp_pileup.gz")
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Facets: snp-pileup'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: stevekm/facets-suite:2.0.8

   inputs:
   - id: vcf_file
     label: vcf_file
     doc: Path to VCF file containing SNP positions
     type: File
     inputBinding:
       prefix: --vcf-file
   - id: normal_bam
     label: normal_bam
     doc: Path to normal sample BAM file
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --normal-bam
   - id: tumor_bam
     label: tumor_bam
     doc: Path to tumor sample BAM file
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --tumor-bam
   - id: output_prefix
     label: output_prefix
     doc: Path to VCF file containing SNP positions
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --output-prefix
   - id: pseudo_snps
     label: pseudo_snps
     doc: Do pileup at every p:th position [default %(default)s]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --pseudo-snps
   - id: max_depth
     label: max_depth
     doc: Maximum read depth [default %(default)s]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-depth

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $((inputs.output_prefix + ".snp_pileup.gz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ''
   - snp-pileup-wrapper.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: FacetsSnpPileup



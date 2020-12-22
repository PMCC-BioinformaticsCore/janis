:orphan:

Facets: snp-pileup
====================================

``FacetsSnpPileup`` · *2 contributors · 3 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.facets.snp_pileup.versions import FacetsSnpPileup_0_5_14_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "facetssnppileup_step",
           FacetsSnpPileup_0_5_14_1(
               vcf_file=None,
               normal=None,
               tumour=None,
           )
       )
       wf.output("out", source=facetssnppileup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FacetsSnpPileup:

.. code-block:: bash

   # user inputs
   janis inputs FacetsSnpPileup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal: normal.bam
       tumour: tumour.bam
       vcf_file: vcf_file.vcf




5. Run FacetsSnpPileup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FacetsSnpPileup





Information
------------

:ID: ``FacetsSnpPileup``
:URL: `https://github.com/vanallenlab/facets <https://github.com/vanallenlab/facets>`_
:Versions: 0.5.14.1, 0.5.14-2, 0.5.14
:Container: vanallenlab/facets:v0.5.14.1
:Authors: mumbler, evanwehi
:Citations: Ronglai Shen, Venkatraman E. Seshan; FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing, Nucleic Acids Research, Volume 44, Issue 16, 19 September 2016, Pages e131,
:DOI: https://doi.org/10.1093/nar/gkw520
:Created: 2019-12-16
:Updated: 2019-12-16


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

================  ==================  ===================  ==========  =======================================================================================================
name              type                prefix                 position  documentation
================  ==================  ===================  ==========  =======================================================================================================
vcf_file          VCF                                              18
normal            IndexedBam                                       20
tumour            IndexedBam                                       21
count_orphans     Optional<Boolean>   --count-orphans               2  Do not discard anomalous read pairs
ignore_overlaps   Optional<Boolean>   --ignore-overlaps             4  Disable read-pair overlap detection.
max_depth         Optional<Integer>   --maxdepth=                   6  Sets the maximum depth. Default is 4000.
min_map_quality   Optional<Integer>   --min-map-quality=            8  Sets the minimum threshold for mapping quality. Default is 0.
min_base_quality  Optional<Integer>   --min-base-quality=          10  Sets the minimum threshold for base quality. Default is 0.
min_read_counts   Optional<String>    --min-read-counts=           12  Comma separated list of minimum read counts for a position to be output. Default is 0.
gzip              Optional<Boolean>   --gzip                       14  Compresses the output file with BGZF.
pseudo_snps       Optional<String>    --pseudo-snps=               16  Every MULTIPLE positions, if there is no SNP,insert a blank record with the total count at theposition.
output_filename   Optional<Filename>                               19
================  ==================  ===================  ==========  =======================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task FacetsSnpPileup {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? count_orphans
       Boolean? ignore_overlaps
       Int? max_depth
       Int? min_map_quality
       Int? min_base_quality
       String? min_read_counts
       Boolean? gzip
       String? pseudo_snps
       File vcf_file
       String? output_filename
       File normal
       File normal_bai
       File tumour
       File tumour_bai
     }
     command <<<
       set -e
       LD_LIBRARY_PATH=/opt/conda/lib /snp-pileup \
         ~{if (defined(count_orphans) && select_first([count_orphans])) then "--count-orphans" else ""} \
         ~{if (defined(ignore_overlaps) && select_first([ignore_overlaps])) then "--ignore-overlaps" else ""} \
         ~{if defined(max_depth) then ("--maxdepth=" + max_depth) else ''} \
         ~{if defined(min_map_quality) then ("--min-map-quality=" + min_map_quality) else ''} \
         ~{if defined(min_base_quality) then ("--min-base-quality=" + min_base_quality) else ''} \
         ~{if defined(min_read_counts) then ("--min-read-counts='" + min_read_counts + "'") else ""} \
         ~{if (defined(gzip) && select_first([gzip])) then "--gzip" else ""} \
         ~{if defined(pseudo_snps) then ("--pseudo-snps='" + pseudo_snps + "'") else ""} \
         '~{vcf_file}' \
         '~{select_first([output_filename, "generated.csv.gz"])}' \
         '~{normal}' \
         '~{tumour}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "vanallenlab/facets:v0.5.14.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([output_filename, "generated.csv.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Facets: snp-pileup'
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: vanallenlab/facets:v0.5.14.1

   inputs:
   - id: count_orphans
     label: count_orphans
     doc: Do not discard anomalous read pairs
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --count-orphans
       position: 2
   - id: ignore_overlaps
     label: ignore_overlaps
     doc: Disable read-pair overlap detection.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-overlaps
       position: 4
   - id: max_depth
     label: max_depth
     doc: Sets the maximum depth. Default is 4000.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --maxdepth=
       position: 6
       separate: false
   - id: min_map_quality
     label: min_map_quality
     doc: Sets the minimum threshold for mapping quality. Default is 0.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-map-quality=
       position: 8
       separate: false
   - id: min_base_quality
     label: min_base_quality
     doc: Sets the minimum threshold for base quality. Default is 0.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-base-quality=
       position: 10
       separate: false
   - id: min_read_counts
     label: min_read_counts
     doc: |-
       Comma separated list of minimum read counts for a position to be output. Default is 0.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --min-read-counts=
       position: 12
       separate: false
   - id: gzip
     label: gzip
     doc: Compresses the output file with BGZF.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --gzip
       position: 14
   - id: pseudo_snps
     label: pseudo_snps
     doc: |-
       Every MULTIPLE positions, if there is no SNP,insert a blank record with the total count at theposition.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pseudo-snps=
       position: 16
       separate: false
   - id: vcf_file
     label: vcf_file
     type: File
     inputBinding:
       position: 18
   - id: output_filename
     label: output_filename
     type:
     - string
     - 'null'
     default: generated.csv.gz
     inputBinding:
       position: 19
   - id: normal
     label: normal
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       position: 20
   - id: tumour
     label: tumour
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       position: 21

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.csv.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - LD_LIBRARY_PATH=/opt/conda/lib /snp-pileup
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: FacetsSnpPileup



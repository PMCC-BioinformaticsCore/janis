:orphan:

Facets: Make plot
=============================

``RunFacets`` · *3 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.facets.run_facets.versions import RunFacets_2_0_8

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "runfacets_step",
           RunFacets_2_0_8(
               counts_file=None,
               directory=None,
               facets_lib_path=None,
           )
       )
       wf.output("out_summary", source=runfacets_step.out_summary)
       wf.output("out_purity_png", source=runfacets_step.out_purity_png)
       wf.output("out_purity_seg", source=runfacets_step.out_purity_seg)
       wf.output("out_purity_rds", source=runfacets_step.out_purity_rds)
       wf.output("out_hisens_png", source=runfacets_step.out_hisens_png)
       wf.output("out_hisens_seg", source=runfacets_step.out_hisens_seg)
       wf.output("out_hisens_rds", source=runfacets_step.out_hisens_rds)
       wf.output("out_arm_level", source=runfacets_step.out_arm_level)
       wf.output("out_gene_level", source=runfacets_step.out_gene_level)
       wf.output("out_qc", source=runfacets_step.out_qc)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for RunFacets:

.. code-block:: bash

   # user inputs
   janis inputs RunFacets > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       counts_file: counts_file




5. Run RunFacets with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       RunFacets

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          RunFacets

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``RunFacets``
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

==============  ==============  ===============
name            type            documentation
==============  ==============  ===============
out_summary     File
out_purity_png  File
out_purity_seg  File
out_purity_rds  File
out_hisens_png  File
out_hisens_seg  File
out_hisens_rds  File
out_arm_level   Optional<File>
out_gene_level  Optional<File>
out_qc          Optional<File>
==============  ==============  ===============


Additional configuration (inputs)
---------------------------------

===============  ==================  =================  ==========  =======================================================================================
name             type                prefix             position    documentation
===============  ==================  =================  ==========  =======================================================================================
counts_file      File                --counts-file                  Merged, gzipped tumor-normal output from snp-pileup
directory        String              --directory                    Output directory to which all output files are written
facets_lib_path  String              --facets-lib-path              path to the facets library. if none provided, uses version available to library(facets)
outputPrefix     Optional<Filename>  --sample-id                    Sample ID, preferrable Tumor_Normal to keep track of the normal used
everything       Optional<Boolean>   --everything                   Run full suite [default False]
genome           Optional<String>    --genome                       Reference genome [default hg19]
cval             Optional<Integer>   --cval                         Segmentation parameter (cval) [default 50]
purity_cval      Optional<Integer>   --purity-cval                  If two pass, purity segmentation parameter (cval)
min_nhet         Optional<Integer>   --min-nhet                     Min. number of heterozygous SNPs required for clustering [default 15]
purity_min_nhet  Optional<Integer>   --purity-min-nhet              If two pass, purity min. number of heterozygous SNPs (cval) [default 15]
snp_window_size  Optional<Integer>   --snp-window-size              Window size for heterozygous SNPs [default 250]
normal_depth     Optional<Integer>   --normal-depth                 Min. depth in normal to keep SNPs [default 35]
dipLogR          Optional<Double>    --dipLogR                      Manual dipLogR
seed             Optional<Integer>   --seed                         Manual seed value [default 100]
legacy_output    Optional<Boolean>   --legacy-output                create legacy output files (.RData and .cncf.txt) [default False]
===============  ==================  =================  ==========  =======================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task RunFacets {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File counts_file
       String? outputPrefix
       String? directory
       Boolean? everything
       String? genome
       Int? cval
       Int? purity_cval
       Int? min_nhet
       Int? purity_min_nhet
       Int? snp_window_size
       Int? normal_depth
       Float? dipLogR
       Int? seed
       Boolean? legacy_output
       String? facets_lib_path
     }

     command <<<
       set -e
        run-facets-wrapper.R \
         --counts-file '~{counts_file}' \
         --sample-id '~{select_first([outputPrefix, "generated"])}' \
         --directory '~{select_first([directory, "."])}' \
         ~{if (defined(everything) && select_first([everything])) then "--everything" else ""} \
         ~{if defined(genome) then ("--genome '" + genome + "'") else ""} \
         ~{if defined(cval) then ("--cval " + cval) else ''} \
         ~{if defined(purity_cval) then ("--purity-cval " + purity_cval) else ''} \
         ~{if defined(min_nhet) then ("--min-nhet " + min_nhet) else ''} \
         ~{if defined(purity_min_nhet) then ("--purity-min-nhet " + purity_min_nhet) else ''} \
         ~{if defined(snp_window_size) then ("--snp-window-size " + snp_window_size) else ''} \
         ~{if defined(normal_depth) then ("--normal-depth " + normal_depth) else ''} \
         ~{if defined(dipLogR) then ("--dipLogR " + dipLogR) else ''} \
         ~{if defined(seed) then ("--seed " + seed) else ''} \
         ~{if (defined(legacy_output) && select_first([legacy_output])) then "--legacy-output" else ""} \
         --facets-lib-path '~{select_first([facets_lib_path, "/usr/local/lib/R/site-library"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "stevekm/facets-suite:2.0.8"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 64, 4])}G"
       preemptible: 2
     }

     output {
       File out_summary = (select_first([outputPrefix, "generated"]) + ".txt")
       File out_purity_png = (select_first([outputPrefix, "generated"]) + "_purity.png")
       File out_purity_seg = (select_first([outputPrefix, "generated"]) + "_purity.seg")
       File out_purity_rds = (select_first([outputPrefix, "generated"]) + "_purity.rds")
       File out_hisens_png = (select_first([outputPrefix, "generated"]) + "_hisens.png")
       File out_hisens_seg = (select_first([outputPrefix, "generated"]) + "_hisens.seg")
       File out_hisens_rds = (select_first([outputPrefix, "generated"]) + "_hisens.rds")
       File? out_arm_level = (select_first([outputPrefix, "generated"]) + ".arm_level.txt")
       File? out_gene_level = (select_first([outputPrefix, "generated"]) + ".gene_level.txt")
       File? out_qc = (select_first([outputPrefix, "generated"]) + ".qc.txt")
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Facets: Make plot'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: stevekm/facets-suite:2.0.8

   inputs:
   - id: counts_file
     label: counts_file
     doc: Merged, gzipped tumor-normal output from snp-pileup
     type: File
     inputBinding:
       prefix: --counts-file
   - id: outputPrefix
     label: outputPrefix
     doc: Sample ID, preferrable Tumor_Normal to keep track of the normal used
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --sample-id
   - id: directory
     label: directory
     doc: Output directory to which all output files are written
     type: string
     default: .
     inputBinding:
       prefix: --directory
   - id: everything
     label: everything
     doc: Run full suite [default False]
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --everything
   - id: genome
     label: genome
     doc: Reference genome [default hg19]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genome
   - id: cval
     label: cval
     doc: Segmentation parameter (cval) [default 50]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cval
   - id: purity_cval
     label: purity_cval
     doc: If two pass, purity segmentation parameter (cval)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --purity-cval
   - id: min_nhet
     label: min_nhet
     doc: Min. number of heterozygous SNPs required for clustering [default 15]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-nhet
   - id: purity_min_nhet
     label: purity_min_nhet
     doc: If two pass, purity min. number of heterozygous SNPs (cval) [default 15]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --purity-min-nhet
   - id: snp_window_size
     label: snp_window_size
     doc: Window size for heterozygous SNPs [default 250]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --snp-window-size
   - id: normal_depth
     label: normal_depth
     doc: Min. depth in normal to keep SNPs [default 35]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --normal-depth
   - id: dipLogR
     label: dipLogR
     doc: Manual dipLogR
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --dipLogR
   - id: seed
     label: seed
     doc: Manual seed value [default 100]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seed
   - id: legacy_output
     label: legacy_output
     doc: create legacy output files (.RData and .cncf.txt) [default False]
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --legacy-output
   - id: facets_lib_path
     label: facets_lib_path
     doc: |-
       path to the facets library. if none provided, uses version available to library(facets)
     type: string
     default: /usr/local/lib/R/site-library
     inputBinding:
       prefix: --facets-lib-path

   outputs:
   - id: out_summary
     label: out_summary
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + ".txt"))
       loadContents: false
   - id: out_purity_png
     label: out_purity_png
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_purity.png"))
       loadContents: false
   - id: out_purity_seg
     label: out_purity_seg
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_purity.seg"))
       loadContents: false
   - id: out_purity_rds
     label: out_purity_rds
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_purity.rds"))
       loadContents: false
   - id: out_hisens_png
     label: out_hisens_png
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_hisens.png"))
       loadContents: false
   - id: out_hisens_seg
     label: out_hisens_seg
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_hisens.seg"))
       loadContents: false
   - id: out_hisens_rds
     label: out_hisens_rds
     type: File
     outputBinding:
       glob: $((inputs.outputPrefix + "_hisens.rds"))
       loadContents: false
   - id: out_arm_level
     label: out_arm_level
     type:
     - File
     - 'null'
     outputBinding:
       glob: $((inputs.outputPrefix + ".arm_level.txt"))
       loadContents: false
   - id: out_gene_level
     label: out_gene_level
     type:
     - File
     - 'null'
     outputBinding:
       glob: $((inputs.outputPrefix + ".gene_level.txt"))
       loadContents: false
   - id: out_qc
     label: out_qc
     type:
     - File
     - 'null'
     outputBinding:
       glob: $((inputs.outputPrefix + ".qc.txt"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ''
   - run-facets-wrapper.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: RunFacets



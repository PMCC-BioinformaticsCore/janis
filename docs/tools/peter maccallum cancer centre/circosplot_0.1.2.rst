:orphan:

Circos Plot
========================

``CircosPlot`` · *1 contributor · 1 version*

Command: Rscript /app/circos_plot/circos_plot_facets_manta.R $tumor_name $normal_name $facets_file $sv_file $out_dir $manta_filter[optional]
    Output: PDF file: ${tumor_name}--${normal_name}.pdf



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.circosplot.versions import CircosPlot_0_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "circosplot_step",
           CircosPlot_0_1_2(
               tumor_name=None,
               normal_name=None,
               facets_file=None,
               sv_file=None,
               genome=None,
           )
       )
       wf.output("out", source=circosplot_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for CircosPlot:

.. code-block:: bash

   # user inputs
   janis inputs CircosPlot > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       facets_file: facets_file
       genome: <value>
       normal_name: <value>
       sv_file: sv_file.vcf.gz
       tumor_name: <value>




5. Run CircosPlot with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       CircosPlot

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          CircosPlot

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``CircosPlot``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.2
:Container: rlupat/pmacutil:latest
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-07-09
:Updated: 2021-08-27


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

============  =================  ========  ==========  ===============
name          type               prefix      position  documentation
============  =================  ========  ==========  ===============
tumor_name    String                                1
normal_name   String                                2
facets_file   File                                  3
sv_file       Gzipped<VCF>                          4
genome        String                                6
output_dir    Optional<String>                      5
manta_filter  Optional<Integer>                     7
============  =================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task CircosPlot {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       String tumor_name
       String normal_name
       File facets_file
       File sv_file
       String? output_dir
       String genome
       Int? manta_filter
     }

     command <<<
       set -e
       Rscript /app/circos_plot/circos_plot_facets_manta.R \
         '~{tumor_name}' \
         '~{normal_name}' \
         '~{facets_file}' \
         '~{sv_file}' \
         ~{if defined(select_first([output_dir, "."])) then ("'" + select_first([output_dir, "."]) + "'") else ""} \
         '~{genome}' \
         ~{manta_filter}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "rlupat/pmacutil:latest"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }

     output {
       File out = "~{select_first([output_dir, "."])}/~{tumor_name}--~{normal_name}.pdf"
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Circos Plot

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: rlupat/pmacutil:latest

   inputs:
   - id: tumor_name
     label: tumor_name
     type: string
     inputBinding:
       position: 1
   - id: normal_name
     label: normal_name
     type: string
     inputBinding:
       position: 2
   - id: facets_file
     label: facets_file
     type: File
     inputBinding:
       position: 3
   - id: sv_file
     label: sv_file
     type: File
     inputBinding:
       position: 4
   - id: output_dir
     label: output_dir
     type: string
     default: .
     inputBinding:
       position: 5
   - id: genome
     label: genome
     type: string
     inputBinding:
       position: 6
   - id: manta_filter
     label: manta_filter
     type:
     - int
     - 'null'
     inputBinding:
       position: 7

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: |-
         "{output_dir}/{tumor_name}--{normal_name}.pdf".replace(/\{output_dir\}/g, inputs.output_dir).replace(/\{tumor_name\}/g, inputs.tumor_name).replace(/\{normal_name\}/g, inputs.normal_name)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: Rscript /app/circos_plot/circos_plot_facets_manta.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: CircosPlot



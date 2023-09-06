:orphan:

Arriba: DrawFusions
=======================================

``ArribaDrawFusions`` · *1 contributor · 2 versions*


Arriba comes with an R script draw_fusions.R that renders publication-quality 
visualizations of the transcripts involved in predicted fusions. It generates 
a PDF file with one page for each predicted fusion. Each page depicts the 
fusion partners, their orientation, the retained exons in the fusion 
transcript, statistics about the number of supporting reads, and - if the 
column fusion_transcript has a value - an excerpt of the sequence around the 
breakpoint.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.arriba.draw_fusions.versions import ArribaDrawFusions_2_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "arribadrawfusions_step",
           ArribaDrawFusions_2_1_0(
               annotation=None,
               fusions=None,
           )
       )
       wf.output("out", source=arribadrawfusions_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for ArribaDrawFusions:

.. code-block:: bash

   # user inputs
   janis inputs ArribaDrawFusions > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       annotation: annotation
       fusions: fusions




5. Run ArribaDrawFusions with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       ArribaDrawFusions

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          ArribaDrawFusions

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``ArribaDrawFusions``
:URL: *No URL to the documentation was provided*
:Versions: 2.1.0, 1.1.0
:Container: quay.io/biocontainers/arriba:2.1.0--hd2e4403_0
:Authors: Jiaan Yu
:Citations: None
:Created: 2021-03-15
:Updated: 2021-03-15


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==========================  ====================  =============================  ==========  ==================
name                        type                  prefix                         position    documentation
==========================  ====================  =============================  ==========  ==================
annotation                  File                  --annotation=                              exonsFile
fusions                     File                  --fusions=                                 fusionsFile
outputFilename              Optional<Filename>    --output=                                  outputFile
alignments                  Optional<IndexedBam>  --alignments=                              alignmentsFile
cytobands                   Optional<File>        --cytobands=                               cytobandsFile
minConfidenceForCircosPlot  Optional<String>      --minConfidenceForCircosPlot=
proteinDomains              Optional<File>        --proteinDomains=                          proteinDomainsFile
squishIntrons               Optional<Boolean>     --squishIntrons=
printExonLabels             Optional<Boolean>     --printExonLabels=
renderThreeDEffect          Optional<Boolean>     --render3dEffect=
pdfWidth                    Optional<Float>       --pdfWidth=
pdfHeight                   Optional<Float>       --pdfHeight=
color_one                   Optional<String>      --color1=
color_two                   Optional<String>      --color2=
mergeDomainsOverlappingBy   Optional<Float>       --mergeDomainsOverlappingBy=
optimizeDomainColors        Optional<Boolean>     --optimizeDomainColors=
fontSize                    Optional<Integer>     --fontSize=
showIntergenicVicinity      Optional<Integer>     --showIntergenicVicinity=
transcriptSelection         Optional<String>      --transcriptSelection=
==========================  ====================  =============================  ==========  ==================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task ArribaDrawFusions {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File annotation
       File fusions
       String? outputFilename
       File? alignments
       File? alignments_bai
       File? cytobands
       String? minConfidenceForCircosPlot
       File? proteinDomains
       Boolean? squishIntrons
       Boolean? printExonLabels
       Boolean? renderThreeDEffect
       Float? pdfWidth
       Float? pdfHeight
       String? color_one
       String? color_two
       Float? mergeDomainsOverlappingBy
       Boolean? optimizeDomainColors
       Int? fontSize
       Int? showIntergenicVicinity
       String? transcriptSelection
     }

     command <<<
       set -e
        /usr/local/bin/draw_fusions.R \
         --annotation='~{annotation}' \
         --fusions='~{fusions}' \
         --output='~{select_first([outputFilename, "generated.pdf"])}' \
         ~{if defined(alignments) then ("--alignments='" + alignments + "'") else ""} \
         ~{if defined(cytobands) then ("--cytobands='" + cytobands + "'") else ""} \
         ~{if defined(minConfidenceForCircosPlot) then ("--minConfidenceForCircosPlot='" + minConfidenceForCircosPlot + "'") else ""} \
         ~{if defined(proteinDomains) then ("--proteinDomains='" + proteinDomains + "'") else ""} \
         ~{if (defined(squishIntrons) && select_first([squishIntrons])) then "--squishIntrons=" else ""} \
         ~{if (defined(printExonLabels) && select_first([printExonLabels])) then "--printExonLabels=" else ""} \
         ~{if (defined(renderThreeDEffect) && select_first([renderThreeDEffect])) then "--render3dEffect=" else ""} \
         ~{if defined(pdfWidth) then ("--pdfWidth=" + pdfWidth) else ''} \
         ~{if defined(pdfHeight) then ("--pdfHeight=" + pdfHeight) else ''} \
         ~{if defined(color_one) then ("--color1='" + color_one + "'") else ""} \
         ~{if defined(color_two) then ("--color2='" + color_two + "'") else ""} \
         ~{if defined(mergeDomainsOverlappingBy) then ("--mergeDomainsOverlappingBy=" + mergeDomainsOverlappingBy) else ''} \
         ~{if (defined(optimizeDomainColors) && select_first([optimizeDomainColors])) then "--optimizeDomainColors=" else ""} \
         ~{if defined(fontSize) then ("--fontSize=" + fontSize) else ''} \
         ~{if defined(showIntergenicVicinity) then ("--showIntergenicVicinity=" + showIntergenicVicinity) else ''} \
         ~{if defined(transcriptSelection) then ("--transcriptSelection='" + transcriptSelection + "'") else ""}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "quay.io/biocontainers/arriba:2.1.0--hd2e4403_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "generated.pdf"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Arriba: DrawFusions'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/arriba:2.1.0--hd2e4403_0

   inputs:
   - id: annotation
     label: annotation
     doc: exonsFile
     type: File
     inputBinding:
       prefix: --annotation=
       separate: false
   - id: fusions
     label: fusions
     doc: fusionsFile
     type: File
     inputBinding:
       prefix: --fusions=
       separate: false
   - id: outputFilename
     label: outputFilename
     doc: outputFile
     type:
     - string
     - 'null'
     default: generated.pdf
     inputBinding:
       prefix: --output=
       separate: false
   - id: alignments
     label: alignments
     doc: alignmentsFile
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: --alignments=
       separate: false
   - id: cytobands
     label: cytobands
     doc: cytobandsFile
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --cytobands=
       separate: false
   - id: minConfidenceForCircosPlot
     label: minConfidenceForCircosPlot
     doc: ''
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --minConfidenceForCircosPlot=
       separate: false
   - id: proteinDomains
     label: proteinDomains
     doc: proteinDomainsFile
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --proteinDomains=
       separate: false
   - id: squishIntrons
     label: squishIntrons
     doc: ''
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --squishIntrons=
       separate: false
   - id: printExonLabels
     label: printExonLabels
     doc: ''
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --printExonLabels=
       separate: false
   - id: renderThreeDEffect
     label: renderThreeDEffect
     doc: ''
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --render3dEffect=
       separate: false
   - id: pdfWidth
     label: pdfWidth
     doc: ''
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --pdfWidth=
       separate: false
   - id: pdfHeight
     label: pdfHeight
     doc: ''
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --pdfHeight=
       separate: false
   - id: color_one
     label: color_one
     doc: ''
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --color1=
       separate: false
   - id: color_two
     label: color_two
     doc: ''
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --color2=
       separate: false
   - id: mergeDomainsOverlappingBy
     label: mergeDomainsOverlappingBy
     doc: ''
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --mergeDomainsOverlappingBy=
       separate: false
   - id: optimizeDomainColors
     label: optimizeDomainColors
     doc: ''
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --optimizeDomainColors=
       separate: false
   - id: fontSize
     label: fontSize
     doc: ''
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --fontSize=
       separate: false
   - id: showIntergenicVicinity
     label: showIntergenicVicinity
     doc: ''
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --showIntergenicVicinity=
       separate: false
   - id: transcriptSelection
     label: transcriptSelection
     doc: ''
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --transcriptSelection=
       separate: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.pdf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ''
   - /usr/local/bin/draw_fusions.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: ArribaDrawFusions



:orphan:

GATK4: Print Reads
====================================

``Gatk4PrintReads`` · *1 contributor · 4 versions*


Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
A common use case is to subset reads by genomic interval using the -L argument. 
Note when applying genomic intervals, the tool is literal and does not retain mates 
of paired-end reads outside of the interval, if any. Data with missing mates will fail 
ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. 
If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

By default, PrintReads applies the WellformedReadFilter at the engine level. 
What this means is that the tool does not print reads that fail the WellformedReadFilter filter. 
You can similarly apply other engine-level filters to remove specific types of reads 
with the --read-filter argument. See documentation category 'Read Filters' for a list of
 available filters. To keep reads that do not pass the WellformedReadFilter, either 
 disable the filter with --disable-read-filter or disable all default filters with 
 ``--disable-tool-default-read-filters``.

The reference is strictly required when handling CRAM files.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.printreads.versions import Gatk4PrintReads_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4printreads_step",
           Gatk4PrintReads_4_1_2(
               bam=None,
           )
       )
       wf.output("out", source=gatk4printreads_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4PrintReads:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4PrintReads > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run Gatk4PrintReads with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4PrintReads





Information
------------

:ID: ``Gatk4PrintReads``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Additional configuration (inputs)
---------------------------------

=================  =======================  ========  ==========  ========================================================================================
name               type                     prefix    position    documentation
=================  =======================  ========  ==========  ========================================================================================
bam                BAM
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                              Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename     Optional<Filename>
=================  =======================  ========  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4PrintReads {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File bam
       String? outputFilename
     }
     command <<<
       set -e
       gatk PrintReads \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}'
       if [ -f $(echo '~{select_first([outputFilename, "generated"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "generated"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "generated"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated"])
       File out_bai = select_first([outputFilename, "generated"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Print Reads'
   doc: |2-

     Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
     A common use case is to subset reads by genomic interval using the -L argument. 
     Note when applying genomic intervals, the tool is literal and does not retain mates 
     of paired-end reads outside of the interval, if any. Data with missing mates will fail 
     ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. 
     If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

     By default, PrintReads applies the WellformedReadFilter at the engine level. 
     What this means is that the tool does not print reads that fail the WellformedReadFilter filter. 
     You can similarly apply other engine-level filters to remove specific types of reads 
     with the --read-filter argument. See documentation category 'Read Filters' for a list of
      available filters. To keep reads that do not pass the WellformedReadFilter, either 
      disable the filter with --disable-read-filter or disable all default filters with 
      ``--disable-tool-default-read-filters``.

     The reference is strictly required when handling CRAM files.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.2.0

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: bam
     label: bam
     type: File
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }
               return [
                       {
                           path: resolveSecondary(self.path, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     outputBinding:
       glob: generated
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - PrintReads
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4PrintReads



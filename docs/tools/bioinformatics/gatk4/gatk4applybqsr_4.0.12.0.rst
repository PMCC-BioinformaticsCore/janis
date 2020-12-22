:orphan:

GATK4: Apply base quality score recalibration
==============================================================

``Gatk4ApplyBQSR`` · *1 contributor · 4 versions*

Apply base quality score recalibration: This tool performs the second pass in a two-stage 
process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
base qualities of the input reads based on the recalibration table produced by the 
BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
that affect the assignment of base quality scores by the sequencer. The first pass consists 
of calculating error empirically and finding patterns in how error varies with basecall 
features over all bases. The relevant observations are written to a recalibration table. 
The second pass consists of applying numerical corrections to each individual basecall 
based on the patterns identified in the first step (recorded in the recalibration table) 
and write out the recalibrated data to a new BAM or CRAM file.

- This tool replaces the use of PrintReads for the application of base quality score 
    recalibration as practiced in earlier versions of GATK (2.x and 3.x).
- You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
- Original qualities can be retained in the output file under the "OQ" tag if desired. 
    See the `--emit-original-quals` argument for details.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.applybqsr.versions import Gatk4ApplyBqsr_4_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4applybqsr_step",
           Gatk4ApplyBqsr_4_0(
               bam=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4applybqsr_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4ApplyBQSR:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4ApplyBQSR > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       reference: reference.fasta




5. Run Gatk4ApplyBQSR with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4ApplyBQSR





Information
------------

:ID: ``Gatk4ApplyBQSR``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
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

=================  =======================  =================  ==========  ========================================================================================
name               type                     prefix               position  documentation
=================  =======================  =================  ==========  ========================================================================================
bam                IndexedBam               -I                         10  The SAM/BAM/CRAM file containing reads.
reference          FastaWithIndexes         -R                             Reference sequence
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                                       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputFilename     Optional<Filename>       -O                             Write output to this file
recalFile          Optional<tsv>            --bqsr-recal-file              Input recalibration table for BQSR
intervals          Optional<bed>            --intervals                    -L (BASE) One or more genomic intervals over which to operate
intervalStrings    Optional<Array<String>>  --intervals                    -L (BASE) One or more genomic intervals over which to operate
tmpDir             Optional<String>         --tmp-dir                  11  Temp directory to use.
=================  =======================  =================  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4ApplyBQSR {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File bam
       File bam_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       File? recalFile
       File? intervals
       Array[String]? intervalStrings
       String? tmpDir
     }
     command <<<
       set -e
       cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
       gatk ApplyBQSR \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         -R '~{reference}' \
         -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' \
         ~{if defined(recalFile) then ("--bqsr-recal-file '" + recalFile + "'") else ""} \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(intervalStrings) && length(select_first([intervalStrings])) > 0) then "--intervals '" + sep("' --intervals '", select_first([intervalStrings])) + "'" else ""} \
         -I '~{bam}' \
         ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--tmp-dir '" + select_first([tmpDir, "/tmp/"]) + "'") else ""}
       if [ -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.0.12.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"])
       File out_bai = select_first([outputFilename, "~{basename(bam, ".bam")}.recalibrated.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Apply base quality score recalibration'
   doc: |-
     Apply base quality score recalibration: This tool performs the second pass in a two-stage 
     process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
     base qualities of the input reads based on the recalibration table produced by the 
     BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

     Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
     that affect the assignment of base quality scores by the sequencer. The first pass consists 
     of calculating error empirically and finding patterns in how error varies with basecall 
     features over all bases. The relevant observations are written to a recalibration table. 
     The second pass consists of applying numerical corrections to each individual basecall 
     based on the patterns identified in the first step (recorded in the recalibration table) 
     and write out the recalibrated data to a new BAM or CRAM file.

     - This tool replaces the use of PrintReads for the application of base quality score 
         recalibration as practiced in earlier versions of GATK (2.x and 3.x).
     - You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
     - Original qualities can be retained in the output file under the "OQ" tag if desired. 
         See the `--emit-original-quals` argument for details.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.0.12.0

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
     doc: The SAM/BAM/CRAM file containing reads.
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
                           location: resolveSecondary(self.location, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     inputBinding:
       prefix: -I
       position: 10
   - id: reference
     label: reference
     doc: Reference sequence
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: -R
   - id: outputFilename
     label: outputFilename
     doc: Write output to this file
     type:
     - string
     - 'null'
     default: generated.recalibrated.bam
     inputBinding:
       prefix: -O
       valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).recalibrated.bam
   - id: recalFile
     label: recalFile
     doc: Input recalibration table for BQSR
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --bqsr-recal-file
   - id: intervals
     label: intervals
     doc: -L (BASE) One or more genomic intervals over which to operate
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --intervals
   - id: intervalStrings
     label: intervalStrings
     doc: -L (BASE) One or more genomic intervals over which to operate
     type:
     - type: array
       inputBinding:
         prefix: --intervals
       items: string
     - 'null'
     inputBinding: {}
   - id: tmpDir
     label: tmpDir
     doc: Temp directory to use.
     type: string
     default: /tmp/
     inputBinding:
       prefix: --tmp-dir
       position: 11

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
       glob: $(inputs.bam.basename.replace(/.bam$/, "")).recalibrated.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - ApplyBQSR
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4ApplyBQSR



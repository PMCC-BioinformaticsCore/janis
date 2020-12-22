:orphan:

GATK4: Base Recalibrator
================================================

``Gatk4BaseRecalibrator`` · *1 contributor · 4 versions*

First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
Since there is a large amount of data one can then calculate an empirical probability of error given the 
particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
table (of the several covariate values, num observations, num mismatches, empirical quality score).


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.baserecalibrator.versions import Gatk4BaseRecalibrator_4_1_4

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4baserecalibrator_step",
           Gatk4BaseRecalibrator_4_1_4(
               bam=None,
               knownSites=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4baserecalibrator_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4BaseRecalibrator:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4BaseRecalibrator > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       knownSites:
       - knownSites_0.vcf.gz
       - knownSites_1.vcf.gz
       reference: reference.fasta




5. Run Gatk4BaseRecalibrator with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4BaseRecalibrator





Information
------------

:ID: ``Gatk4BaseRecalibrator``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=================  =======================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name               type                     prefix           position  documentation
=================  =======================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bam                IndexedBam               -I                      6  BAM/SAM/CRAM file containing reads
knownSites         Array<Gzipped<VCF>>      --known-sites          28  **One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.** This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference, so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
reference          FastaWithIndexes         -R                      5  Reference sequence file
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                                   Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
tmpDir             Optional<String>         --tmp-dir                  Temp directory to use.
outputFilename     Optional<Filename>       -O                      8  **The output recalibration table filename to create.** After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate. Use '/dev/stdout' to print to standard out.
intervals          Optional<bed>            --intervals                -L (BASE) One or more genomic intervals over which to operate
intervalStrings    Optional<Array<String>>  --intervals                -L (BASE) One or more genomic intervals over which to operate
=================  =======================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4BaseRecalibrator {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       String? tmpDir
       File bam
       File bam_bai
       Array[File] knownSites
       Array[File] knownSites_tbi
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       File? intervals
       Array[String]? intervalStrings
     }
     command <<<
       set -e
       cp -f '~{bam_bai}' $(echo '~{bam}' | sed 's/\.[^.]*$//').bai
       gatk BaseRecalibrator \
         --java-options '-Xmx~{((select_first([runtime_memory, 16, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--tmp-dir '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         ~{if (defined(intervalStrings) && length(select_first([intervalStrings])) > 0) then "--intervals '" + sep("' --intervals '", select_first([intervalStrings])) + "'" else ""} \
         -R '~{reference}' \
         -I '~{bam}' \
         -O '~{select_first([outputFilename, "~{basename(bam, ".bam")}.table"])}' \
         ~{if length(knownSites) > 0 then "--known-sites '" + sep("' --known-sites '", knownSites) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.4.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(bam, ".bam")}.table"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Base Recalibrator'
   doc: |-
     First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
     The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

     This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
     that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
     We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
     Since there is a large amount of data one can then calculate an empirical probability of error given the 
     particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
     table (of the several covariate values, num observations, num mismatches, empirical quality score).

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.4.0

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
   - id: tmpDir
     label: tmpDir
     doc: Temp directory to use.
     type: string
     default: /tmp/
     inputBinding:
       prefix: --tmp-dir
   - id: bam
     label: bam
     doc: BAM/SAM/CRAM file containing reads
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
       position: 6
   - id: knownSites
     label: knownSites
     doc: |-
       **One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.** This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference, so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
     type:
       type: array
       inputBinding:
         prefix: --known-sites
       items: File
     inputBinding:
       position: 28
   - id: reference
     label: reference
     doc: Reference sequence file
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
       position: 5
   - id: outputFilename
     label: outputFilename
     doc: |-
       **The output recalibration table filename to create.** After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate. Use '/dev/stdout' to print to standard out.
     type:
     - string
     - 'null'
     default: generated.table
     inputBinding:
       prefix: -O
       position: 8
       valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).table
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

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.bam.basename.replace(/.bam$/, "")).table
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - BaseRecalibrator
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 16, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4BaseRecalibrator



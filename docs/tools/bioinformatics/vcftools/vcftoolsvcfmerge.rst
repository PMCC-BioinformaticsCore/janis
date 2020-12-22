:orphan:

VcfTools: VcfMerge
=====================================

``VcfToolsVcfMerge`` · *1 contributor · 1 version*

Merges two or more VCF files into one so that, for example, if two source files had one column each, on output will be printed a file with two columns. See also vcf-concat for concatenating VCFs split by chromosome.

vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

Note that this script is not intended for concatenating VCF files. For this, use vcf-concat instead.
Note: A fast htslib C version of this tool is now available (see bcftools merge).


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcftools.vcfmerge.versions import VcfToolsVcfMerge_0_1_16

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcftoolsvcfmerge_step",
           VcfToolsVcfMerge_0_1_16(
               vcfTabix=None,
           )
       )
       wf.output("out", source=vcftoolsvcfmerge_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for VcfToolsVcfMerge:

.. code-block:: bash

   # user inputs
   janis inputs VcfToolsVcfMerge > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcfTabix:
       - vcfTabix_0.vcf.gz
       - vcfTabix_1.vcf.gz




5. Run VcfToolsVcfMerge with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       VcfToolsVcfMerge





Information
------------

:ID: ``VcfToolsVcfMerge``
:URL: `http://vcftools.sourceforge.net/perl_module.html#vcf-merge <http://vcftools.sourceforge.net/perl_module.html#vcf-merge>`_
:Versions: 0.1.16
:Container: biocontainers/vcftools:v0.1.16-1-deb_cv1
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-05-21
:Updated: 2020-05-21


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

================  =======================  ===================  ==========  ===================================================================================================================================================================
name              type                     prefix                 position  documentation
================  =======================  ===================  ==========  ===================================================================================================================================================================
vcfTabix          Array<Gzipped<VCF>>                                   10
collapse          Optional<String>         -c                               treat as identical sites with differing alleles [any] <snps|indels|both|any|none>
removeDuplicates  Optional<Boolean>        --remove-duplicates              If there should be two consecutive rows with the same chr:pos, print only the first one.
vcfHeader         Optional<File>           --vcf-header                     Use the provided VCF header
regionsList       Optional<Array<String>>  --regions                        Do only the given regions (comma-separated list).
regionsFile       Optional<File>           --regions                        Do only the given regions (one region per line in a file).
refForMissing     Optional<String>         --ref-for-missing                Use the REF allele instead of the default missing genotype. Because it is not obvious what ploidy should be used, a user-defined string is used instead (e.g. 0/0).
silent            Optional<Boolean>        --silent                         Try to be a bit more silent, no warnings about duplicate lines.
trimALTs          Optional<Boolean>        --trim-ALTs                      If set, redundant ALTs will be removed
================  =======================  ===================  ==========  ===================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task VcfToolsVcfMerge {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? collapse
       Boolean? removeDuplicates
       File? vcfHeader
       Array[String]? regionsList
       File? regionsFile
       String? refForMissing
       Boolean? silent
       Boolean? trimALTs
       Array[File] vcfTabix
       Array[File] vcfTabix_tbi
     }
     command <<<
       set -e
        vcf-merge \
         ~{if defined(collapse) then ("-c '" + collapse + "'") else ""} \
         ~{if (defined(removeDuplicates) && select_first([removeDuplicates])) then "--remove-duplicates" else ""} \
         ~{if defined(vcfHeader) then ("--vcf-header '" + vcfHeader + "'") else ""} \
         ~{if (defined(regionsList) && length(select_first([regionsList])) > 0) then "--regions '" + sep("','", select_first([regionsList])) + "'" else ""} \
         ~{if defined(regionsFile) then ("--regions '" + regionsFile + "'") else ""} \
         ~{if defined(refForMissing) then ("--ref-for-missing '" + refForMissing + "'") else ""} \
         ~{if (defined(silent) && select_first([silent])) then "--silent" else ""} \
         ~{if (defined(trimALTs) && select_first([trimALTs])) then "--trim-ALTs" else ""} \
         ~{if length(vcfTabix) > 0 then "'" + sep("' '", vcfTabix) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/vcftools:v0.1.16-1-deb_cv1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'VcfTools: VcfMerge'
   doc: |-
     Merges two or more VCF files into one so that, for example, if two source files had one column each, on output will be printed a file with two columns. See also vcf-concat for concatenating VCFs split by chromosome.

     vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

     Note that this script is not intended for concatenating VCF files. For this, use vcf-concat instead.
     Note: A fast htslib C version of this tool is now available (see bcftools merge).

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/vcftools:v0.1.16-1-deb_cv1

   inputs:
   - id: collapse
     label: collapse
     doc: |-
       treat as identical sites with differing alleles [any] <snps|indels|both|any|none> 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -c
   - id: removeDuplicates
     label: removeDuplicates
     doc: |-
       If there should be two consecutive rows with the same chr:pos, print only the first one.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --remove-duplicates
   - id: vcfHeader
     label: vcfHeader
     doc: Use the provided VCF header
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --vcf-header
   - id: regionsList
     label: regionsList
     doc: Do only the given regions (comma-separated list).
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --regions
       itemSeparator: ','
   - id: regionsFile
     label: regionsFile
     doc: Do only the given regions (one region per line in a file).
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --regions
   - id: refForMissing
     label: refForMissing
     doc: |-
       Use the REF allele instead of the default missing genotype. Because it is not obvious what ploidy should be used, a user-defined string is used instead (e.g. 0/0).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --ref-for-missing
   - id: silent
     label: silent
     doc: Try to be a bit more silent, no warnings about duplicate lines.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --silent
   - id: trimALTs
     label: trimALTs
     doc: If set, redundant ALTs will be removed
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --trim-ALTs
   - id: vcfTabix
     label: vcfTabix
     type:
       type: array
       items: File
     inputBinding:
       position: 10

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ''
   - vcf-merge
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: VcfToolsVcfMerge



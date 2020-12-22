:orphan:

VcfTools: VcfConcat
=======================================

``VcfToolsVcfConcat`` · *1 contributor · 1 version*

Concatenates VCF files (for example split by chromosome). Note that the input and output VCFs will have the same number of columns, the script does not merge VCFs by position (see also vcf-merge).

In the basic mode it does not do anything fancy except for a sanity check that all files have the same columns. When run with the -s option, it will perform a partial merge sort, looking at limited number of open files simultaneously.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcftools.vcfconcat.versions import VcfToolsVcfConcat_0_1_16

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcftoolsvcfconcat_step",
           VcfToolsVcfConcat_0_1_16(
               vcfTabix=None,
           )
       )
       wf.output("out", source=vcftoolsvcfconcat_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for VcfToolsVcfConcat:

.. code-block:: bash

   # user inputs
   janis inputs VcfToolsVcfConcat > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcfTabix:
       - vcfTabix_0.vcf.gz
       - vcfTabix_1.vcf.gz




5. Run VcfToolsVcfConcat with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       VcfToolsVcfConcat





Information
------------

:ID: ``VcfToolsVcfConcat``
:URL: `http://vcftools.sourceforge.net/perl_module.html#vcf-concat <http://vcftools.sourceforge.net/perl_module.html#vcf-concat>`_
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

============  ===================  ============  ==========  =============================================================================
name          type                 prefix          position  documentation
============  ===================  ============  ==========  =============================================================================
vcfTabix      Array<Gzipped<VCF>>                        10
checkColumns  Optional<Boolean>    -c                        Do not concatenate, only check if the columns agree.
padMissing    Optional<Boolean>    -p                        Write '.' in place of missing columns. Useful for joining chrY with the rest.
mergeSort     Optional<Integer>    --merge-sort              Allow small overlaps in N consecutive files.
============  ===================  ============  ==========  =============================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task VcfToolsVcfConcat {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? checkColumns
       Boolean? padMissing
       Int? mergeSort
       Array[File] vcfTabix
       Array[File] vcfTabix_tbi
     }
     command <<<
       set -e
        vcf-concat \
         ~{if (defined(checkColumns) && select_first([checkColumns])) then "-c" else ""} \
         ~{if (defined(padMissing) && select_first([padMissing])) then "-p" else ""} \
         ~{if defined(mergeSort) then ("--merge-sort " + mergeSort) else ''} \
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
   label: 'VcfTools: VcfConcat'
   doc: |-
     Concatenates VCF files (for example split by chromosome). Note that the input and output VCFs will have the same number of columns, the script does not merge VCFs by position (see also vcf-merge).

     In the basic mode it does not do anything fancy except for a sanity check that all files have the same columns. When run with the -s option, it will perform a partial merge sort, looking at limited number of open files simultaneously.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/vcftools:v0.1.16-1-deb_cv1

   inputs:
   - id: checkColumns
     label: checkColumns
     doc: Do not concatenate, only check if the columns agree.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -c
   - id: padMissing
     label: padMissing
     doc: Write '.' in place of missing columns. Useful for joining chrY with the rest.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
   - id: mergeSort
     label: mergeSort
     doc: Allow small overlaps in N consecutive files.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --merge-sort
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
   - vcf-concat
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: VcfToolsVcfConcat



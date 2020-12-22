:orphan:

VcfLib: Vcf Filter
==============================

``vcffilter`` · *1 contributor · 1 version*

Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields

Operators can be any of: =, !, <, >, |, &

Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"

For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcffilter.versions import VcfFilter_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcffilter_step",
           VcfFilter_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcffilter_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcffilter:

.. code-block:: bash

   # user inputs
   janis inputs vcffilter > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcffilter with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcffilter





Information
------------

:ID: ``vcffilter``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-06-04
:Updated: 2020-06-04


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>  Filtered VCF
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

===============  =============================  =================  ==========  ===================================================================================================================================================================
name             type                           prefix               position  documentation
===============  =============================  =================  ==========  ===================================================================================================================================================================
vcf              VCF                                                        1  VCF to filter
info_filter      Optional<String>               --info-filter                  (-f) specifies a filter to apply to the info fields of records, removes alleles which do not pass the filter
genotype_filter  Optional<String>               --genotype-filter              (-g) specifies a filter to apply to the genotype fields of records
keep_info        Optional<Boolean>              --keep-info                    (-k) used in conjunction with '-g', keeps variant info, but removes genotype
filter_sites     Optional<Boolean>              --filter-sites                 (-s) filter entire records, not just alleles
tag_pass         Optional<String>               --tag-pass                     (-t) tag vcf records as positively filtered with this tag, print all records
tag_fail         Optional<String>               --tag-fail                     (-F) tag vcf records as negatively filtered with this tag, print all records
append_filter    Optional<Boolean>              --append-filter                (-A) append the existing filter tag, don't just replace it
allele_tag       Optional<String>               --allele-tag                   (-a) apply -t on a per-allele basis. adds or sets the corresponding INFO field tag
invert           Optional<Boolean>              --invert                       (-v) inverts the filter, e.g. grep -v
use_logical_or   Optional<Boolean>              --or                           (-o) use logical OR instead of AND to combine filters
region           Optional<Array<Gzipped<bed>>>  --region                       (-r) specify a region on which to target the filtering, requires a BGZF compressed file which has been indexed with tabix.  any number of regions may be specified.
===============  =============================  =================  ==========  ===================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcffilter {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? info_filter
       String? genotype_filter
       Boolean? keep_info
       Boolean? filter_sites
       String? tag_pass
       String? tag_fail
       Boolean? append_filter
       String? allele_tag
       Boolean? invert
       Boolean? use_logical_or
       Array[File]? region
       Array[File]? region_tbi
     }
     command <<<
       set -e
       vcffilter \
         ~{if defined(info_filter) then ("--info-filter '" + info_filter + "'") else ""} \
         ~{if defined(genotype_filter) then ("--genotype-filter '" + genotype_filter + "'") else ""} \
         ~{if (defined(keep_info) && select_first([keep_info])) then "--keep-info" else ""} \
         ~{if (defined(filter_sites) && select_first([filter_sites])) then "--filter-sites" else ""} \
         ~{if defined(tag_pass) then ("--tag-pass '" + tag_pass + "'") else ""} \
         ~{if defined(tag_fail) then ("--tag-fail '" + tag_fail + "'") else ""} \
         ~{if (defined(append_filter) && select_first([append_filter])) then "--append-filter" else ""} \
         ~{if defined(allele_tag) then ("--allele-tag '" + allele_tag + "'") else ""} \
         ~{if (defined(invert) && select_first([invert])) then "--invert" else ""} \
         ~{if (defined(use_logical_or) && select_first([use_logical_or])) then "--or" else ""} \
         ~{if (defined(region) && length(select_first([region])) > 0) then "--region '" + sep("' '", select_first([region])) + "'" else ""} \
         '~{vcf}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "shollizeck/vcflib:1.0.1"
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
   label: 'VcfLib: Vcf Filter'
   doc: |-
     Filter the specified vcf file using the set of filters.
     Filters are specified in the form "<ID> <operator> <value>:
      -f "DP > 10"  # for info fields
      -g "GT = 1|1" # for genotype fields
      -f "CpG"  # for 'flag' fields

     Operators can be any of: =, !, <, >, |, &

     Any number of filters may be specified.  They are combined via logical AND
     unless --or is specified on the command line.  Obtain logical negation through
     the use of parentheses, e.g. "! ( DP = 10 )"

     For convenience, you can specify "QUAL" to refer to the quality of the site, even
     though it does not appear in the INFO fields.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/vcflib:1.0.1

   inputs:
   - id: vcf
     label: vcf
     doc: VCF to filter
     type: File
     inputBinding:
       position: 1
   - id: info_filter
     label: info_filter
     doc: |-
       (-f) specifies a filter to apply to the info fields of records, removes alleles which do not pass the filter
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --info-filter
       separate: true
   - id: genotype_filter
     label: genotype_filter
     doc: (-g) specifies a filter to apply to the genotype fields of records
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genotype-filter
       separate: true
   - id: keep_info
     label: keep_info
     doc: (-k) used in conjunction with '-g', keeps variant info, but removes genotype
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep-info
       separate: true
   - id: filter_sites
     label: filter_sites
     doc: (-s) filter entire records, not just alleles
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --filter-sites
       separate: true
   - id: tag_pass
     label: tag_pass
     doc: (-t) tag vcf records as positively filtered with this tag, print all records
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tag-pass
       separate: true
   - id: tag_fail
     label: tag_fail
     doc: (-F) tag vcf records as negatively filtered with this tag, print all records
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tag-fail
       separate: true
   - id: append_filter
     label: append_filter
     doc: (-A) append the existing filter tag, don't just replace it
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --append-filter
       separate: true
   - id: allele_tag
     label: allele_tag
     doc: |-
       (-a) apply -t on a per-allele basis. adds or sets the corresponding INFO field tag
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --allele-tag
       separate: true
   - id: invert
     label: invert
     doc: (-v) inverts the filter, e.g. grep -v
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --invert
       separate: true
   - id: use_logical_or
     label: use_logical_or
     doc: (-o) use logical OR instead of AND to combine filters
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --or
       separate: true
   - id: region
     label: region
     doc: |-
       (-r) specify a region on which to target the filtering, requires a BGZF compressed file which has been indexed with tabix.  any number of regions may be specified.
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --region
       separate: true

   outputs:
   - id: out
     label: out
     doc: Filtered VCF
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcffilter
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcffilter



:orphan:

FilterVep
=========

``FilterVep`` · *1 contributor · 1 version*

#------------#
# filter_vep #
#------------#
http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.ensembl.filtervep.versions import FilterVep_98_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "filtervep_step",
           FilterVep_98_3(

           )
       )
       wf.output("out", source=filtervep_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FilterVep:

.. code-block:: bash

   # user inputs
   janis inputs FilterVep > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run FilterVep with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FilterVep





Information
------------

:ID: ``FilterVep``
:URL: *No URL to the documentation was provided*
:Versions: 98.3
:Container: quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0
:Authors: Michael Franklin
:Citations: None
:Created: 2020-05-26
:Updated: 2020-05-26


Outputs
-----------

======  ========  ===============
name    type      documentation
======  ========  ===============
out     TextFile
======  ========  ===============


Additional configuration (inputs)
---------------------------------

===============  =======================  =================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================
name             type                     prefix             position    documentation
===============  =======================  =================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================
input_file       Optional<File>           --input_file                   (-i) Specify the input file (i.e. the VEP results file). If no input file is specified, the script will attempt to read from STDIN. Input may be gzipped - to force the script to read a file as gzipped, use --gz
format           Optional<String>         --format                       [vcf|tab] Specify input file format (tab for any tab-delimited format, including default VEP output format)
outputFilename   Optional<Filename>       --output_file                  (-o) Specify the output file to write to. If no output file is specified, the script will write to STDOUT
force_overwrite  Optional<Boolean>        --force_overwrite              Force the script to overwrite the output file if it already exists
filter           Optional<Array<String>>  --filter                       (-f) Add filter. Multiple --filter flags may be used, and are treated as logical ANDs, i.e. all filters must pass for a line to be printed
list             Optional<Array<String>>  --list                         (-l) List allowed fields from the input file
count            Optional<Boolean>        --count                        (-c) Print only a count of matched lines
only_matched     Optional<Boolean>        --only_matched                 In VCF files, the CSQ field that contains the consequence data will often contain more than  one 'block' of consequence data, where each block corresponds to a variant/feature overlap. Using  filters. By default, the script prints out the entire VCF line if any of the blocks pass the filters.
vcf_info_field   Optional<String>         --vcf_info_field               With VCF input files, by default filter_vep expects to find VEP annotations encoded in the CSQ INFO key; VEP itself can be configured to write to a different key (with the equivalent --vcf_info_field flag). Use this flag to change the INFO key VEP expects to decode.
ontology         Optional<Boolean>        --ontology                     (-y) Use Sequence Ontology to match consequence terms. Use with operator 'is' to match against all child terms of your value. e.g. 'Consequence is coding_sequence_variant' will match missense_variant, synonymous_variant etc. Requires database connection; defaults to connecting to ensembldb.ensembl.org. Use --host, --port, --user, --version) connection parameters.
help             Optional<Boolean>        --help                         -h Print usage message and exit
===============  =======================  =================  ==========  =============================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task FilterVep {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? input_file
       String? format
       String? outputFilename
       Boolean? force_overwrite
       Array[String]? filter
       Array[String]? list
       Boolean? count
       Boolean? only_matched
       String? vcf_info_field
       Boolean? ontology
       Boolean? help
     }
     command <<<
       set -e
       filter_vep \
         ~{if defined(input_file) then ("--input_file '" + input_file + "'") else ""} \
         ~{if defined(format) then ("--format '" + format + "'") else ""} \
         --output_file '~{select_first([outputFilename, "generated.txt"])}' \
         ~{if (defined(force_overwrite) && select_first([force_overwrite])) then "--force_overwrite" else ""} \
         ~{if (defined(filter) && length(select_first([filter])) > 0) then "--filter '" + sep("' --filter '", select_first([filter])) + "'" else ""} \
         ~{if (defined(list) && length(select_first([list])) > 0) then "--list '" + sep("' '", select_first([list])) + "'" else ""} \
         ~{if (defined(count) && select_first([count])) then "--count" else ""} \
         ~{if (defined(only_matched) && select_first([only_matched])) then "--only_matched" else ""} \
         ~{if defined(vcf_info_field) then ("--vcf_info_field '" + vcf_info_field + "'") else ""} \
         ~{if (defined(ontology) && select_first([ontology])) then "--ontology" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: FilterVep
   doc: |
     #------------#
     # filter_vep #
     #------------#
     http://www.ensembl.org/info/docs/tools/vep/script/vep_filter.html

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0

   inputs:
   - id: input_file
     label: input_file
     doc: |-
       (-i) Specify the input file (i.e. the VEP results file). If no input file is specified, the script will attempt to read from STDIN. Input may be gzipped - to force the script to read a file as gzipped, use --gz
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --input_file
       separate: true
   - id: format
     label: format
     doc: |-
       [vcf|tab] Specify input file format (tab for any tab-delimited format, including default VEP output format)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --format
       separate: true
   - id: outputFilename
     label: outputFilename
     doc: |-
       (-o) Specify the output file to write to. If no output file is specified, the script will write to STDOUT
     type:
     - string
     - 'null'
     default: generated.txt
     inputBinding:
       prefix: --output_file
       separate: true
   - id: force_overwrite
     label: force_overwrite
     doc: Force the script to overwrite the output file if it already exists
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force_overwrite
       separate: true
   - id: filter
     label: filter
     doc: |-
       (-f) Add filter. Multiple --filter flags may be used, and are treated as logical ANDs, i.e. all filters must pass for a line to be printed
     type:
     - type: array
       inputBinding:
         prefix: --filter
         separate: true
       items: string
     - 'null'
     inputBinding: {}
   - id: list
     label: list
     doc: (-l) List allowed fields from the input file
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --list
       separate: true
   - id: count
     label: count
     doc: (-c) Print only a count of matched lines
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --count
       separate: true
   - id: only_matched
     label: only_matched
     doc: |-
       In VCF files, the CSQ field that contains the consequence data will often contain more than  one 'block' of consequence data, where each block corresponds to a variant/feature overlap. Using  filters. By default, the script prints out the entire VCF line if any of the blocks pass the filters.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --only_matched
       separate: true
   - id: vcf_info_field
     label: vcf_info_field
     doc: |-
       With VCF input files, by default filter_vep expects to find VEP annotations encoded in the CSQ INFO key; VEP itself can be configured to write to a different key (with the equivalent --vcf_info_field flag). Use this flag to change the INFO key VEP expects to decode.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --vcf_info_field
       separate: true
   - id: ontology
     label: ontology
     doc: |-
       (-y) Use Sequence Ontology to match consequence terms. Use with operator 'is' to match against all child terms of your value. e.g. 'Consequence is coding_sequence_variant' will match missense_variant, synonymous_variant etc. Requires database connection; defaults to connecting to ensembldb.ensembl.org. Use --host, --port, --user, --version) connection parameters.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ontology
       separate: true
   - id: help
     label: help
     doc: -h Print usage message and exit
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
       separate: true

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - filter_vep
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: FilterVep



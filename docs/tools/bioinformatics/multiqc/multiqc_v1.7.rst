:orphan:

Multiqc
=================

``MultiQC`` · *1 contributor · 1 version*

Usage: multiqc [OPTIONS] <analysis directory>
MultiQC aggregates results from bioinformatics analyses across many samples into a single report.
It searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, 
perfect for summarising the output from numerous bioinformatics tools.
To run, supply with one or more directory to scan for analysis results. To run here, use 'multiqc .'

Author: Phil Ewels (http://phil.ewels.co.uk)



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.multiqc.versions import Multiqc_v1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "multiqc_step",
           Multiqc_v1_7(
               directory=None,
           )
       )

    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for MultiQC:

.. code-block:: bash

   # user inputs
   janis inputs MultiQC > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       directory: null




5. Run MultiQC with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       MultiQC





Information
------------

:ID: ``MultiQC``
:URL: `http://multiqc.info <http://multiqc.info>`_
:Versions: v1.7
:Container: ewels/multiqc:v1.7
:Authors: Michael Franklin
:Citations: None
:Created: 2019-10-24
:Updated: 2019-10-24


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  =======================  ==================  ==========  ============================================================================================
name            type                     prefix                position  documentation
==============  =======================  ==================  ==========  ============================================================================================
directory       Directory                                             1  searches a given directory for analysis logs and compiles a HTML report
force           Optional<Boolean>        --force                         (-f) Overwrite any existing reports
dirs            Optional<String>         --dirs                          (-d) Prepend directory to sample names
dirsDepth       Optional<Integer>        --dirs-depth                    (-dd) Prepend [INT] directories to sample names. Negative number to take from start of path.
fullnames       Optional<Boolean>        --fullnames                     (-s) Do not clean the sample names (leave as full file name)
title           Optional<String>         --title                         (-i) Report title. Printed as page header, used for filename if not otherwise specified.
comment         Optional<String>         --comment                       (-b) Custom comment, will be printed at the top of the report.
filename        Optional<Filename>       --filename                      (-n) Report filename. Use 'stdout' to print to standard out.
outdir          Optional<Filename>       --outdir                        (-o) Create report in the specified output directory.
template        Optional<String>         --template                      (-t)  Report template to use.
tag             Optional<String>         --tag                           Use only modules which tagged with this keyword, eg. RNA
view_tags       Optional<Boolean>        --view_tags                     (--view-tags) View the available tags and which modules they load
ignore          Optional<Boolean>        --ignore                        (-x) Ignore analysis files (glob expression)
ignoreSamples   Optional<Boolean>        --ignore-samples                Ignore sample names (glob expression)
ignoreSymlinks  Optional<Boolean>        --ignore-symlinks               Ignore symlinked directories and files
sampleNames     Optional<File>           --sample-names                  File containing alternative sample names
exclude         Optional<Array<String>>  --exclude                       (-e) Do not use this module. Can specify multiple times.
module          Optional<Array<String>>  --module                        (-m) Use only this module. Can specify multiple times.
dataDir         Optional<Boolean>        --data-dir                      Force the parsed data directory to be created.
noDataDir       Optional<Boolean>        --no-data-dir                   Prevent the parsed data directory from being created.
dataFormat      Optional<String>         --data-format                   (-k)  Output parsed data in a different format. Default: tsv
export          Optional<Boolean>        --export                        (-p) Export plots as static images in addition to the report
flat            Optional<Boolean>        --flat                          (-fp) Use only flat plots (static images)
interactive     Optional<Boolean>        --interactive                   (-ip) Use only interactive plots (HighCharts Javascript)
lint            Optional<Boolean>        --lint                          Use strict linting (validation) to help code development
pdf             Optional<Boolean>        --pdf                           Creates PDF report with 'simple' template. Requires Pandoc to be installed.
noMegaqcUpload  Optional<Boolean>        --no-megaqc-upload              Don't upload generated report to MegaQC, even if MegaQC options are found
config          Optional<File>           --config                        (-c) Specific config file to load, after those in MultiQC dir / home dir / working dir.
cl_config       Optional<File>           --cl_config                     (--cl-config) Specify MultiQC config YAML on the command line
verbose         Optional<Boolean>        --verbose                       (-v) Increase output verbosity.
quiet           Optional<Boolean>        --quiet                         (-q) Only show log warnings
==============  =======================  ==================  ==========  ============================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task MultiQC {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Directory directory
       Boolean? force
       String? dirs
       Int? dirsDepth
       Boolean? fullnames
       String? title
       String? comment
       String? filename
       String? outdir
       String? template
       String? tag
       Boolean? view_tags
       Boolean? ignore
       Boolean? ignoreSamples
       Boolean? ignoreSymlinks
       File? sampleNames
       Array[String]? exclude
       Array[String]? module
       Boolean? dataDir
       Boolean? noDataDir
       String? dataFormat
       Boolean? export
       Boolean? flat
       Boolean? interactive
       Boolean? lint
       Boolean? pdf
       Boolean? noMegaqcUpload
       File? config
       File? cl_config
       Boolean? verbose
       Boolean? quiet
     }
     command <<<
       set -e
       multiqc \
         ~{if (defined(force) && select_first([force])) then "--force" else ""} \
         ~{if defined(dirs) then ("--dirs '" + dirs + "'") else ""} \
         ~{if defined(dirsDepth) then ("--dirs-depth " + dirsDepth) else ''} \
         ~{if (defined(fullnames) && select_first([fullnames])) then "--fullnames" else ""} \
         ~{if defined(title) then ("--title '" + title + "'") else ""} \
         ~{if defined(comment) then ("--comment '" + comment + "'") else ""} \
         --filename '~{select_first([filename, "generated"])}' \
         --outdir '~{select_first([outdir, "generated"])}' \
         ~{if defined(template) then ("--template '" + template + "'") else ""} \
         ~{if defined(tag) then ("--tag '" + tag + "'") else ""} \
         ~{if (defined(view_tags) && select_first([view_tags])) then "--view_tags" else ""} \
         ~{if (defined(ignore) && select_first([ignore])) then "--ignore" else ""} \
         ~{if (defined(ignoreSamples) && select_first([ignoreSamples])) then "--ignore-samples" else ""} \
         ~{if (defined(ignoreSymlinks) && select_first([ignoreSymlinks])) then "--ignore-symlinks" else ""} \
         ~{if defined(sampleNames) then ("--sample-names '" + sampleNames + "'") else ""} \
         ~{if (defined(exclude) && length(select_first([exclude])) > 0) then "--exclude '" + sep("' --exclude '", select_first([exclude])) + "'" else ""} \
         ~{if (defined(module) && length(select_first([module])) > 0) then "--module '" + sep("' --module '", select_first([module])) + "'" else ""} \
         ~{if (defined(dataDir) && select_first([dataDir])) then "--data-dir" else ""} \
         ~{if (defined(noDataDir) && select_first([noDataDir])) then "--no-data-dir" else ""} \
         ~{if defined(dataFormat) then ("--data-format '" + dataFormat + "'") else ""} \
         ~{if (defined(export) && select_first([export])) then "--export" else ""} \
         ~{if (defined(flat) && select_first([flat])) then "--flat" else ""} \
         ~{if (defined(interactive) && select_first([interactive])) then "--interactive" else ""} \
         ~{if (defined(lint) && select_first([lint])) then "--lint" else ""} \
         ~{if (defined(pdf) && select_first([pdf])) then "--pdf" else ""} \
         ~{if (defined(noMegaqcUpload) && select_first([noMegaqcUpload])) then "--no-megaqc-upload" else ""} \
         ~{if defined(config) then ("--config '" + config + "'") else ""} \
         ~{if defined(cl_config) then ("--cl_config '" + cl_config + "'") else ""} \
         ~{if (defined(verbose) && select_first([verbose])) then "--verbose" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         '~{directory}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "ewels/multiqc:v1.7"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Multiqc
   doc: |
     Usage: multiqc [OPTIONS] <analysis directory>
     MultiQC aggregates results from bioinformatics analyses across many samples into a single report.
     It searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, 
     perfect for summarising the output from numerous bioinformatics tools.
     To run, supply with one or more directory to scan for analysis results. To run here, use 'multiqc .'

     Author: Phil Ewels (http://phil.ewels.co.uk)

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ewels/multiqc:v1.7

   inputs:
   - id: directory
     label: directory
     doc: searches a given directory for analysis logs and compiles a HTML report
     type: Directory
     inputBinding:
       position: 1
   - id: force
     label: force
     doc: (-f) Overwrite any existing reports
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force
       separate: true
   - id: dirs
     label: dirs
     doc: (-d) Prepend directory to sample names
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --dirs
       separate: true
   - id: dirsDepth
     label: dirsDepth
     doc: |-
       (-dd) Prepend [INT] directories to sample names. Negative number to take from start of path.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --dirs-depth
       separate: true
   - id: fullnames
     label: fullnames
     doc: (-s) Do not clean the sample names (leave as full file name)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --fullnames
       separate: true
   - id: title
     label: title
     doc: |-
       (-i) Report title. Printed as page header, used for filename if not otherwise specified.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --title
       separate: true
   - id: comment
     label: comment
     doc: (-b) Custom comment, will be printed at the top of the report.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --comment
       separate: true
   - id: filename
     label: filename
     doc: (-n) Report filename. Use 'stdout' to print to standard out.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --filename
       separate: true
   - id: outdir
     label: outdir
     doc: (-o) Create report in the specified output directory.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --outdir
       separate: true
   - id: template
     label: template
     doc: (-t)  Report template to use.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --template
       separate: true
   - id: tag
     label: tag
     doc: Use only modules which tagged with this keyword, eg. RNA
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tag
       separate: true
   - id: view_tags
     label: view_tags
     doc: (--view-tags) View the available tags and which modules they load
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --view_tags
       separate: true
   - id: ignore
     label: ignore
     doc: (-x) Ignore analysis files (glob expression)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore
       separate: true
   - id: ignoreSamples
     label: ignoreSamples
     doc: Ignore sample names (glob expression)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-samples
       separate: true
   - id: ignoreSymlinks
     label: ignoreSymlinks
     doc: Ignore symlinked directories and files
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-symlinks
       separate: true
   - id: sampleNames
     label: sampleNames
     doc: File containing alternative sample names
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --sample-names
       separate: true
   - id: exclude
     label: exclude
     doc: (-e) Do not use this module. Can specify multiple times.
     type:
     - type: array
       inputBinding:
         prefix: --exclude
         separate: true
       items: string
     - 'null'
     inputBinding: {}
   - id: module
     label: module
     doc: (-m) Use only this module. Can specify multiple times.
     type:
     - type: array
       inputBinding:
         prefix: --module
         separate: true
       items: string
     - 'null'
     inputBinding: {}
   - id: dataDir
     label: dataDir
     doc: Force the parsed data directory to be created.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --data-dir
       separate: true
   - id: noDataDir
     label: noDataDir
     doc: Prevent the parsed data directory from being created.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-data-dir
       separate: true
   - id: dataFormat
     label: dataFormat
     doc: '(-k)  Output parsed data in a different format. Default: tsv'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --data-format
       separate: true
   - id: export
     label: export
     doc: (-p) Export plots as static images in addition to the report
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --export
       separate: true
   - id: flat
     label: flat
     doc: (-fp) Use only flat plots (static images)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --flat
       separate: true
   - id: interactive
     label: interactive
     doc: (-ip) Use only interactive plots (HighCharts Javascript)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --interactive
       separate: true
   - id: lint
     label: lint
     doc: Use strict linting (validation) to help code development
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --lint
       separate: true
   - id: pdf
     label: pdf
     doc: Creates PDF report with 'simple' template. Requires Pandoc to be installed.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pdf
       separate: true
   - id: noMegaqcUpload
     label: noMegaqcUpload
     doc: Don't upload generated report to MegaQC, even if MegaQC options are found
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-megaqc-upload
       separate: true
   - id: config
     label: config
     doc: |-
       (-c) Specific config file to load, after those in MultiQC dir / home dir / working dir.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --config
       separate: true
   - id: cl_config
     label: cl_config
     doc: (--cl-config) Specify MultiQC config YAML on the command line
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --cl_config
       separate: true
   - id: verbose
     label: verbose
     doc: (-v) Increase output verbosity.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbose
       separate: true
   - id: quiet
     label: quiet
     doc: (-q) Only show log warnings
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet
       separate: true

   outputs: []
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - multiqc
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: MultiQC



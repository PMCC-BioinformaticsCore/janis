:orphan:

Multiqc
=================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


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
:URL: *No URL to the documentation was provided*
:Versions: v1.7
:Container: ewels/multiqc:v1.7
:Authors: 
:Citations: None
:Created: None
:Updated: None



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

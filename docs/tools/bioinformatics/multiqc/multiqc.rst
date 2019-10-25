:orphan:

Multiqc
=================

0 contributors · 1 version

:ID: ``MultiQC``
:Python: ``janis_bioinformatics.tools.multiqc.versions import Multiqc_v1_7``
:Versions: v1.7
:Container: ewels/multiqc:v1.7
:Authors: 
:Citations: None
:Created: None
:Updated: None
:Required inputs:
   - ``directory: Directory``
:Outputs: 


Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

------

Additional configuration (inputs)
---------------------------------

==============  =======================  ============================================================================================
name            type                     documentation
==============  =======================  ============================================================================================
directory       Directory                searches a given directory for analysis logs and compiles a HTML report
force           Optional<Boolean>        (-f) Overwrite any existing reports
dirs            Optional<String>         (-d) Prepend directory to sample names
dirsDepth       Optional<Integer>        (-dd) Prepend [INT] directories to sample names. Negative number to take from start of path.
fullnames       Optional<Boolean>        (-s) Do not clean the sample names (leave as full file name)
title           Optional<String>         (-i) Report title. Printed as page header, used for filename if not otherwise specified.
comment         Optional<String>         (-b) Custom comment, will be printed at the top of the report.
filename        Optional<Filename>       (-n) Report filename. Use 'stdout' to print to standard out.
outdir          Optional<Filename>       (-o) Create report in the specified output directory.
template        Optional<String>         (-t)  Report template to use.
tag             Optional<String>         Use only modules which tagged with this keyword, eg. RNA
view_tags       Optional<Boolean>        (--view-tags) View the available tags and which modules they load
ignore          Optional<Boolean>        (-x) Ignore analysis files (glob expression)
ignoreSamples   Optional<Boolean>        Ignore sample names (glob expression)
ignoreSymlinks  Optional<Boolean>        Ignore symlinked directories and files
sampleNames     Optional<File>           File containing alternative sample names
exclude         Optional<Array<String>>  (-e) Do not use this module. Can specify multiple times.
module          Optional<Array<String>>  (-m) Use only this module. Can specify multiple times.
dataDir         Optional<Boolean>        Force the parsed data directory to be created.
noDataDir       Optional<Boolean>        Prevent the parsed data directory from being created.
dataFormat      Optional<String>         (-k)  Output parsed data in a different format. Default: tsv
export          Optional<Boolean>        (-p) Export plots as static images in addition to the report
flat            Optional<Boolean>        (-fp) Use only flat plots (static images)
interactive     Optional<Boolean>        (-ip) Use only interactive plots (HighCharts Javascript)
lint            Optional<Boolean>        Use strict linting (validation) to help code development
pdf             Optional<Boolean>        Creates PDF report with 'simple' template. Requires Pandoc to be installed.
noMegaqcUpload  Optional<Boolean>        Don't upload generated report to MegaQC, even if MegaQC options are found
config          Optional<File>           (-c) Specific config file to load, after those in MultiQC dir / home dir / working dir.
cl_config       Optional<File>           (--cl-config) Specify MultiQC config YAML on the command line
verbose         Optional<Boolean>        (-v) Increase output verbosity.
quiet           Optional<Boolean>        (-q) Only show log warnings
==============  =======================  ============================================================================================


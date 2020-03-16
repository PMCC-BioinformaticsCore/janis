:orphan:

CNVKit
======

*0 contributors · 1 version*

:ID: ``CNVKit``
:Python: ``janis_bioinformatics.tools.ucsf.cnvkit.cnvkit_0_9_6 import CNVKit_0_9_6``
:Versions: 0.9.6
:Container: etal/cnvkit:0.9.6
:Authors: 
:Citations: Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014). CNVkit: Genome-wide copy number detection and visualization from targeted sequencing. PLOS Computational Biology 12(4):e1004873
:DOI: 10.1371/journal.pcbi.1004873
:Created: 2019-07-03 00:00:00
:Updated: 2019-07-03 00:00:00
:Required inputs:
   - ``reference: File``
:Outputs: 


Documentation
-------------

URL: `https://github.com/etal/cnvkit <https://github.com/etal/cnvkit>`_


        A command-line toolkit and Python library for detecting copy number variants 
        and alterations genome-wide from high-throughput sequencing.

------

None

Additional configuration (inputs)
---------------------------------

===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
name             type                prefix               position    documentation
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
reference        File                --reference                      REFERENCE Copy number reference file (.cnn).
outputDirectory  Optional<Filename>  --output-dir                     DIRECTORY Output directory.
method           Optional<String>    --method                         (-m) {hybrid,amplicon,wgs} Sequencing protocol: hybridization capture ('hybrid'), targeted amplicon sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines whether and how to use antitarget bins. [Default: hybrid]
maleReference    Optional<String>    --male-reference                 (-y, --haploid-x-reference) Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).
countReads       Optional<String>    --count-reads                    (-c) Get read depths by counting read midpoints within each bin. (An alternative algorithm).
dropLowCoverage  Optional<String>    --drop-low-coverage              Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
processes        Optional<String>    --processes                      (-p) [PROCESSES] Number of subprocesses used to running each of the BAM files in parallel. Without an argument, use the maximum number of available CPUs. [Default: process each BAM in serial]
rscriptPath      Optional<String>    --rscript-path                   Path to the Rscript excecutable to use for running R code. Use this option to specify a non-default R installation. [Default: Rscript]
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================


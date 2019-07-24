
.. include:: cnvkit_0.9.6

CNVKit
======

Description
-------------

Tool identifier: ``CNVKit``

Tool path: ``janis_bioinformatics.tools.ucsf.cnvkit.cnvkit_0_9_6 import CNVKit_0_9_6``

Version: 0.9.6

Docker: ``etal/cnvkit:0.9.6``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  ======  ===========  ==========  ============================================
name       type    prefix       position    documentation
=========  ======  ===========  ==========  ============================================
reference  File    --reference              REFERENCE Copy number reference file (.cnn).
=========  ======  ===========  ==========  ============================================

Optional inputs
***************

===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
name             type                prefix               position    documentation
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
outputDirectory  Optional<Filename>  --output-dir                     DIRECTORY Output directory.
method           Optional<String>    --method                         (-m) {hybrid,amplicon,wgs} Sequencing protocol: hybridization capture ('hybrid'), targeted amplicon sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines whether and how to use antitarget bins. [Default: hybrid]
maleReference    Optional<String>    --male-reference                 (-y, --haploid-x-reference) Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).
countReads       Optional<String>    --count-reads                    (-c) Get read depths by counting read midpoints within each bin. (An alternative algorithm).
dropLowCoverage  Optional<String>    --drop-low-coverage              Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
processes        Optional<String>    --processes                      (-p) [PROCESSES] Number of subprocesses used to running each of the BAM files in parallel. Without an argument, use the maximum number of available CPUs. [Default: process each BAM in serial]
rscriptPath      Optional<String>    --rscript-path                   Path to the Rscript excecutable to use for running R code. Use this option to specify a non-default R installation. [Default: Rscript]
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*CNVKit was last updated on **Unknown***.
*This page was automatically generated on 2019-07-24*.

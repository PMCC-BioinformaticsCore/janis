
GATK4: Apply base quality score recalibration
==============================================================
Tool identifier: ``GATK4ApplyBQSR``

Tool path: ``from janis_bioinformatics.tools.gatk4 import Gatk4ApplyBqsr_4_0``

Documentation
-------------

Docker
******
``broadinstitute/gatk:4.0.12.0``

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php>`_

Docstring
*********
Apply base quality score recalibration: This tool performs the second pass in a two-stage 
process called Base Quality Score Recalibration (BQSR). Specifically, it recalibrates the 
base qualities of the input reads based on the recalibration table produced by the 
BaseRecalibrator tool, and outputs a recalibrated BAM or CRAM file.

Summary of the BQSR procedure: The goal of this procedure is to correct for systematic bias 
that affect the assignment of base quality scores by the sequencer. The first pass consists 
of calculating error empirically and finding patterns in how error varies with basecall 
features over all bases. The relevant observations are written to a recalibration table. 
The second pass consists of applying numerical corrections to each individual basecall 
based on the patterns identified in the first step (recorded in the recalibration table) 
and write out the recalibrated data to a new BAM or CRAM file.

- This tool replaces the use of PrintReads for the application of base quality score 
    recalibration as practiced in earlier versions of GATK (2.x and 3.x).
- You should only run ApplyBQSR with the covariates table created from the input BAM or CRAM file(s).
- Original qualities can be retained in the output file under the "OQ" tag if desired. 
    See the `--emit-original-quals` argument for details.

Outputs
-------
======  =======  ===============
name    type     documentation
======  =======  ===============
out     BamPair
======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ========  ==========  =======================================
name       type           prefix      position  documentation
=========  =============  ========  ==========  =======================================
bam        BamPair        -I                10  The SAM/BAM/CRAM file containing reads.
reference  FastaWithDict  -R                    Reference sequence
=========  =============  ========  ==========  =======================================

Optional inputs
***************

==============  ==================  =================  ==========  =============================================================
name            type                prefix               position  documentation
==============  ==================  =================  ==========  =============================================================
outputFilename  Optional<Filename>  -O                             Write output to this file
recalFile       Optional<tsv>       --bqsr-recal-file              Input recalibration table for BQSR
intervals       Optional<bed>       --intervals                    -L (BASE) One or more genomic intervals over which to operate
tmpDir          Optional<String>    --tmp-dir                  11  Temp directory to use.
==============  ==================  =================  ==========  =============================================================


Metadata
********

Author: Michael Franklin


*GATK4: Apply base quality score recalibration was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-09*.

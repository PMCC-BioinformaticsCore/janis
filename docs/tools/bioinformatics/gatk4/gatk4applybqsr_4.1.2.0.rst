:orphan:

GATK4: Apply base quality score recalibration
==============================================================

*1 contributor · 4 versions*

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

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.applybqsr.versions import Gatk4ApplyBqsr_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4applybqsr_step",
           Gatk4ApplyBQSR(
               bam=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4applybqsr_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4ApplyBQSR:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4ApplyBQSR > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       reference: reference.fasta




5. Run Gatk4ApplyBQSR with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4ApplyBQSR





Information
------------


:ID: ``Gatk4ApplyBQSR``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  =================  ==========  =============================================================
name            type                prefix               position  documentation
==============  ==================  =================  ==========  =============================================================
bam             IndexedBam          -I                         10  The SAM/BAM/CRAM file containing reads.
reference       FastaWithIndexes    -R                             Reference sequence
outputFilename  Optional<Filename>  -O                             Write output to this file
recalFile       Optional<tsv>       --bqsr-recal-file              Input recalibration table for BQSR
intervals       Optional<bed>       --intervals                    -L (BASE) One or more genomic intervals over which to operate
tmpDir          Optional<String>    --tmp-dir                  11  Temp directory to use.
==============  ==================  =================  ==========  =============================================================

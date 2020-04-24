:orphan:

GATK4: MuTect2
=============================

*1 contributor · 4 versions*

Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV) 
and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic 
genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 
for a step-by-step description of the workflow and Article#11127 for an overview of what traditional 
somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory. 
Although we present the tool for somatic calling, it may apply to other contexts, 
such as mitochondrial variant calling.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.mutect2.versions import GatkMutect2_4_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4mutect2_step",
           GatkMutect2_4_0(
               tumor=None,
               tumorName=None,
               normal=None,
               normalName=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4mutect2_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4Mutect2:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4Mutect2 > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal: normal.bam
       normalName: <value>
       reference: reference.fasta
       tumor: tumor.bam
       tumorName: <value>




5. Run Gatk4Mutect2 with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4Mutect2





Information
------------


:ID: ``Gatk4Mutect2``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  ====================  =================
name    type                  documentation
======  ====================  =================
out     CompressedIndexedVCF  To determine type
======  ====================  =================



Additional configuration (inputs)
---------------------------------

========================  ====================  ===============================  ==========  ==============================================================================================================================================================
name                      type                  prefix                             position  documentation
========================  ====================  ===============================  ==========  ==============================================================================================================================================================
tumor                     IndexedBam            -I                                        6  BAM/SAM/CRAM file containing reads
tumorName                 String                -tumor                                    6  BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode.
normal                    IndexedBam            -I                                        5  BAM/SAM/CRAM file containing reads
normalName                String                -normal                                   6  BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode.
reference                 FastaWithIndexes      -R                                        8  Reference sequence file
intervals                 Optional<bed>         -L                                        7  One or more genomic intervals over which to operate
outputFilename            Optional<Filename>    -O                                       20
germlineResource          Optional<IndexedVCF>  --germline-resource                      10
afOfAllelesNotInResource  Optional<Float>       --af-of-alleles-not-in-resource          11  Population allele fraction assigned to alleles not found in germline resource. Please see docs/mutect/mutect2.pdf fora derivation of the default value.
panelOfNormals            Optional<IndexedVCF>  --panel-of-normals                       10  A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
========================  ====================  ===============================  ==========  ==============================================================================================================================================================

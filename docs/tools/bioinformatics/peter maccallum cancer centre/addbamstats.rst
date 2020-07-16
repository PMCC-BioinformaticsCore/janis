:orphan:

Add Bam Statistics to Vcf
=======================================

*1 contributor · 2 versions*

usage: add_bam_stats.py [-h] -i I -o O --type {germline,somatic}
                        [--mpileup MPILEUP] [--normal_mpileup NORMAL_MPILEUP]
                        [--tumor_mpileup TUMOR_MPILEUP]
                        [--normal_id NORMAL_ID] [--tumor_id TUMOR_ID]

Get stats from bam file and write to vcf

required arguments:
  -i I                  input vcf
  -o O                  output vcf
  --type {germline,somatic}
                        must be either germline or somatic
  --mpileup MPILEUP     mpileup file extracted from bam file
  --normal_mpileup NORMAL_MPILEUP
                        mpileup file extracted from the normal sample bam,
                        required if input is somatic vcf
  --tumor_mpileup TUMOR_MPILEUP
                        mpileup file extracted from the tumor sample, required
                        if input is somatic vcf
  --normal_id NORMAL_ID
                        Normal sample id, required if input is somatic vcf
  --tumor_id TUMOR_ID   Tumor sample id, required if input is somatic vcf

optional arguments:
  -h, --help            show this help message and exit
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addbamstats.versions import AddBamStats_dev

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addbamstats_step",
           AddBamStats_dev(
               inputVcf=None,
               type=None,
           )
       )
       wf.output("out", source=addbamstats_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for addBamStats:

.. code-block:: bash

   # user inputs
   janis inputs addBamStats > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputVcf: inputVcf.vcf
       type: <value>




5. Run addBamStats with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       addBamStats





Information
------------


:ID: ``addBamStats``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_
:Versions: dev, 0.0.7
:Container: jyu/pmacutil:dev
:Authors: Jiaan Yu
:Citations: None
:Created: None
:Updated: 2020-05-20 00:00:00



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ================  ==========  ===================================================================================
name            type                prefix            position    documentation
==============  ==================  ================  ==========  ===================================================================================
inputVcf        VCF                 -i                            input vcf
type            String              --type                        must be either germline or somatic
mpileup         Optional<File>      --mpileup                     mpileup file extracted from bam file
normalMpileup   Optional<File>      --normal_mpileup              mpileup file extracted from the normal sample bam, required if input is somatic vcf
tumorMpileup    Optional<File>      --tumor_mpileup               mpileup file extracted from the tumor sample bam, required if input is somatic vcf
normalID        Optional<String>    --normal_id                   normal sample id, required if input is somatic vcf
tumorID         Optional<String>    --tumor_id                    tumor sample id, required if input is somatic vcf
outputFilename  Optional<Filename>  -o                            output vcf name
==============  ==================  ================  ==========  ===================================================================================

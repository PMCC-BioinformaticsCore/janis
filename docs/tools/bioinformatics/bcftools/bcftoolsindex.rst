:orphan:

BCFTools: Index
===============================

*0 contributors · 1 version*

Index bgzip compressed VCF/BCF files for random access.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.index.versions import BcfToolsIndex_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolsindex_step",
           BcfToolsIndex_1_9(
               vcf=None,
           )
       )
       wf.output("out", source=bcftoolsindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolsIndex:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolsIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run bcftoolsIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolsIndex





Information
------------


:ID: ``bcftoolsIndex``
:URL: `https://samtools.github.io/bcftools/bcftools.html#norm <https://samtools.github.io/bcftools/bcftools.html#norm>`_
:Versions: v1.9
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
:Authors: 
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: 2019-01-24
:Updated: None



Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
out     CompressedIndexedVCF
======  ====================  ===============



Additional configuration (inputs)
---------------------------------

========  =================  ===========  ==========  ============================================================
name      type               prefix         position  documentation
========  =================  ===========  ==========  ============================================================
vcf       CompressedVCF                            1
csi       Optional<Boolean>  --csi                    (-c) generate CSI-format index for VCF/BCF files [default]
force     Optional<Boolean>  --force                  (-f) overwrite index if it already exists
minShift  Optional<Integer>  --min-shift              (-m) set minimal interval size for CSI indices to 2^INT [14]
tbi       Optional<Boolean>  --tbi                    (-t) generate TBI-format index for VCF files
threads   Optional<Integer>  --threads                sets the number of threads [0]
nrecords  Optional<Boolean>  --nrecords               (-n) print number of records based on existing index file
stats     Optional<Boolean>  --stats                  (-s) print per contig stats based on existing index file
========  =================  ===========  ==========  ============================================================

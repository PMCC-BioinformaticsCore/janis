:orphan:

Gene Coverage Per Sample
================================================

*1 contributor · 3 versions*

usage: gene_coverage_per_sample.py [-h] [-l LIST] [-n NAME] [-p PATH] [-b BED]
                                   [-g GENE] [-r REGION] [-f FOLDS] [-d]
                                   [-t THREADS]

Gene or region coverage of bam

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  List file: A tsv file contains SampleName
                        PathToBedtoolsOutput on each line
  -n NAME, --name NAME  Sample name if list not used
  -p PATH, --path PATH  Path to bedtools output if list not used
  -b BED, --bed BED     (Deprecated option) Bed file
  -g GENE, --gene GENE  Output gene file
  -r REGION, --region REGION
                        Output region file
  -f FOLDS, --folds FOLDS
                        Folds, quoted and commna sepparated, default
                        1,10,20,100
  -d, --remove_duplicates
                        (Deprecated option) Remove marked duplicates in
                        analysis, default:false
  -t THREADS, --threads THREADS
                        number of threads, default:32
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.genecovpersample.versions import GeneCoveragePerSample_dev

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "genecoveragepersample_step",
           GeneCoveragePerSample_dev(

           )
       )
       wf.output("geneFileOut", source=genecoveragepersample_step.geneFileOut)
       wf.output("regionFileOut", source=genecoveragepersample_step.regionFileOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for geneCoveragePerSample:

.. code-block:: bash

   # user inputs
   janis inputs geneCoveragePerSample > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run geneCoveragePerSample with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       geneCoveragePerSample





Information
------------


:ID: ``geneCoveragePerSample``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance>`_
:Versions: dev, 0.0.8, 0.0.7
:Container: jyu/pmacutil:dev
:Authors: Jiaan Yu
:Citations: None
:Created: None
:Updated: 2020-04-03 00:00:00



Outputs
-----------

=============  ========  ===============
name           type      documentation
=============  ========  ===============
geneFileOut    TextFile
regionFileOut  TextFile
=============  ========  ===============



Additional configuration (inputs)
---------------------------------

==================  ==================  =========  ==========  ========================================================
name                type                prefix     position    documentation
==================  ==================  =========  ==========  ========================================================
listFile            Optional<File>      --list                 List file: A tsv file contains SampleName	PathToBedtoolsOutput on each line
sampleName          Optional<String>    --name                 Sample name if list not used
bedtoolsOutputPath  Optional<File>      --path                 Path to bedtools output if list not used
outputGeneFile      Optional<Filename>  --gene                 Output gene file
outputRegionFile    Optional<Filename>  --region               Output region file
fold                Optional<String>    --fold                 Folds, quoted and commna sepparated, default 1,10,20,100
threads             Optional<Integer>   --threads              number of threads, default:32
==================  ==================  =========  ==========  ========================================================

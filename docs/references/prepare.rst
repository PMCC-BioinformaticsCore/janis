Janis Prepare
###############

- # * = - .

``prepare`` is functionality of Janis to improve the process of running pipelines. Specifically, Janis prepare will perform a few key actions:

- Downloads reference datasets
- Performs simple transforms data into the correct format (eg: VCF -> gVCF)
- Checks quality of some inputs:
    - for example, contigs in a bed file match what's in a declared reference

.. note::

   **Prepare** only works with Janis pipelines


Quickstart
************

The ``janis prepare`` command line is almost exactly the same as the ``janis run``. You should supply it with  inputs (either through an inputs yaml or on the command line) and any other configuration options. You must supply an output directory (or declare an ``output_dir`` in your janis config), the job file and a run script is written to this directory. The job file is also always written to stdout.

For example:

.. code-block:: bash

   # Write inputs file
   cat <<EOT >> inputs.yaml
   sampleName: NA12878
   fastqs:
   - - /<fastqdata>/WGS_30X_R1.fastq.gz
     - /<fastqdata>//WGS_30X_R2.fastq.gz
   EOT

   # run janis prepare
   janis prepare \
       -o "$HOME/janis/WGSGermlineGATK-run/" \
       --inputs inputs.yaml \
       --source-hint hg38 \
       WGSGermlineGATK

This will:
- Write an inputs file to disk,
- download all the hg38 reference files
- transform any data types that might need to be transformed, eg:
    - ``gridss blacklist`` requires a bed, but the source hint gives a gzipped bed (``ENCFF001TDO.bed.gz``)
    - ``snps_dbsnp`` wants a compressed and tabix indexed VCF, but the source hint gives a regular VCF.
    - ``reference`` build the appropriate indexes for the hg38 assembly (as these aren't downloaded by default)
- Perform some sanity checks on the data you've provided, eg:
    - the contigs in the gridss_blacklist will be checked against those found in the assembly's reference.
- Write a job file and run script into the output directory.

Downloading reference datasets
*******************************

A ``source`` is declared on an input through the input documentation, there are two methods for doing this:

1. Single file - this file is always used regardless of the inputs
2. A dictionary of *source hints*. This might be useful for specifying a reference type, eg hg38 or hg19. However this does involve the pipeline author finding public datasets for each *hint*.

Example:

.. code-block:: python

   import janis_core as j

   myworkflow = j.WorkflowBuilder("wf")

   # Input that can be used for ALL workflows
   myworkflow.input(
       "inp",
       File,
       doc=j.InputDocumentation(
           source="https://janis.readthedocs.io/en/latest/references/prepare.html"
       )
   )

   # Input that is only localised if the hint is specified
   mworkflow.input(
       "gridss_blacklist",
       File,
       doc=j.InputDocumentation(
           source={
               "hg19": "https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz",
               "hg38": "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz",
           }
       ),
   )

Secondary files and localising sources
========================================

By default, secondary files are downloaded with the primary file. The pipeline author can optionally request that secondary files be not downloaded with ``skip_sourcing_secondary_files=True``. In particular, this is important for the exemplar WGS pipelines, because the BWA indexes (``".amb", ".ann", ".bwt", ".pac", ".sa"``) were generated with a different version of BWA.

Although the Janis Transformation step of the janis prepare will perform the re-index, this could take a few hours. If you can speed up this step, please `raise an issue <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues/new>` or open a Pull Request!!

Types of reference paths
=========================

Janis can localise remote paths through its `FileScheme` mechanic. Currently, Janis supports the following fileschemes:

- Local (useful for keeping local references to a pipeline definition when sharing a pipeline internally)
- HTTP ``http:// OR http://`` prefix (using a GET request)
- GCS ``gs://`` prefix (public buckets only)

Requires more work:
- S3: `Implementation required <https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/fbcefc665a9a23073830acb8f87cc41ea58bcf8f/janis_assistant/management/filescheme.py#L483-L499>`_


Input Documentation Source
============================

An ``InputDocumentation`` class can be used to document the different options about an input. See the `:class:janis.InputDocumentation` initialiser below. You can supply it directly to the doc field of an input, but you can also provide a dictionary which will get converted into an InputDocumentation, for example:

.. code-block:: python

   w = j.WorkflowBuilder("wf")

   # using the InputDocumentation class
   w.input("inp1", str, doc=j.InputDocumentation(doc="This is inp1", quality=j.InputQualityType.user))

   # Use a dictionary
   w.input("inp2", str, doc={"doc": "This is inp2", "quality": "user"})

.. autoclass:: janis.InputDocumentation
   :members: __init__


Transformations
****************

A strong benefit of the rich data types that Janis provides, means that we can do perform basic transformations to prepare
your data for the pipeline. Simply, we've taught Janis how to perform some basic transformations (using tools in Janis!),
and then Janis can determine how to convert your data, for example:

If you provide a VCF, and the pipeline requires a gzipped and tabix'd VCF, we construct a prepare stage that performs this for you:

.. code-block:: text

   VCF -> VcfGZ -> VcfGzTabix

More information

- `The Conversion Problem <https://janis.readthedocs.io/en/latest/casestudy/conversionproblem.html?highlight=transformation>`_ (Janis guide)
- Bioinformatics transformations `github:janis-bioinformatics/transformations/__init__.py <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/blob/master/janis_bioinformatics/transformations/__init__.py>`_

.. note::

   These Janis transformations are performed AFTER the reference localisation, so janis can transform your downloaded file to the correct format if possible.

Quality Checks
*****************

These are a number of fairly custom checks to catch some frequent errors that we've seen:

- Input checker: makes sure your inputs can be found (only works for local paths)
- `Contig checker <https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/master/janis_assistant/modifiers/contigchecker.py>`_: If you define a single input of ``FASTA`` type and any number of inputs with ``BED`` type, Janis will check that the contigs you declare in the ``BEDs`` are found in the ``.fasta.fai`` index. This check only returns warnings, and will NOT stop you from running a workflow.

Janis performs some of these checks on every run, some are only performed on the `janis prepare`.

.. raw:: html

   <hr />

Developer notes
***********************

The Janis prepare steps are all implemented using ``PipelineModifiers``:

- ``FileFinderLocatorModifier``
- ``InputFileQualifierModifier``
- ``InputTransformerModifier``
- ``InputChecker``
- ``ContigChecker``

If entering through the `prepare` cli method, the ``run_prepare_processing`` flag is set which initialises a custom set of pipeline modifiers to execute while preparing the job file.

Janis Transformations
========================

Janis transformations are fairly simple, they're defined in the relevant tool registry (eg: janis-bioinformatics), and exposed through the ``janis.datatype_transformations`` entrypoint. These JanisTransformations are added to a graph, and then we just perform a breadth first search on this graph looking for the shortest number of steps to connect two data types. Transformations are directional, and no logic is performed to evaluate the *weight* or *effort* of a step.

`The Conversion Problem <https://janis.readthedocs.io/en/latest/casestudy/conversionproblem.html?highlight=transformation>`_ (Janis guide) describes JanisTransformations in a blog sort of style.


Job File
==========

Although not strictly related to Janis Prepare, it was an important change that was made for janis-prepare to work correctly. Functionally, a ``PreparedJob`` describes everything needed to run a workflow (except the workflow). It's an amalgamation of different sections of the janis configuration, adjusted inputs and runtime configurations. It's serializable, which means it can automatically be written to disk and parsed back in.

A ``run.sh`` script can be generated, which just calls the job file with the workflow reference, something like:

.. code-block:: text

   # This script was automatically generated by Janis on 2021-01-13 11:14:40.121673.

   janis run \
       -j /Users/franklinmichael/janis/hello/20210113_111440/janis/job.yaml \
       hello
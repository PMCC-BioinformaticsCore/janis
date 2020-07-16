:orphan:

Sequenza: bam2seqz
=================================

*2 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.sequenza.bam2seqz.versions import SequenzaBam2Seqz_3_0_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "seqzbam2seqz_step",
           SequenzaBam2Seqz_3_0_0(
               normal=None,
               tumour=None,
               wiggle_file=None,
               fasta_reference=None,
           )
       )
       wf.output("out", source=seqzbam2seqz_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SeqzBam2Seqz:

.. code-block:: bash

   # user inputs
   janis inputs SeqzBam2Seqz > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fasta_reference: fasta_reference.fasta
       normal: normal.bam
       tumour: tumour.bam
       wiggle_file: wiggle_file




5. Run SeqzBam2Seqz with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SeqzBam2Seqz





Information
------------


:ID: ``SeqzBam2Seqz``
:URL: `http://www.cbs.dtu.dk/biotools/sequenza/ <http://www.cbs.dtu.dk/biotools/sequenza/>`_
:Versions: 3.0.0, 2.2.0.9000
:Container: sequenza/sequenza:3.0.0
:Authors: mumbler, evanwehi
:Citations: None
:Created: 2019-12-10
:Updated: 2019-12-10



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============



Additional configuration (inputs)
---------------------------------

===============  ==================  ========  ==========  ==============================================================================================
name             type                prefix      position  documentation
===============  ==================  ========  ==========  ==============================================================================================
normal           IndexedBam          --normal           2  Name of the BAM/pileup file from the reference/normal sample
tumour           IndexedBam          --tumor            4  Name of the BAM/pileup file from the reference/normal sample
wiggle_file      File                -gc                6  The GC-content wiggle file
fasta_reference  FastaFai            --fasta            8  The reference FASTA file used to generate the intermediate pileup. Required when input are BAM
output_filename  Optional<Filename>  --output          10  Name of the output file. To use gzip compression name the file ending in .gz. Default STDOUT.
===============  ==================  ========  ==========  ==============================================================================================

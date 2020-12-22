
Fastq
=====

FASTQ files are text files containing sequence data with quality score, there are different typeswith no standard: https://www.drive5.com/usearch/manual/fastq_files.html



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fastq import Fastq

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastq", Fastq(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.

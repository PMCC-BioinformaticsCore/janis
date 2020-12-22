
FastqGz
=======

FastqGz files are compressed sequence data with quality score, there are different typeswith no standard: https://en.wikipedia.org/wiki/FASTQ_format



Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.data_types.fastq import FastqGz

   w = WorkflowBuilder("my_workflow")

   w.input("input_fastqgz", FastqGz(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.

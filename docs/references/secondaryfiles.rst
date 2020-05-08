Secondary / Accessory Files
=============================

In some domains (looking specifically at Bioinformatics here), a single file isn't enough to contain all the information. Janis borrows the concept of secondary files from the CWL specification, and in fact we use the same pattern for grouping these files.

Often a *secondary* or *accessory* file is used to provide additional information, potentially a quick access index. These files are attached to the original file by a specific file pattern. We'll talk about that more soon.


For this reason, Janis allows data_types that inherit from a ``File`` to specify a ``secondary_file`` list of files to be bundled with.


Secondary file pattern
-----------------------

As earlier mentioned, we follow the `Common Workflow Language secondary file pattern <https://github.com/common-workflow-language/common-workflow-language/blob/a062055fddcc7d7d9dbc53d28288e3ccb9a800d8/v1.0/Process.yml#L465-L473>`_:

1. If string begins with one or more caret ``^`` characters, for each caret, remove the last file extension from the path (the last period ``.`` and all following characters).  If there are no file extensions, the path is unchanged.
2. Append the remainder of the string to the end of the file path.

Examples
*********

- `IndexedBam <https://janis.readthedocs.io/en/latest/datatypes/indexedbam.html>`_
    - Pattern: ``.bai``
    - Files:
        - Base: ``myfile.bam``
        - ``myfile.bam.bai``

- `FastaWithIndexes <https://janis.readthedocs.io/en/latest/datatypes/fastawithindexes.html>`_:
    - Pattern: ``.amb``, ``.ann``, ``.bwt``, ``.pac``, ``.sa,`` ``.fai``, ``^.dict``
    - Files:
        - Base: ``reference.fasta``
        - ``reference.fasta.amb``
        - ``reference.fasta.ann``
        - ``reference.fasta.bwt``
        - ``reference.fasta.pac``
        - ``reference.fasta.sa,``
        - ``reference.fasta.fai``
        - ``reference.dict``


Proposed
*********

Implement optional secondary files as per `CWL v1.1 <https://github.com/common-workflow-language/cwl-v1.1/blob/565eb43ac8836505c5fd57b0e8378f2988f61626/Process.yml#L525-L535>`_.



Implementation note
----------------------

CWL
*******

As we mimic the CWL secondary file pattern, we don't need to do any extra work except by providing this pattern to a:

- `CommandInputParameter <https://www.commonwl.org/v1.0/CommandLineTool.html#CommandInputParameter>`_
- `CommandOutputParameter <https://www.commonwl.org/v1.0/CommandLineTool.html#CommandOutputParameter>`_
- `InputParameter <https://www.commonwl.org/v1.0/Workflow.html#InputParameter>`_
- `WorkflowOutputParameter <https://www.commonwl.org/v1.0/Workflow.html#WorkflowOutputParameter>`_

If you use the ``secondaries_present_as`` on a :class:`janis.ToolInput` or :class:`janis.ToolOutput`, a CWL expression is generated to rename the secondary file expression. More information can be found about this in `common-workflow-language/cwltool#1232 <https://github.com/common-workflow-language/cwltool/issues/1232>`_.

Issues
+++++++++

CWLTool has an issue when attempting to scatter using multiple fields (using the ``dotproduct`` or ``*_crossproduct`` methods), more information can be found on `common-workflow-language/cwltool#1208 <https://github.com/common-workflow-language/cwltool/issues/1208>`_.


WDL
********

The translation for WDL to implement secondary files was one of the most challenging aspects of the translation. Notably, WDL has no concept of secondary files. There are a few things we had to consider:

- Every file needs to be individually localised.
- A data type with secondary files can be used in an array of inputs
- Secondary files may need to be globbed if used as an Output data type
- An array of files with secondaries can be scattered on (including scattered by multiple fields)
- Janis should fill the input job with these secondary files (with the correct extension)



Implementation
+++++++++++++++

Let's just break this down into different sections

Case 1: Simple index
########################

The following ``workflow.input("my_bam", BamBai)``, definition when connected to a tool might look like the following

.. code-block:: none

   workflow WGSGermlineGATK {
     input {
       File my_bam
       File my_bam_bai

     }
     call my_tool {
       input:
         bam=my_bam
         bam_bai=my_bam_bai
     }
     output {
       File out_bam = my_tool.out
       File out_bam_bai = my_tool.out_bai
     }
   }



Note the extra annotations and mappings fot the ``bai`` type.

Case 2: Array of inputs with simple scatter
#############################################

This is modification of the first example, nb: this isn't full functional workflow code:

.. code-block:: python

   workflow.input("my_bams", Array(BamBai))

   workflow.step(
       "my_step",
       MyTool(bam=workflow.my_bams),
       scatter="bam"
   )

Might result in the following workflow:

.. code-block:: none

   workflow WGSGermlineGATK {
     input {
       Array[File] my_bams
       Array[File] my_bams_bai

     }
     scatter (Q in zip(my_bams, my_bams_bai)) {
       call my_tool as my_step {
         input:
           bam=Q.left
           bam_bai=Q.right
       }
     }

     output {
       Array[File] out_bams = my_tool_that_accepts_array.out
       Array[File] out_bams_bai = my_tool_that_accepts_array.out_bai
     }
   }



Case 3: Multiple array inputs, scattering by multiple fields
#################################################################

Consider the following workflow:

.. code-block:: python

    workflow.input("my_bams", Array(BamBai))
    workflow.input("my_references", Array(FastaBwa))

    workflow.step(
        "my_step",
        ToolTypeThatAcceptsMultipleBioinfTypes(
            bam=workflow.my_bams, reference=workflow.my_references
        ),
        scatter=["bam", "reference"],
    )

    workflow.output("out_bam", source=workflow.my_step.out_bam)
    workflow.output("out_reference", source=workflow.my_step.out_reference)

This gets complicated quickly:


.. code-block:: none
   
   workflow scattered_bioinf_complex {
     input {
       Array[File] my_bams
       Array[File] my_bams_bai
       Array[File] my_references
       Array[File] my_references_amb
       Array[File] my_references_ann
       Array[File] my_references_bwt
       Array[File] my_references_pac
       Array[File] my_references_sa
     }
     scatter (Q in zip(transpose([my_bams, my_bams_bai]), transpose([my_references, my_references_amb, my_references_ann, my_references_bwt, my_references_pac, my_references_sa]))) {
        call MyTool as my_step {
         input:
           bam=Q.left[0],
           bam_bai=Q.left[1],
           reference=Q.right[0],
           reference_amb=Q.right[1],
           reference_ann=Q.right[2],
           reference_bwt=Q.right[3],
           reference_pac=Q.right[4],
           reference_sa=Q.right[5]
       }
     }
     output {
       Array[File] out_bam = my_step.out_bam
       Array[File] out_reference = my_step.out_reference
     }
   }


Known limitations
+++++++++++++++++++++++


- There is no namespace collision:
    - Two files with similar prefixes but differences in punctuation will clash
    - A second input that is suffixed with the secondary’s extension will clash: eg: mybam_bai will clash with mybam with a secondary of .bai.
- Globbing a secondary file might not be possible when the original file extension is unknown. There are 2 considerations for this:
    - Subclasses of File should caller super() with the expected extension
    - Globbing based on a generated Filename (through InputSelector), will consider the extension property.

Relevant WDL issues:

- `broadinstitute/cromwell#2269 (Secondary index files and directories in WDL)  <https://github.com/broadinstitute/cromwell/issues/2269>`_
- `openwdl/wdl#289 (File Bundles and Secondary / Accessory Files) <https://github.com/openwdl/wdl/issues/289>`_
- `GATK Forums#9299 (Secondary index files and directories in WDL) <https://gatkforums.broadinstitute.org/wdl/discussion/9299/secondary-index-files-and-directories-in-wdl>`_

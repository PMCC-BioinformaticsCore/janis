Collecting tool outputs
#########################

This guide will roughly walk you through the different ways of collecting outputs from a CommandTool.
The bread and butter of this tutorial is the ToolOutput, and how you use various `Selectors / Operators <https://janis.readthedocs.io/en/latest/references/listoperators.html>`_.

Examples:

- :ref:`ex-stderrstdout`
- :ref:`ex-wildcardsglobs`
- :ref:`ex-namedoutputs`

.. autoclass:: janis.ToolOutput
   :members: __init__

You should note that there are key differences between how strings are coerced into Files / Directories in CWL and WDL.

- In WDL, a string is automatically coercible to a file, where the path is relative to the execution directory
- In CWL, a path is NOT automatically coercible, and instead a FILE object (``{path: "<path>", class: "File / Directory"}``) must be created. Janis shortcuts this instead, by inserting your strings as a globs, and letting CWL do this. There may be unintended side effects of this process.



Convention
=============

We'll presume in this workflow that you've imported Janis like the following:

.. code-block:: python

   import janis_core as j

Examples
===========

.. _ex-stderrstdout:

Stderr / Stdout
*****************

Collecting stdout and stderr can be done by simply annotating the types. This is functionally equivalent to type File, and using Stderr / Stdout as a selector:

.. code-block:: python

   outputs=[
       # stdout
       j.ToolOutput("out_stdout_1", j.Stdout()),
       j.ToolOutput("out_stdout_2", j.File(), selector=j.Stdout()),
       # stderr
       j.ToolOutput("out_stderr_1", j.Stderr()),
       j.ToolOutput("out_stderr_2", j.File(), selector=j.Stderr()),
   ]

.. _ex-wildcardsglobs:

Wildcard / glob outputs
*************************

If it's not practical or impossible to determine the names of the outputs, you can use a :class:`janis.WildcardSelector` to find all the files that match a particular pattern. This glob pattern is not transformed, and differences may occur between CWL / WDL depending on what glob syntax they use - please refer to their individual documentation for more information

- CWL: Globs
- WDL: Globs

You can use a glob in Janis with:

.. code-block:: python

   outputs=[
       j.ToolOutput("out_text_files", j.Array(j.File), selector=j.WildcardSelector("*.txt")),
       # the next two are functionally equivalent
       j.ToolOutput("out_single_text_file_1", j.Array(j.File), selector=j.WildcardSelector("*.txt", select_first=True)),
        j.ToolOutput("out_single_text_file_2", j.Array(j.File), selector=j.WildcardSelector("*.txt")[0])
   ]


Roughly, this is translated to the following:

WDL:

.. code-block:: text

   Array[File] out_txt_files = glob("*.txt")
   Array[File] out_single_txt_file_1 = glob("*.txt")[0]
   Array[File] out_single_txt_file_2 = glob("*.txt")[0]

CWL:

.. code-block:: text
   - id: out_optional_glob
     label: out_optional_glob
     type:
     - File
     - 'null'
     outputBinding:
       glob: '*.txt'
       loadContents: false
   - id: out_single_csv_file_1
     label: out_single_csv_file_1
     type:
       type: array
       items: File
     outputBinding:
       glob: '*.csv'
       loadContents: false
   - id: out_single_csv_file_2
     label: out_single_csv_file_2
     type:
       type: array
       items: File
     outputBinding:
       glob: '*.csv'
       outputEval: $(self[0])
       loadContents: false

.. _ex-namedoutputs:

Named outputs
**********************

Often we'll use a string or a Filename generator to name an output of a tool. For example, ``samtools sort`` accepts an argument ``-o`` which is an output filename, on top of the regular ``"bam"`` input. We want our output to be of type Bam, so we'll use a :class:`janis.Filename` class, this accepts a few arguments: ``prefix``, `suffix`` and ``extension``, and will generate a filename based on these attributes.

We want our filename to be based on the input bam to keep consistency in our naming, so let's choose the following attributes:

- ``prefix`` - will be the Bam, but want the file extension removed (this will automatically take the basename)
- ``suffix`` - `".sorted"``
- ``extension`` - ``.bam``

We can create the following ToolInput value to match this (we use a ToolInput as it means we could override it later):

.. code-block:: python

   ToolInput(
       "outputFilename",
       Filename(
           prefix=InputSelector("bam", remove_file_extension=True),
           suffix=".sorted",
           extension=".bam",
       ),
       position=1,     # Ensure it appears before the "bam" input
       prefix="-o",    # Prefix, eg: '-o <output-filename>'
   )

Then, on the ToolOutput, we can use the selector ``selector=InputSelector("outputFilename")`` to get this value. This results in the final tool:

.. code-block:: python

   SamtoolsSort_1_9_0 = CommandToolBuilder(
       tool="SamToolsSort",
       base_command=["samtools", "sort"],
       inputs=[
           ToolInput("bam", input_type=Bam(), position=2),
           ToolInput(
               "outputFilename",
               Filename(
                   prefix=InputSelector("bam", remove_file_extension=True),
                   suffix=".sorted",
                   extension=".bam",
               ),
               position=1,
               prefix="-o",
           ),
       ],
       outputs=[ToolOutput("out", Bam(), selector=InputSelector("outputFilename"))],
       container="quay.io/biocontainers/samtools:1.9--h8571acd_11",
       version="1.9.0",
   )


Looking at the relevant WDL:

.. code-block:: text

   task SamToolsSort {
     input {
       File bam
       String? outputFilename
     }
     command <<<
       set -e
       samtools sort \
         -o '~{select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])}' \
         '~{bam}'
     >>>
     output {
       File out = select_first([outputFilename, "~{basename(bam, ".bam")}.sorted.bam"])
     }
   }


And CWL:

.. code-block:: text

   class: CommandLineTool
   id: SamToolsSort
   baseCommand:
   - samtools
   - sort

   inputs:
   - id: bam
   type: File
   inputBinding:
   position: 2
   - id: outputFilename
   type:
   - string
   - 'null'
   default: generated.sorted.bam
   inputBinding:
   prefix: -o
   position: 1
   valueFrom: $(inputs.bam.basename.replace(/.bam$/, "")).sorted.bam

   outputs:
   - id: out
   type: File
   outputBinding:
   glob: $(inputs.bam.basename.replace(/.bam$/, "")).sorted.bam

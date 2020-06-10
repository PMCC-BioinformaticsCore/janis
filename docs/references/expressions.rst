Expressions
#############

.. note::

   This feature is only available in Janis-Core v0.9.x and above.

Expressions in Janis refer to a set of features that improve how you can manipulate and massage inputs into an appropriate structure.

In this set of changes, we consider two important types:

- :class:`janis.Selector` A placeholder used for selecting other values (without modification).
- :class:`janis.Operator` Analogous to a function call.

To remain abstract, Janis builds up a set of these operations, all of which have a CWL and WDL analogue. We'll talk more about how these are converted later.

Selectors
==============

As initially described, a Selector can be seen as a _placeholder_ for referencing other values. Some common types of selectors:

- InputSelector - references the input by a string

    - MemorySelector - Special InputSelector for selecting memory in GBs
    - CpuSelector - Special InputSelector for getting number of CPUs that a task will request.
- InputNodeSelector - specifically references the input node, only available in a workflow
- StepOutputSelector - references the output of a step, only available in a workflow.
- WildcardSelector - Collect a number of files with a pattern. Currently the use of this pattern isn't well defined, but is passed without modification to CWL and WDL. Only available on a ToolOutput.

Selectors participate in the typesystem, and as such have a dynamic return type.

Operators
====================

Operators are equivalent to the concept of functions. They take a set of inputs, and generate some output.
Operators inherit from selectors, The inputs and outputs are typed, and hence operators and selectors particpate in the typing system.

Both are more flexible with the types that they provide, though some selectors (notable InputSelector) may be unable to give a definite type.

In terms of their implementation, it's the operator's responsibility to convert itself to CWL and WDL. There will be an example below.

In some cases (where possible), we've overridden the standard python operators to give the Janis analogue.

- ``__neg__``
- ``__and__`` NB: this is the `bitwise AND <https://docs.python.org/3/library/operator.html#operator.and_>`_ and not the logical ``and``.
- ``__rand__`` NB: similar to __and__
- ``__or__`` NB: this is the `bitwise OR <https://docs.python.org/3/library/operator.html#operator.or_>`_ and not the logical ``or``.
- ``__ror__``
- ``__add__``
- ``__radd__``
- ``__sub__``
- ``__rsub__``
- ``__mul__``
- ``__rmul__``
- ``__truediv__``
- ``__rtruediv__``
- ``__eq__``
- ``__ne__``
- ``__gt__``
- ``__ge__``
- ``__lt__``
- ``__le__``
- ``__len__``
- ``__getitem__``

- ``to_string_formatter()``
- ``as_str()``
- ``as_bool()``
- ``as_int()``
- ``op_and()``
- ``op_or()``
- ``basename()``


List of operators
++++++++++++++++++

.. note::

   See the `List of operators <https://janis.readthedocs.io/en/latest/references/listoperators.html>`_ guide for more information.

Example usage
===================

An operator's usage should be as you'd expect, let's see an example use of the FileSizeOperator, highlighting:

- Import the operator from ``janis_core.operators.standard``
- Applying the operation directly on the``echo.inp`` step input field.
- Converting the transformed input to a string with ``.as_str()``

.. code:: python

   from janis_core import WorkflowBuilder, File
   from janis_core.operators.standard import FileSizeOperator
   from janis_unix.tools import Echo

   w = WorkflowBuilder("sizetest")
   w.input("fileInp", File)

   w.step("print",
       Echo(inp=(FileSizeOperator(w.fileInp) * 1024).as_str())
   )
   w.output("out", source=w.print)


Before we go any further, let's look at the WDL and CWL translations:

WDL
+++++

We can see how the step input expression is converted directly inline.

.. code-block:: none

   version development

   import "tools/echo.wdl" as E

   workflow sizetest {
     input {
       File fileInp
     }
     call E.echo as print {
       input:
         inp=((1024 * size(fileInp, "MB")))
     }
     output {
       File out = print.out
     }
   }

CWL
+++++

The CWL translation is a little bit trickier due to the way scope of ``valueFrom`` expressions. We can see
that our variable is passed into the scope (prefixed by an underscore), and then the valueFrom contains the
expression that we produced - this is the value that the ``print`` step will see.

.. code-block:: yaml

   #!/usr/bin/env cwl-runner
   class: Workflow
   cwlVersion: v1.0

   requirements:
     InlineJavascriptRequirement: {}
     StepInputExpressionRequirement: {}

   inputs:
  inp:
       type: File

   outputs:
  out:
       type: File
       outputSource: print/out

  steps:
    print:
      label: Echo
      in:
        _print_inp_fileInp:
          source: fileInp
        inp:
          valueFrom: $(String((1024 * (inputs._print_inp_fileInp.size / 1048576))))
      run: tools/echo.cwl
      out:
      - out
  id: sizetest

Implementation notes
=====================

Let's first look at the implementation of the ``janis_core.operators.standard.FileSizeOperator``:

For example, we could consider the implementation of the ``FileSizeOperator``:

.. code:: python

   class FileSizeOperator(Operator):
       """
       Returned in MB: Note that this does NOT include the reference files (yet)
       """

       def argtypes(self):
           return [File()]

       def returntype(self):
           return Float

       def __str__(self):
           f = self.args[0]
           return f"file_size({f})"

       def to_wdl(self, unwrap_operator, *args):
           f = unwrap_operator(self.args[0])
           return f'size({f}, "MB")'

       def to_cwl(self, unwrap_operator, *args):
           f = unwrap_operator(self.args[0])
           return f"({f}.size / 1048576)"

- The argtypes returns an array of the types of the arguments that we expect. FileSizeOperator only expects one arg of a File type
- The ReturnType is a singular field which depicts the return type of the function.
- ``to_wdl`` is the function that builds our WDL function, it calls ``unwrap_operator`` on the argument (to ensure that any tokens are unwrapped), and then builds the command line using string interpolation.
- ``to_cwl`` operates exactly the same, except we use ES5 javascript to build our operation. The ``/ 1048576`` is to ensure that the value we receive in Bytes is converted to Megabytes (MB).



PROPOSED
========

In order to keep the current spec scoped, there is additional functionality that is planned:


- Using operators to build up the CPU, Memory and future resource values
- Custom python operations
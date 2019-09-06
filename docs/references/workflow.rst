Workflow
========

*Manages the connections between tools*

.. note::
	There have been large changes, backwards incompatible changes to the syntax for workflows in v0.6.0.

Declaration
###########

.. autoclass:: janis.Workflow


Overview
########

The :class:`janis.Workflow` class exposes inputs, and manages the connections between these inputs, tools and exposes some outputs.


A workflow does not directly execute, but declares what inputs a `janis.CommandTool` should receive.

A representation of a workflow can be exported to ``cwl`` or ``wdl``
through the :method:`janis.Workflow.translate()` function.

------------

Translating
***********

Currently Janis supports two translation targets:

1. `Common Workflow Language <cwl.html>`_
2. `Workflow Description Language <wdl.html>`_

.. automethod:: janis.Workflow.translate

------------

Structure of a workflow
***********************

A workflow has the following _nodes_:

- Inputs - :class:`janis.Workflow.input()`
- Steps - :class:`janis.Workflow.step()`
- Outputs - :class:`janis.Workflow.output()`

Once an node has been added to the workflow, it may be referenced through dot-notation on the workflow. For this reason, identifiers have certain naming restrictions.

Creating an input
*********************

An input requires a unique identifier (string) and a :class:`janis.DataType`.

.. automethod:: janis.Workflow.input

The input node is returned from this function, and is also available as a property on a workflow (accessible through dot-notation OR index notation).

.. code-block:: python

    import janis as j

    w = j.Workflow("myworkflow")
    myInput = w.input("myInput", String)
    myInput == w.myInput == w["myInput"] # True

Creating a step
*********************

A step requires a unique identifier (string), a tool (either a :class:`janis.CommandTool` or :class:`janis.Workflow`), scattering information (if required), and a mapping of all the inputs.

.. automethod:: janis.Workflow.step

Janis will throw an error if all the *required* inputs are not provided. There may

.. code-block:: python

    from janis.unix.tools.echo import Echo

    # Echo has the required input: "inp": String
    # https://janis.readthedocs.io/en/latest/tools/unix/echo.html

    echoStep = w.step("echoStep", Echo, inp=w.myInput)
    echoStep == w.echoStep == w["echoStep"] # True

Creating an output
*********************

An output requires a unique identifier (string), an output source and an *optional* :class:`janis.DataType`. If a data type is provided, it is type-checked against the output source.

.. automethod:: janis.Workflow.output

You are unable to connect an input node to an output node, and an output node cannot be referenced as a step input.


.. code-block:: python

    # w.echoStep added to workflow
    w.output("out", source=w.echoStep)


Examples
..........

    The `Echo <https://janis.readthedocs.io/en/latest/tools/unix/echo.html>`_ tool
    has one inputs ``inp`` of type string, and one output ``out``.

.. code-block:: python

    import janis as j
    from janis.unix.tools.echo import Echo

    w = j.Workflow("myworkflow")
    w.input("myInput", String)
    echoStep = w.step("echoStep", Echo, inp=w.myInput)
    w.output("out", source=w.echoStep)

    # Will print the CWL, input file and relevant tools to the console
    w.translate("cwl", to_disk=False)  # or "wdl"

------------

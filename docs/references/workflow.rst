Workflow
========

*Manages the connections between inputs, steps and outputs.*

Declaration
###########

.. autoclass:: janis.Workflow


Overview
########

The ``Workflow`` class manages the connections between
an :class:`janis.Input`, :class:`janis.Step` and an :class:`janis.Output`,
called an :class:`janis.Edge`. These edges are automatically wrapped
in a :class:`janis.Node` to become a :class:`janis.InputNode`,
:class:`janis.StepNode` and :class:`janis.OutputNode`.


A workflow does not directly execute,
but declares what inputs a `janis.CommandTool` should receive.

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

Creating connections
*********************

An edge requires a *start* and *finish* tuple. Janis is flexible in the way
this *start* and *finish* can be specified. Note that a step can have
multiple inputs and outputs, so it's important to make sure these are fully
qualified. Workflow level inputs and outputs do not need further qualifiers.

The most efficient way is to provide a *start* or *finish* is as the actual
input, step or output classes. The step input can be qualified through dot
notation. By providing these nodes directly to the ``add_edge`` function, you
do **not** need to separately add them to the Workflow first.

    *This dot notation is possible due to a Python* ``__getattr__`` *override
    on the* :class:`janis.Step` to return the correct format.

It's also possible to provide a '/' (forward slash) delimited string,
eg: ``"{stepId}/{tag}"`` or simple ``"{inputId}"``, however you must
ensure the inputs, steps and outputs are added to your Workflow first
through the :class:`janis.Workflow.add_items` method.

.. automethod:: janis.Workflow.add_edge
.. automethod:: janis.Workflow.add_edges

Examples
..........

    The `Echo <https://janis.readthedocs.io/en/latest/tools/unix/echo.html>`_ tool
    has one inputs ``inp`` of type string, and one output ``out``.

Example: https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/edges.py

.. code-block:: python

    from janis.unix.tools.echo import Echo

    w = Workflow("workflowId")

    inp = Input("inputId", data_type=String())
    stp = Step("stepId", tool=Echo())
    out = Output("outputId")

    # Method 1 - Preferred
    w.add_edge(inp, stp.inp)
    w.add_edges([
        (stp.out, out)
    ])

    # Method 2 - Prone to issues when IDs change
    w.add_items(inp, stp, out)
    w.add_edge("inputId", "stepId/inp")
    w.add_edges([
        ("stepId/out", "outputId")
    ])


------------

Adding items to Workflow
*********************

In order to provide hints, a Workflow must know about all the items being
added to it. The ``add_edge`` method returns an ordered list of nodes,
corresponding to the arguments.

.. automethod:: janis.Workflow.add_items
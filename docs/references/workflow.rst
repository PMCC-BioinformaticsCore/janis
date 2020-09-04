Workflow
========

*Manages the connections between tools*

Declaration
###########

There are two major ways to construct a workflow:

- Inline using the :class:`janis.WorkflowBuilder` class,
- or Inheriting from the :class:`janis.Workflow` class and implementing the required methods.

.. autoclass:: janis.WorkflowBuilder

.. autoclass:: janis.Workflow

Advanced Workflows
*******************

Janis allows you to dynamically create workflows based on inputs. More information can be found on the `Dynamic Workflows <dynamicworkflows.html>`_ page.

Overview
########

The :class:`janis.Workflow` and :class:`janis.WorkflowBuilder` classes exposes inputs, and manages the connections between these inputs, tools and exposes some outputs.

A :class:`janis.WorkflowBuilder` is the class used inline to declare workflows. The :class:`janis.Workflow` class should only be inherited through subclasses.

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

Once an node has been added to the workflow, it may be referenced through dot-notation on the workflow. For this reason, identifiers have certain naming restrictions. In the following examples we're going to create an *inline* workflow using the ``WorkflowBuilder`` class.

Creating an input
*********************

An input requires a unique identifier (string) and a :class:`janis.DataType`.

.. automethod:: janis.Workflow.input

The input node is returned from this function, and is also available as a property on a workflow (accessible through dot-notation OR index notation).

.. code-block:: python

    import janis as j

    w = j.WorkflowBuilder("myworkflow")
    myInput = w.input("myInput", String)
    myInput == w.myInput == w["myInput"] # True


.. note::

   **Default vs Value**: The input


Creating a step
*********************

A step requires a unique identifier (string), a mapped tool (either a :class:`janis.CommandTool` or :class:`janis.Workflow` called with it's inputs), scattering information (if required).

.. automethod:: janis.Workflow.step

Janis will throw an error if all the *required* inputs are not provided. You can provide the parameter ``ignore_missing=True`` to the step function to skip this check.

.. code-block:: python

    from janis.unix.tools.echo import Echo

    # Echo has the required input: "inp": String
    # https://janis.readthedocs.io/en/latest/tools/unix/echo.html

    echoStep = w.step("echoStep", Echo(inp=w.myInput))
    echoStep == w.echoStep == w["echoStep"] # True

Creating an output
*********************

An output requires a unique identifier (string), an output source and an *optional* :class:`janis.DataType`. If a data
type is provided, it is type-checked against the output source. Don't be put off by the automatically generated
interface for the output method, it's there to be exhaustive for the type definitions.

Here is the (simplified) method definition:

.. code-block:: python

   def output(
       self,
       identifier: str,
       datatype: Optional[ParseableType] = None,
       source: Union[Selector, ConnectionSource]=None # or List[Selector, ConnectionSource]
       output_folder: Union[str, Selector, List[Union[str, Selector]]] = None,
       output_name: Union[bool, str, Selector, ConnectionSource] = True, # let janis decide output name
       extension: Optional[str] = None, # file extension if janis names file
       doc: Union[str, OutputDocumentation] = None,
   ):


.. automethod:: janis.Workflow.output

You are unable to connect an input node directly to an output node, and an output node cannot be referenced as a step input.


.. code-block:: python

    # w.echoStep added to workflow
    w.output("out", source=w.echoStep)


Subclassing Workflow
***********************

Instead of creating inline workflows, it's possible to subclass :class:`janis.Workflow`, implement the required methods which allows a tool to have documentation automatically generated.

Required methods:

.. automethod:: janis.Workflow.id

.. automethod:: janis.Workflow.friendly_name

.. automethod:: janis.Workflow.constructor

Within the ``constructor`` method, you have access to ``self`` to add inputs, steps and outputs.

OPTIONAL:
..........

.. automethod:: janis.Workflow.bind_metadata


Examples
*********


Inline example
................

    The `Echo <https://janis.readthedocs.io/en/latest/tools/unix/echo.html>`_ tool
    has one inputs ``inp`` of type string, and one output ``out``.

.. code-block:: python

    import janis as j
    from janis.unix.tools.echo import Echo

    w = j.WorkflowBuilder("my_workflow")
    w.input("my_input", String)
    echoStep = w.step("echo_step", Echo(inp=w.my_input))
    w.output("out", source=w.echo_step)

    # Will print the CWL, input file and relevant tools to the console
    w.translate("cwl", to_disk=False)  # or "wdl"

Subclass example
.................

.. code-block:: python

    import janis as j
    from janis.unix.tools.echo import Echo

    class MyWorkflow(j.Workflow):

        def id(self):
            return "my_workflow"

        def friendly_name(self):
            return "My workflow"

        def constructor(self):
            self.input("my_input", String)
            echoStep = w.step("echo_step", Echo(inp=self.my_input))
            self.output("out", source=self.echo_step)

        # optional

        def metadata(self):
            self.metadata.author = "Michael Franklin"
            self.metadata.version = "v1.0.0"
            self.metadata.documentation = "A tool that echos the input to standard_out


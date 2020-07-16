Workflow
========

*Dynamically build workflow definitions based on the inputs and hints*

Declaration
############

You might want to dynamically generate your workflow specification based on the inputs and hints.

.. note::

   Dynamic workflows produce a static workflow definition and cannot be used to conditionally run steps based on the outputs of another step. For that, please see the `*conditionals* <conditionals>`_ documentation.

Other important comments:

- ``DynamicWorkflows`` can only be run directly - they cannot be imported from a registry.
- You must override the ``constructor`` argument.
- You could achieve the similar functionality with a ``WorkflowBuilder``, but the DynamicWorkflow class allows `JanisAssistant` to process the inputs and submit to an engine for you.

.. autoclass:: janis.DynamicWorkflow
   :members: constructor, modify_inputs


Example
########

We'll create a workflow that takes an input called ``inp`` that takes a string or an array of strings and we'll create a workflow based on the input.

.. code-block:: python

   class MyFirstDynamicWorkflow(DynamicWorkflow):
       def friendly_name(self):
           return "My first dynamic workflow"

       def id(self) -> str:
           return "MyFirstDynamicWorkflow"

       def constructor(self, inputs: Dict[str, any], hints: Dict[str, str]):

           my_value_to_print = inputs.get("inp")

           if isinstance(my_value_to_print, list):
               # print each element if 'inp' is a list
               self.input("inp_list", Array(String))
               self.step("print", Echo(inp=self.inp_list), scatter="inp")
           else:
               # print only once
               self.input("inp", String)
               self.step("print", Echo(inp=self.inp))

           # Janis will derive the output from the source (Array(String) | String)
           self.output("out", source=self.print)

       def modify_inputs(self, inputs, hints: Dict[str, str]):
           my_value_to_print = inputs.get("inp")
           if isinstance(my_value_to_print, list):
               return {**inputs, "inp_list": my_value_to_print}
           else:
               return inputs



We'll run the following command to generate WDL of the above workflow:

.. code-block:: bash

   janis translate --inputs input.yaml dynamicworkflow.py wdl

Using an ``input.yaml`` with single string, we get the following WDL:

.. code-block:: none

   workflow MyFirstDynamicWorkflow {
     input {
       String inp
     }
     call E.echo as print {
       input:
         inp=inp
     }
     output {
       File out = print.out
     }
   }

Using an ``input.yaml`` with an array of strings, we get the following WDL:

.. code-block:: none

   workflow MyFirstDynamicWorkflow {
     input {
       Array[String] inp_list
     }
     scatter (i in inp_list) {
        call E.echo as print {
         input:
           inp=i
       }
     }
     output {
       Array[File] out = print.out
     }
   }

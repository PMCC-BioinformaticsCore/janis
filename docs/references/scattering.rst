Scattering
==========

*Improving workflow performance with embarrassingly parallel tasks*

Janis support scattering by field when constructing a :class:`janis.Workflow.step()` through the `scatter=Union[str, `:class:`janis.ScatterDescription```]`` parameter.

For example, let's presume you have the tool ``MyTool`` which accepts a single string input on the ``myToolInput`` field.


.. code-block:: python

   w = Workflow("mywf")
   w.input("arrayInp", Array(String))
   w.step("stp", scatter="myToolinput", myToolInput=w.arrayInp)
   # equivalent to
   w.step("stp", scatter=ScatterDescription(fields=["myToolInput"]), myToolInput=w.arrayInp)


Scattering by more than one field
*********************************

Janis only officially supports scattering by 1 field, however there is work being done to support ``dot`` and ``scatter`` methods in WDL through subworkflows.

When it's supported, you will be able to include a ``method=ScatterMethods.{method}`` within the ``ScatterDescription`` constructor. If you scatter by more than field, you WILL need to

WDL
*******
Calls in WDL are explicitly scattered with essentially a for-loop. This calls for
intermediary aliasing of the variable. It attempts to intelligently do this based
on the field name to scatter on where possible.

.. code-block:: wdl

   scatter (b in bams) {
      call G.gatk4sortsam as sortsam {
       input:
         bam=b,
     }
   }



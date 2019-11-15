Scattering
==========

*Improving workflow performance with embarrassingly parallel tasks*

Janis support scattering by field when constructing a :class:`janis.Workflow.step()` through the ``scatter=Union[str,`` :class:`janis.ScatterDescription` ``]`` parameter.


.. autoclass:: janis.ScatterDescription
   :members: __init__

.. autoclass:: janis.ScatterMethods
   :members: dot, cross
   :undoc-members:

Simple scatter
****************

To simply scatter by a single field, you can simple provide the ``scatter="fieldname"`` parameter to the ``janis.Workflow.step()``  method.

For example, let's presume you have the tool ``MyTool`` which accepts a single string input on the ``myToolInput`` field.

.. code-block:: python

   w = Workflow("mywf")
   w.input("arrayInp", Array(String))
   w.step("stp", MyTool(inp1=w.arrayInp), scatter="inp1")
   # equivalent to
   w.step("stp", MyTool(inp1=w.arrayInp), scatter=ScatterDescription(fields=["inp1"]))


Scattering by more than one field
*********************************

Janis supports scattering by multiple fields by the ``dot`` and ``scatter`` methods, you will need to use a :class:`janis.ScatterDescription` and :class:`janis.ScatterMethods`:

Example:

.. code-block:: python

   from janis import ScatterDescription, ScatterMethods
   # OR
   from janis_core import ScatterDescription, ScatterMethods

   w = Workflow("mywf")
   w.input("arrayInp1", Array(String))
   w.input("arrayInp2", Array(String))
   w.step(
     "stp", 
     MyTool(inp1=w.arrayInp1, inp2=w.arrayInp2), 
     scatter=ScatterDescription(fields=["inp1", "inp2"], method=ScatterMethods.dot)
   )

from datetime import date
from inspect import isclass

from janis_core import DataType


def prepare_data_type(dt: DataType):
    dt_name = dt.name()
    secondary = ""

    name = dt.__name__ if isclass(dt) else dt.__class__.__name__

    if dt.secondary_files():
        mapped_secs = "\n".join(f"- ``{s}``" for s in dt.secondary_files())
        secondary = f"""\
Secondary Files
---------------

{mapped_secs}

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__
"""

    return f"""
{dt_name}
{"=" * len(dt_name)}

{dt.doc()}

{secondary}

Quickstart
-----------

.. code-block:: python

   from {dt.__module__} import {name}

   w = WorkflowBuilder("my_workflow")

   w.input("input_{name.lower()}", {name}(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""

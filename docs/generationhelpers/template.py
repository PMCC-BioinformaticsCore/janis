from textwrap import dedent

import ruamel.yaml
from tabulate import tabulate

from janis_assistant.templates import EnvironmentTemplate, get_schema_for_template


def prepare_template(name: str, template: EnvironmentTemplate):
    schema = get_schema_for_template(template)

    tname = name.replace("_", " ").title()

    required_params = [s for s in schema if not s.optional]

    required_cli_params = "".join(
        " \\\n       --" + s.id() + ("" if s.type == bool else " <value>")
        for s in required_params
    )

    required_yaml = {}
    for s in required_params:
        if s.type == bool:
            required_yaml[s.id()] = True
        else:
            required_yaml[s.id()] = "<value>"
    required_yaml_str = (
        ""
        if not required_yaml
        else "".join(
            " " * 5 + s
            for s in ruamel.yaml.dump(
                required_yaml, default_flow_style=False
            ).splitlines(True)
        )
    )

    optional = [(s.id(), s.type, s.default, s.doc or "") for s in schema if s.optional]

    # mapped

    required_section = ""

    if required_params:
        required = [(s.id(), s.type, s.doc or "") for s in required_params]
        requiredf = tabulate(required, ["ID", "Type", "Documentation"], tablefmt="rst")
        required_section = f"""\
**Required**

{requiredf}"""

    optionalf = tabulate(
        optional, ["ID", "Type", "Default", "Documentation"], tablefmt="rst"
    )

    template_doc = dedent(template.__doc__ or "")

    return f"""\
{tname}
{"=" * len(tname)}

Template ID: ``{name}``

{template_doc}

Quickstart
-----------

Take note below how to configure the template. This quickstart only includes the fields you absolutely require. \
This will write a new configuration to ``~/.janis.conf``. See `configuring janis <https://janis.readthedocs.io/en/latest/references/configuration.html>`__ for more information.

.. code-block:: bash

   janis init {name}{required_cli_params}
   
   # or to find out more information
   janis init {name} --help

OR you can insert the following lines into your template:

.. code-block:: yaml

   template:
     id: {name}
{required_yaml_str}


Fields
-------

{required_section}

**Optional**

{optionalf}

"""

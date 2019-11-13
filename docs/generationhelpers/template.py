from tabulate import tabulate

from janis_assistant.templates import EnvironmentTemplate, get_schema_for_template


def prepare_template(name: str, template: EnvironmentTemplate):
    schema = get_schema_for_template(template)

    tname = name.title()

    required = [(s.id(), s.type) for s in schema if not s.optional]
    optional = [(s.id(), s.type, s.default) for s in schema if s.optional]

    # mapped

    requiredf = tabulate(required, ["ID", "Type"], tablefmt="rst")

    optionalf = tabulate(optional, ["ID", "Type", "Default"], tablefmt="rst")

    nl = "\n"

    return f"""\
{tname}
{"=" * len(tname)}

Fields
-------

**Required**

{requiredf}

**Optional**

{optionalf}

"""

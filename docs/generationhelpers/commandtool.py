from textwrap import indent
from typing import List
from tabulate import tabulate


from janis_core import CommandTool, ToolMetadata

from docs.generationhelpers.utils import (
    prepare_quickstart,
    prepare_container_warning_for_commandtool,
)
from .utils import prepare_byline, format_rst_link, get_tool_url


def prepare_commandtool_page(tool: CommandTool, versions: List[str]):
    if not tool:
        return None

    # tool_modules = tool.__module__.split(".") # janis._bioinformatics_.tools.$toolproducer.$toolname.$version

    metadata: ToolMetadata = tool.bind_metadata() or tool.metadata

    if not tool.friendly_name():
        raise Exception(
            f"Tool '{type(tool).__name__}' ({tool.id()}) did not provide the required 'friendly_name' for the docs"
        )

    fn = tool.friendly_name() if tool.friendly_name() else tool.id()
    en = f" ({tool.id()})" if fn != tool.id() else ""
    tn = fn + en

    has_container = tool.container() is not None

    onelinedescription = prepare_byline(
        tool.id(), metadata.short_documentation, metadata.contributors, versions
    )

    formatted_url = (
        format_rst_link(metadata.documentationUrl, metadata.documentationUrl)
        if metadata.documentationUrl
        else "*No URL to the documentation was provided*"
    )

    toolmetadata = [
        t
        for t in [
            ("ID", f"``{tool.id()}``"),
            ("URL", formatted_url),
            ("Versions", ", ".join(versions)),
            ("Container", tool.container()),
            ("Authors", ", ".join(metadata.contributors)),
            ("Citations", metadata.citation),
            ("DOI", metadata.doi) if metadata.doi else None,
            ("Created", str(metadata.dateCreated)),
            ("Updated", str(metadata.dateUpdated)),
        ]
        if t and len(t) == 2
    ]

    input_headers = ["name", "type", "prefix", "position", "documentation"]
    argument_headers = ["value", "prefix", "position", "documentation"]

    required_input_tuples = [
        [i.id(), i.input_type.id(), i.prefix, i.position, i.doc.doc if i.doc else ""]
        for i in tool.inputs()
        if not i.input_type.optional
    ]
    optional_input_tuples = [
        [i.id(), i.input_type.id(), i.prefix, i.position, i.doc.doc if i.doc else ""]
        for i in tool.inputs()
        if i.input_type.optional
    ]

    formatted_args = None
    args = tool.arguments()
    if args:
        argument_tuples = [
            [str(a.value), a.prefix, a.position, a.doc.doc if a.doc else ""]
            for a in tool.arguments()
        ]
        fargs = tabulate(argument_tuples, headers=argument_headers, tablefmt="rst")
        formatted_args = "Arguments\n----------\n\n" + fargs

    formatted_inputs = tabulate(
        required_input_tuples + optional_input_tuples, input_headers, tablefmt="rst"
    )

    formatted_toolversions_array = []
    formatted_toolincludes_array = []
    for v in versions:
        link = get_tool_url(tool.id(), v)
        formatted_toolincludes_array.append(".. include:: " + link)
        if v == tool.version():
            formatted_toolversions_array.append(
                f"- {v} (current)"
            )  # + format_rst_link(v + " (current)", link))
        else:
            formatted_toolversions_array.append(
                "- " + format_rst_link(v, link + ".html")
            )

    output_headers = ["name", "type", "documentation"]
    output_tuples = [[o.id(), o.outtype.id(), o.doc] for o in tool.tool_outputs()]
    formatted_outputs = tabulate(output_tuples, output_headers, tablefmt="rst")

    cwl = tool.translate("cwl", to_console=False, allow_empty_container=True)
    wdl = tool.translate("wdl", to_console=False, allow_empty_container=True)

    tool_prov = ""
    if tool.tool_provider() is None:
        print("Tool :" + tool.id() + " has no company")
    else:
        tool_prov = "." + tool.tool_provider().lower()

    nl = "\n"

    return f"""\
:orphan:

{fn}
{"=" * len(tn)}

{onelinedescription}

{metadata.documentation if metadata.documentation else "No documentation was provided: " + format_rst_link(
    "contribute one", f"https://github.com/PMCC-BioinformaticsCore/janis-{tool.tool_module()}")}

{prepare_container_warning_for_commandtool(tool)}
{prepare_quickstart(tool)}

Information
------------

{nl.join(f":{key}: {value}" for key, value in toolmetadata)}


Outputs
-----------

{formatted_outputs}


Additional configuration (inputs)
---------------------------------

{formatted_inputs}

Workflow Description Language
------------------------------

.. code-block:: text

{indent(wdl, "   ")}

Common Workflow Language
-------------------------

.. code-block:: text

{indent(cwl, "   ")}

"""

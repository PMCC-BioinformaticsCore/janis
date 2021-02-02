from typing import List

from janis_core import CodeTool, ToolMetadata
from tabulate import tabulate

from docs.generationhelpers.utils import (
    prepare_byline,
    format_rst_link,
    get_tool_url,
    prepare_quickstart,
    prepare_container_warning_for_commandtool,
)


def prepare_code_tool_page(tool: CodeTool, versions: List[str]):
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

    input_headers = ["name", "type", "documentation"]

    required_input_tuples = [
        [i.id(), i.intype.id(), i.doc.doc if i.doc else ""]
        for i in tool.tool_inputs()
        if not i.intype.optional
    ]
    optional_input_tuples = [
        [i.id(), i.intype.id(), i.doc.doc if i.doc else ""]
        for i in tool.tool_inputs()
        if i.intype.optional
    ]

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
    """

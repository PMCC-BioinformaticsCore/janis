from typing import List
from tabulate import tabulate


from janis_core import CommandTool, ToolMetadata


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

    onelinedescription = prepare_byline(
        metadata.short_documentation, metadata.contributors, versions
    )

    toolmetadata = [
        t
        for t in [
            ("ID", f"``{tool.id()}``"),
            ("Python", f"``{tool.__module__} import {tool.__class__.__name__}``"),
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

    formatted_url = (
        format_rst_link(metadata.documentationUrl, metadata.documentationUrl)
        if metadata.documentationUrl
        else "*No URL to the documentation was provided*"
    )

    input_headers = ["name", "type", "documentation"]

    required_input_tuples = [
        [i.id(), i.intype.id(), i.doc]
        for i in tool.tool_inputs()
        if not i.intype.optional
    ]
    optional_input_tuples = [
        [i.id(), i.intype.id(), i.doc] for i in tool.tool_inputs() if i.intype.optional
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

{nl.join(f":{key}: {value}" for key, value in toolmetadata)}
:Required inputs:
{(2 * nl).join(f"   - ``{ins.id()}: {ins.input_type.id()}``" for ins in tool.inputs() if (not ins.input_type.optional and ins.default is None))}
:Outputs: 
{(2 * nl).join(f"   - ``{out.id()}: {out.output_type.id()}``" for out in tool.outputs())}

Documentation
-------------

URL: {formatted_url}

{metadata.documentation if metadata.documentation else "No documentation was provided: " + format_rst_link(
    "contribute one", f"https://github.com/PMCC-BioinformaticsCore/janis-{tool.tool_module()}")}

------

Additional configuration (inputs)
---------------------------------

{formatted_inputs}

"""

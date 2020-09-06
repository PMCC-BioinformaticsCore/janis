from textwrap import indent
from typing import List

from janis_core import WorkflowMetadata, Workflow
from requests.utils import requote_uri
from tabulate import tabulate

from docs.generationhelpers.utils import prepare_run_instructions
from .utils import prepare_byline, format_rst_link, get_tool_url, version_html


def generate_pipeline_box(workflow: Workflow, leading_space=""):

    meta: WorkflowMetadata = workflow.bind_metadata() or workflow.metadata

    tag_component = lambda tag: f'<span class="no-select tagelement">{tag}</span>'
    href = f"{workflow.id().lower()}.html"

    tags = "".join(tag_component(t) for t in (meta.keywords or []))
    # date: datetime = meta.dateUpdated or meta.dateCreated

    contributors = meta.contributors or []
    max_display_conts = 5
    if len(contributors) < max_display_conts:
        contributorstr = ", ".join(meta.contributors or ["None"])
    else:
        nothers = len(contributors) - max_display_conts + 2
        contributorstr = (
            ", ".join(meta.contributors[: max_display_conts - 2])
            + f" and {nothers} others"
        )

    return "\n".join(
        leading_space + l
        for l in f"""
<div class="col-6" style="margin: 10px; padding: 20px; border: 1px solid #e3e3e3; border-radius: 5px;">
    <h4 style="margin-bottom: 10px"><a href="{href}">{workflow.friendly_name()}</a></h4>
    {f'<p style="overflow-x: scroll;">{tags}</p>' if tags else ""}
    <p>{meta.short_documentation or "<em>Short documentation required</em>"}</p>
    <p style="margin-bottom: 10px; font-size: 12px; margin-top: 0px;">Contributors: {contributorstr}</p>
    {version_html(workflow.version() or meta.version)}
</div>
""".split(
            "\n"
        )
    )


def prepare_published_pipeline_page(workflow: Workflow, versions: List[str]):
    if not workflow:
        return None

    # tool_modules = tool.__module__.split(".") # janis._bioinformatics_.tools.$toolproducer.$toolname.$version

    metadata: WorkflowMetadata = workflow.bind_metadata() or workflow.metadata

    if not workflow.friendly_name():
        raise Exception(
            f"Tool '{type(workflow).__name__}' ({workflow.id()}) did not provide the required 'friendly_name' for the docs"
        )

    fn = workflow.friendly_name() if workflow.friendly_name() else workflow.id()
    en = f" ({workflow.id()})" if fn != workflow.id() else ""
    tn = fn + en

    onelinedescription = prepare_byline(
        workflow.id(), metadata.short_documentation, metadata.contributors, versions
    )

    citation = "\n\n".join([el for el in [metadata.citation, metadata.doi] if el])

    toolmetadata = [
        ("ID", f"``{workflow.id()}``"),
        ("Versions", ", ".join(str(s) for s in versions[::-1]) if versions else ""),
        ("Authors", ", ".join(metadata.contributors)),
        ("Citations", citation),
        ("Created", str(metadata.dateCreated)),
        ("Updated", str(metadata.dateUpdated)),
    ]

    cwl = workflow.translate("cwl", to_console=False, allow_empty_container=True)[0]
    wdl = workflow.translate("wdl", to_console=False, allow_empty_container=True)[0]

    formatted_url = (
        format_rst_link(metadata.documentationUrl, metadata.documentationUrl)
        if metadata.documentationUrl
        else "*No URL to the documentation was provided*"
    )

    embeddedtoolsraw = {
        f"{s.tool.id()}/{s.tool.version()}": s.tool
        for s in workflow.step_nodes.values()
    }
    embeddedtools = tabulate(
        [
            [tool.friendly_name(), f"``{key}``"]
            for key, tool in embeddedtoolsraw.items()
        ],
        tablefmt="rst",
    )

    input_headers = ["name", "type", "documentation"]

    required_input_tuples = [
        [i.id(), i.intype.id(), i.doc.doc]
        for i in workflow.tool_inputs()
        if not i.intype.optional
    ]
    optional_input_tuples = [
        [i.id(), i.intype.id(), i.doc.doc]
        for i in workflow.tool_inputs()
        if i.intype.optional
    ]

    formatted_inputs = tabulate(
        required_input_tuples + optional_input_tuples, input_headers, tablefmt="rst"
    )

    formatted_toolversions_array = []
    formatted_toolincludes_array = []
    for v in versions:
        link = get_tool_url(workflow.id(), v)
        formatted_toolincludes_array.append(".. include:: " + link)
        if v == workflow.version():
            formatted_toolversions_array.append(
                f"- {v} (current)"
            )  # + format_rst_link(v + " (current)", link))
        else:
            formatted_toolversions_array.append(
                "- " + format_rst_link(v, link + ".html")
            )

    output_headers = ["name", "type", "documentation"]
    output_tuples = [
        [o.id(), o.outtype.id(), o.doc.doc] for o in workflow.tool_outputs()
    ]
    formatted_outputs = tabulate(output_tuples, output_headers, tablefmt="rst")

    tool_prov = ""
    if workflow.tool_provider() is None:
        print("Tool :" + workflow.id() + " has no company")
    else:
        tool_prov = "." + workflow.tool_provider().lower()

    workflow_image = requote_uri(workflow.versioned_id()) + ".dot.png"

    nl = "\n"

    return f"""\
:orphan:

{fn}
{"=" * len(tn)}

{onelinedescription}

{metadata.documentation if metadata.documentation else "No documentation was provided: " + format_rst_link(
    "contribute one", f"https://github.com/PMCC-BioinformaticsCore/janis-{workflow.tool_module()}")}

Quickstart
-----------

{prepare_run_instructions(workflow)}

Outputs
-----------

{formatted_outputs}

Workflow
--------

.. image:: {workflow_image}


Information
------------


{nl.join(f":{key}: {value}" for key, value in toolmetadata)}

Embedded Tools
~~~~~~~~~~~~~~~~~

{embeddedtools}


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

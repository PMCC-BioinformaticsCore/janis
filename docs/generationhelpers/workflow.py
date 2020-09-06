from textwrap import indent

from requests.utils import requote_uri
from tabulate import tabulate
from typing import List

from janis_core import Workflow, WorkflowMetadata
from .utils import (
    prepare_byline,
    format_rst_link,
    get_tool_url,
    prepare_quickstart,
    prepare_container_warning_for_workflow,
)

SHOW_WORKFLOW_IMAGE = False


def prepare_workflow_page(workflow: Workflow, versions: List[str]):
    if not workflow:
        return None

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

    formatted_url = (
        format_rst_link(metadata.documentationUrl, metadata.documentationUrl)
        if metadata.documentationUrl
        else "*No URL to the documentation was provided*"
    )

    toolmetadata = [
        ("ID", f"``{workflow.id()}``"),
        ("URL", formatted_url),
        ("Versions", ", ".join(str(s) for s in versions[::-1]) if versions else ""),
        ("Authors", ", ".join(metadata.contributors)),
        ("Citations", citation),
        ("Created", str(metadata.dateCreated)),
        ("Updated", str(metadata.dateUpdated)),
    ]

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
        [i.id(), i.intype.id(), i.doc.doc if i.doc else ""]
        for i in workflow.tool_inputs()
        if not i.intype.optional
    ]
    optional_input_tuples = [
        [i.id(), i.intype.id(), i.doc.doc if i.doc else ""]
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

    cwl = workflow.translate("cwl", to_console=False, allow_empty_container=True)[0]
    wdl = workflow.translate("wdl", to_console=False, allow_empty_container=True)[0]

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

{prepare_container_warning_for_workflow(workflow)}
{prepare_quickstart(workflow)}

Information
------------

URL: {formatted_url}

{nl.join(f":{key}: {value}" for key, value in toolmetadata)}



Outputs
-----------

{formatted_outputs}


Workflow
--------

.. image:: {workflow_image}

Embedded Tools
***************

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

from tabulate import tabulate
from ruamel.yaml import comments, round_trip_dump

from typing import List, Dict, Union

from janis_core import (
    WorkflowMetadata,
    Workflow,
    Array,
    File,
    DataType,
    String,
    Int,
    Float,
    Boolean,
    InputQualityType,
)
from janis_core.translations import CwlTranslator

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


def prepare_run_instructions(workflow: Workflow, metadata: WorkflowMetadata):

    has_array_of_arrays_inps = True
    bla = any(
        (isinstance(i.intype, Array) and isinstance(i.intype.subtype(), Array))
        for i in workflow.tool_inputs()
    )

    static_input_tuples = [
        [i.id(), i.intype.id(), i.doc.example, i.doc.doc]
        for i in workflow.tool_inputs()
        if i.doc.quality == InputQualityType.static
    ]

    reference_information = ""

    if len(static_input_tuples) > 0:
        static_input_headers = ["Name", "Type", "Example", "Description"]
        reference_information = tabulate(
            static_input_tuples, headers=static_input_headers, tablefmt="rst"
        )

    overrides = metadata.sample_input_overrides or {}
    inps = {}
    for i in workflow.tool_inputs():
        if i.intype.optional or i.default:
            continue
        inps[i.id()] = (
            overrides.get(i.id())
            if i.id() in overrides
            else prepare_default_for_type(i.id(), i.intype)
        )
        # if i.doc:
        #     inps.yaml_set_comment_before_after_key(i.id(), i.doc)

    if has_array_of_arrays_inps:
        return prepare_run_instructions_input_file(
            workflow, inps, reference_information
        )
    else:
        return prepare_run_instructions_cli(workflow, inps, reference_information)


def prepare_default_for_type(identifier: str, t: DataType, idx=None):
    if isinstance(t, Array):
        st = t.subtype()
        return [
            prepare_default_for_type(identifier, st, 0),
            prepare_default_for_type(identifier, st, 1),
        ]
    elif isinstance(t, File):
        return (
            identifier
            + ("_" + str(idx) if idx is not None else "")
            + str(t.extension or "")
        )

    elif isinstance(t, String):
        return "<value>"
    elif isinstance(t, Int):
        return 0
    elif isinstance(t, Float):
        return 0.0
    elif isinstance(t, Boolean):
        return True
    return None


def prepare_run_instructions_input_file(
    workflow: Workflow, ins, reference_information: str
):
    yaml_inp = CwlTranslator.stringify_translated_inputs(ins)
    indented = "".join(" " * 7 + s for s in yaml_inp.splitlines(True))

    return f"""\
1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.

{reference_information}

4. Generate an inputs file for {workflow.id()}:

.. code-block:: bash
   
   janis inputs {workflow.id()} > inputs.yaml

**inputs.yaml**

.. code-block:: yaml

{indented}

5. Run the {workflow.id()} pipeline with:

.. code-block:: bash

   janis run [...workflow options] --inputs inputs.yaml {workflow.id()}

"""


def prepare_run_instructions_cli(workflow: Workflow, ins: dict):

    sp = 7 * " "
    cli = []
    for k, v in ins.items():
        vv = v
        if isinstance(v, bool):
            cli.append("--" + k)
        else:
            if isinstance(v, list):
                vv = " ".join(str(v))
            cli.append("--" + k + " " + str(vv))

    clis = (" \\\n" + sp).join(cli)

    return f"""\
In your terminal, you can run the workflow with the following command:

.. code-block:: bash

   janis run [...workflow options] {workflow.id()} \\
       {clis}

See below for documentation for each input.
"""


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
        metadata.short_documentation, metadata.contributors, versions
    )

    citation = "\n\n".join([el for el in [metadata.citation, metadata.doi] if el])

    toolmetadata = [
        ("ID", f"``{workflow.id()}``"),
        ("Python", f"``{workflow.__module__} import {workflow.__class__.__name__}``"),
        ("Versions", ", ".join(str(s) for s in versions[::-1]) if versions else ""),
        ("Authors", ", ".join(metadata.contributors)),
        ("Citations", citation),
        ("Created", str(metadata.dateCreated)),
        ("Updated", str(metadata.dateUpdated)),
    ]

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

    workflow_image = """
Workflow
--------

.. raw:: html

   <script src="https://cdnjs.cloudflare.com/ajax/libs/vue/2.6.10/vue.min.js"></script>
   <script src="https://unpkg.com/vue-cwl/dist/index.js"></script>
   <div id="vue" style="width: 800px; height: 500px; border-radius: 5px; overflow: hidden;">
          <cwl cwl-url="https://unpkg.com/cwl-svg@2.1.5/cwl-samples/fastqc.json"></cwl>
   </div>
   <script>
   new Vue({
       el: '#vue',
       components: {
           cwl: vueCwl.default
       }
   });
   </script>
    """

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

{prepare_run_instructions(workflow, metadata)}

Outputs
-----------

{formatted_outputs}


Information
------------

{nl.join(f":{key}: {value}" for key, value in toolmetadata)}

Embedded Tools
~~~~~~~~~~~~~~~~~

{embeddedtools}


Additional configuration (inputs)
---------------------------------

{formatted_inputs}
"""

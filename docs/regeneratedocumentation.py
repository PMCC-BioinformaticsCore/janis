"""
This class regenerates all the tool definitions

It's a bit of a random collection of things that should be refactored:

    - NestedDictionary (store values in nested structure (with root key))
    - RST Helpers
"""
from requests.utils import requote_uri
from distutils.version import StrictVersion
from inspect import isfunction, ismodule, isabstract, isclass
from janis_runner.templates import (
    templates,
    get_schema_for_template,
    EnvironmentTemplate,
)

import janis
from constants import PROJECT_ROOT_DIR
from janis_core import (
    Array,
    DataType,
    Workflow,
    CommandTool,
    Tool,
    Metadata,
    Logger,
    JanisShed,
    WorkflowMetadata,
    ToolTypes,
)

# Import modules here so that the tool registry knows about them

# Output settings
docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"
dt_dir = docs_dir + "datatypes/"
templates_dir = docs_dir + "templates/"
pipelines_dir = docs_dir + "pipelines/"

modules = [janis.bioinformatics, janis.unix]

###### Shouldn't need to touch below this line #####

import os
from typing import List, Set, Type, Tuple, Dict

import tabulate
from datetime import date, datetime


class NestedDictionaryTypeException(Exception):
    def __init__(self, key, error_type, keys=None):
        keychain = (" for keys: " + ".".join(keys)) if keys else ""
        super(NestedDictionaryTypeException, self).__init__(
            f"Incorrect type '{error_type}' for key '{key}'{keychain}"
        )
        self.keys = keys
        self.key = key
        self.error_type = error_type


def format_rst_link(text, link):
    return f"`{text} <{link}>`_"


def get_tool_url(toolname, version):
    return requote_uri(toolname + "_" + version).lower()


def prepare_tool(tool: Tool, toolversions: List[str], isorphan: bool):
    # Stuff to list on the documentation page:
    #   - Versions of tools
    #   - Generated command
    #   - Cool if it grouped the tools by vendor
    #   -

    if not tool:
        return None

    # tool_modules = tool.__module__.split(".") # janis._bioinformatics_.tools.$toolproducer.$toolname.$version

    metadata = tool.bind_metadata() or tool.metadata

    if not tool.friendly_name():
        raise Exception(
            f"Tool '{type(tool).__name__}' ({tool.id()}) did not provide the required 'friendly_name' for the docs"
        )

    fn = tool.friendly_name() if tool.friendly_name() else tool.id()
    en = f" ({tool.id()})" if fn != tool.id() else ""
    tn = fn + en

    formatted_url = (
        format_rst_link(metadata.documentationUrl, metadata.documentationUrl)
        if metadata.documentationUrl
        else "*No URL to the documentation was provided*"
    )

    input_headers = ["name", "type", "prefix", "position", "documentation"]

    required_input_tuples = [
        [i.id(), i.input_type.id(), i.prefix, i.position, i.doc]
        for i in tool.inputs()
        if not i.input_type.optional
    ]
    optional_input_tuples = [
        [i.id(), i.input_type.id(), i.prefix, i.position, i.doc]
        for i in tool.inputs()
        if i.input_type.optional
    ]

    formatted_required_inputs = tabulate.tabulate(
        required_input_tuples, input_headers, tablefmt="rst"
    )
    formatted_optional_inputs = tabulate.tabulate(
        optional_input_tuples, input_headers, tablefmt="rst"
    )

    versiontext = ""

    formatted_toolversions_array = []
    formatted_toolincludes_array = []
    for v in toolversions:
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

    formatted_toolincludes = "\n".join(formatted_toolincludes_array)
    if len(formatted_toolincludes_array) > 1:

        versiontext = "Versions\n*********\n\n" + "\n".join(
            formatted_toolversions_array
        )

    output_headers = ["name", "type", "documentation"]
    output_tuples = [[o.id(), o.output_type.id(), o.doc] for o in tool.outputs()]
    formatted_outputs = tabulate.tabulate(output_tuples, output_headers, tablefmt="rst")

    container_tag = ""
    if isinstance(tool, CommandTool):
        container_tag = f"Container: ``{tool.container()}``"

    tool_prov = ""
    if tool.tool_provider() is None:
        print("Tool :" + tool.id() + " has no company")
    else:
        tool_prov = "." + tool.tool_provider().lower()

    return f"""\
{':orphan:' if isorphan else ''}
{formatted_toolincludes if not isorphan else ''}

{fn}
{"=" * len(tn)}

Description
-------------

Tool identifier: ``{tool.id()}``

Tool path: ``{tool.__module__} import {tool.__class__.__name__}``

Version: {tool.version()}

{container_tag}

{versiontext}

Documentation
-------------

URL
******

{formatted_url}

Tool documentation
******************

{metadata.documentation if metadata.documentation else "*No documentation was provided: " + format_rst_link(
        "contribute one", f"https://github.com/PMCC-BioinformaticsCore/janis-{tool.tool_module()}") + "*"}

Outputs
-------
{formatted_outputs}

Inputs
------


Required inputs
***************

{formatted_required_inputs}

Optional inputs
***************

{formatted_optional_inputs}


Metadata
********

Contributors: {", ".join(metadata.contributors) if metadata.contributors else "**Unknown**"}


*{fn} was last updated on {metadata.dateUpdated if metadata.dateUpdated else "**Unknown**"}*.
*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""


def prepare_workflow_page(workflow: Workflow, versions: List[str]):
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
        ("Versions", "\n".join(versions)),
        ("Contributors", "\n".join(metadata.contributors)),
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
    embeddedtools = tabulate.tabulate(
        [
            [tool.friendly_name(), f"``{key}``"]
            for key, tool in embeddedtoolsraw.items()
        ],
        tablefmt="rst",
    )

    input_headers = ["name", "type", "documentation"]

    required_input_tuples = [
        [i.id(), i.input_type.id(), i.doc]
        for i in workflow.inputs()
        if not i.input_type.optional
    ]
    optional_input_tuples = [
        [i.id(), i.input_type.id(), i.doc]
        for i in workflow.inputs()
        if i.input_type.optional
    ]

    formatted_inputs = tabulate.tabulate(
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
    output_tuples = [[o.id(), o.output_type.id(), o.doc] for o in workflow.outputs()]
    formatted_outputs = tabulate.tabulate(output_tuples, output_headers, tablefmt="rst")

    tool_prov = ""
    if workflow.tool_provider() is None:
        print("Tool :" + workflow.id() + " has no company")
    else:
        tool_prov = "." + workflow.tool_provider().lower()

    workflow_image = (
        "."
        or """
Workflow
--------

.. raw:: html

   <script src="https://cdnjs.cloudflare.com/ajax/libs/vue/2.6.10/vue.min.js"></script>
   <script src="https://unpkg.com/vue-cwl/dist/index.js"></script>
   <div id="vue" style="width: 800px; height: 500px; border-radius: 5px; overflow: hidden;">
          <cwl cwl-url="https://unpkg.com/cwl-svg@2.1.5/cwl-samples/fastqc.json"></cwl>
   </div>
   <script>
   new Vue({{
       el: '#vue',
       components: {{
           cwl: vueCwl.default
       }}
   }});
   </script>
    """
    )

    nl = "\n"

    return f"""\
:orphan:

{fn}
{"=" * len(tn)}

{onelinedescription}

{nl.join(f":{key}: {value}" for key, value in toolmetadata)}
:Required inputs:
{(2 * nl).join(f"   - ``{ins.id()}: {ins.datatype.id()}``" for ins in workflow.input_nodes.values() if (not ins.datatype.optional and ins.default is None))}
:Outputs: 
{(2 * nl).join(f"   - ``{out.id()}: {out.datatype.id()}``" for out in workflow.output_nodes.values())}

Documentation
-------------

URL: {formatted_url}

{metadata.documentation if metadata.documentation else "No documentation was provided: " + format_rst_link(
    "contribute one", f"https://github.com/PMCC-BioinformaticsCore/janis-{workflow.tool_module()}")}

Embedded Tools
***************

{embeddedtools}

------

Inputs
------

{formatted_inputs}

{workflow_image}
"""


def prepare_byline(
    shortdocumentation: str, contributors: List[str], versions: List[str]
):
    short_prepared = (
        shortdocumentation[:-1]
        if shortdocumentation and shortdocumentation[-1] == "."
        else shortdocumentation
    )

    contributors, versions = contributors or [], versions or []
    contributorstr = (
        f"{len(contributors)} contributor{'s' if len(contributors) != 1 else ''}"
    )
    versionstr = f"{len(versions)} version{'s' if len(versions) != 1 else ''}"
    return f"{short_prepared} · {contributorstr} · {versionstr}"


def prepare_data_type(dt: DataType):
    dt_name = dt.name()
    secondary = ""

    if dt.secondary_files():
        secondary = "Secondary files: " + ", ".join(
            f"``{s}``" for s in dt.secondary_files()
        )

    return f"""
{dt_name}
{"=" * len(dt_name)}

{secondary}

Documentation
-------------

{dt.doc()}

*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""


def nested_keys_append_with_root(d, keys: List[str], value, root_key):
    if len(keys) == 0:
        if root_key in d:
            d[root_key].append(value)
        else:
            d[root_key] = [value]
        return d
    try:
        key = keys[0]
        if key in d:
            if not isinstance(d[key], dict):
                raise NestedDictionaryTypeException(
                    key=key, error_type=type(d[key]), keys=keys
                )
            nested_keys_append_with_root(d[key], keys[1:], value, root_key=root_key)
            return d
        else:
            d[key] = nested_keys_append_with_root(
                {}, keys[1:], value, root_key=root_key
            )
    except NestedDictionaryTypeException as de:
        raise NestedDictionaryTypeException(
            key=de.key, error_type=de.error_type, keys=keys
        )

    return d


def nested_keys_add(d, keys: List[str], value) -> bool:
    if len(keys) == 0:
        raise Exception(
            f"Couldn't add {value} to nested dictionary as no keys were provided"
        )
    key = keys[0]
    if len(keys) == 1:
        if key in d:
            return False
        d[key] = value
        return True
    try:
        if key not in d:
            d[key] = {}
        elif not isinstance(d[key], dict):
            raise NestedDictionaryTypeException(
                key=key, error_type=type(d[key]), keys=keys
            )
        return nested_keys_add(d[key], keys[1:], value)
    except NestedDictionaryTypeException as de:
        raise NestedDictionaryTypeException(
            key=de.key, error_type=de.error_type, keys=keys
        )


def get_toc(title, intro_text, subpages, caption="Contents", max_depth=1):
    prepared_subpages = "\n".join(
        "   " + m.lower() for m in sorted(subpages, key=lambda l: l.lower())
    )
    return f"""
{title.replace('{title}', title)}
{"=" * len(title)}

{intro_text}

.. toctree::
   :maxdepth: {max_depth}
   :caption: {caption}:

{prepared_subpages}

*This page was auto-generated on {date.today().strftime(
        "%d/%m/%Y")}. Please do not directly alter the contents of this page.*
"""


def sort_tool_versions(versions: List[str]) -> List[str]:
    try:
        return sorted(versions, key=StrictVersion, reverse=True)
    except:
        return sorted(versions, reverse=True)


def prepare_all_tools():
    JanisShed.hydrate(modules=[janis.unix, janis.bioinformatics])

    data_types = JanisShed.get_all_datatypes()
    tools = {
        ts[0].id(): {t.version(): t for t in ts} for ts in JanisShed.get_all_tools()
    }

    Logger.info(f"Preparing documentation for {len(tools)} tools")
    Logger.info(f"Preparing documentation for {len(data_types)} data_types")

    tool_module_index = {}
    dt_module_index = {}
    ROOT_KEY = "root"

    for toolname, toolsbyversion in tools.items():
        # tool = tool_vs[0][0]()
        tool_versions = sort_tool_versions(list(toolsbyversion.keys()))
        default_version = tool_versions[0]
        Logger.log(
            f"Preparing {toolname}, found {len(tool_versions)} version[s] ({','.join(tool_versions)})"
        )

        defaulttool = toolsbyversion[default_version]
        tool_path_components = list(
            filter(
                lambda a: bool(a),
                [defaulttool.tool_module(), defaulttool.tool_provider()],
            )
        )

        # (toolURL, tool, isPrimary)
        toolurl_to_tool = [(toolname.lower(), defaulttool, True)] + [
            (get_tool_url(toolname, v), toolsbyversion[v], False) for v in tool_versions
        ]

        path_components = "/".join(tool_path_components)
        output_dir = f"{tools_dir}/{path_components}/".lower()
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for (toolurl, tool, isprimary) in toolurl_to_tool:
            output_str = prepare_tool(tool, tool_versions, not isprimary)
            output_filename = output_dir + toolurl + ".rst"
            with open(output_filename, "w+") as tool_file:
                tool_file.write(output_str)

        nested_keys_append_with_root(
            tool_module_index, tool_path_components, toolname, root_key=ROOT_KEY
        )

        Logger.log("Prepared " + toolname)

    for d in data_types:
        # tool = tool_vs[0][0]()
        if issubclass(d, Array):
            Logger.info("Skipping Array DataType")
            continue
        try:
            dt = d()
        except:
            print(d.__name__ + " failed to instantiate")
            continue
        did = dt.name().lower()
        Logger.log("Preparing " + dt.name())
        output_str = prepare_data_type(dt)

        dt_path_components = []
        # dt_path_components = list(filter(
        #     lambda a: bool(a),
        #     [, tool.tool_provider()]
        # ))

        path_components = "/".join(dt_path_components)
        output_dir = f"{dt_dir}{path_components}/"
        output_filename = output_dir + did + ".rst"

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        nested_keys_append_with_root(
            dt_module_index, dt_path_components, did, root_key=ROOT_KEY
        )

        with open(output_filename, "w+") as dt_file:
            dt_file.write(output_str)

        Logger.log("Prepared " + did)

    def prepare_modules_in_index(contents, title, dir, max_depth=1):
        module_filename = dir + "/index.rst"
        module_tools = sorted(set(contents[ROOT_KEY] if ROOT_KEY in contents else []))
        submodule_keys = sorted(m for m in contents.keys() if m != ROOT_KEY)
        indexed_submodules_tools = [m.lower() + "/index" for m in submodule_keys]

        with open(module_filename, "w+") as module_file:
            module_file.write(
                get_toc(
                    title=title,
                    intro_text="Automatically generated index page for {title}",
                    subpages=indexed_submodules_tools + module_tools,
                    max_depth=max_depth,
                )
            )

        for submodule in submodule_keys:
            prepare_modules_in_index(
                contents=contents[submodule], title=submodule, dir=f"{dir}/{submodule}/"
            )

    prepare_modules_in_index(tool_module_index, title="Tools", dir=tools_dir)
    prepare_modules_in_index(
        dt_module_index, title="Data Types", dir=dt_dir, max_depth=1
    )


def prepare_runner_templates():
    """
    Templates can probably 
    :return:
    """

    os.makedirs(templates_dir, exist_ok=True)

    for tkey, template in templates.items():
        with open(os.path.join(templates_dir, tkey + ".rst"), "w+") as f:
            f.write(prepare_template(tkey, template))

    with open(os.path.join(templates_dir, "index.rst"), "w+") as f:
        f.write(
            get_toc(
                "Templates",
                intro_text="List of templates for ``janis-runner``",
                subpages=list(templates.keys()),
            )
        )


def prepare_template(name: str, template: EnvironmentTemplate):

    schema = get_schema_for_template(template)

    tname = name.title()

    required = [(s.id(), s.type) for s in schema if not s.optional]
    optional = [(s.id(), s.type, s.default) for s in schema if s.optional]

    # mapped

    requiredf = tabulate.tabulate(required, ["ID", "Type"], tablefmt="rst")

    optionalf = tabulate.tabulate(optional, ["ID", "Type", "Default"], tablefmt="rst")

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


def generate_pipelines_page():
    print("Generating pipelines page...")

    modules = [janis.pipelines]
    workflows: List[Workflow] = []

    for module in modules:

        workflows.extend(
            cls()
            for n, cls in list(module.__dict__.items())
            if not n.startswith("__")
            and type(cls) != type
            and isclass(cls)
            and hasattr(cls, "type")
            and cls.type() == ToolTypes.Workflow
        )

    wf_strings = "\n".join(
        generate_pipeline_box(w, leading_space=" " * 8) for w in workflows
    )

    page = f"""\
Pipelines
=========

.. raw:: html

    <div id="box" style="display: flex; flex-flow: wrap;">
{wf_strings}
    </div>
"""

    os.makedirs(pipelines_dir, exist_ok=True)

    with open(os.path.join(pipelines_dir, "index.rst"), "w+") as f:
        f.write(page)

    # Write all the pages
    for w in workflows:
        toolstr = prepare_workflow_page(w, [w.version()])
        with open(os.path.join(pipelines_dir, w.id().lower() + ".rst"), "w+") as f:
            f.write(toolstr)


def generate_pipeline_box(workflow: Workflow, leading_space=""):

    meta: WorkflowMetadata = workflow.bind_metadata() or workflow.metadata

    tag_component = lambda tag: f'<span class="no-select tagelement">{tag}</span>'
    href = f"{workflow.id().lower()}.html"

    version_component = (
        lambda version: f'<p style="margin-bottom: 10px"><a class="version-button" href="{href}">'
        f'    Version <b>{version[1:] if (version and version.startswith("v")) else version}</b>'
        f"</a></p>"
    )

    tags = "".join(tag_component(t) for t in (meta.keywords or []))
    date: datetime = meta.dateUpdated or meta.dateCreated

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
    {f'<p>{tags}</p>' if tags else ""}
    <p>{meta.short_documentation or "<em>Short documentation required</em>"}</p>
    {version_component(workflow.version() or meta.version)}
    <p style="margin-bottom: 0px; font-size: 12px">Contributors: {contributorstr}</p>
</div>""".split(
            "\n"
        )
    )


if __name__ == "__main__":
    # prepare_all_tools()
    # prepare_runner_templates()
    generate_pipelines_page()

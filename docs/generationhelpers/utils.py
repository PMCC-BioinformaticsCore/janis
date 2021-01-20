from textwrap import indent

from janis_core.workflow.workflow import InputNode
from tabulate import tabulate
from typing import List, Dict
from datetime import date
from distutils.version import StrictVersion

from janis_core import (
    Tool,
    Metadata,
    DataType,
    Array,
    File,
    String,
    Int,
    Float,
    Boolean,
    InputQualityType,
    CommandTool,
    Workflow,
    WorkflowBase,
)
from janis_core.translations import CwlTranslator
from requests.utils import requote_uri


class TocObject:
    def __init__(self, title, description, url):
        self.title = title
        self.description = description
        self.url = url


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
    if not (toolname and version):
        return ""
    return requote_uri(toolname + "_" + str(version)).lower()


def prepare_byline(
    toolid, shortdocumentation: str, contributors: List[str], versions: List[str]
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
    return (
        f"``{toolid}``"
        + " · *"
        + " · ".join(t for t in [short_prepared, contributorstr, versionstr] if t)
        + "*"
    )


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


def get_tool_toc(
    alltoolsmap: Dict[str, Dict[str, Tool]],
    title,
    intro_text,
    subpages,
    tools,
    caption="Contents",
    max_depth=1,
):

    pd = " " * 5
    mappedtools = "\n".join(
        "\n".join(pd + r for r in get_tool_row(alltoolsmap[tool]).split("\n"))
        for tool in tools
    )

    mappedmodules = "\n\n".join(
        f'     <li><a href="{m.lower()}/index.html">{m}</a></li>'
        for m in sorted(subpages)
    )

    return f"""
:orphan:

{title.replace('{title}', title)}
{"=" * len(title)}

{intro_text}

.. raw:: html

   <ul>
{mappedmodules}
   </ul>
{mappedtools}

"""


def prepare_container_warning_for_commandtool(tool: CommandTool):
    if tool.container() is not None and len(tool.container()) > 0:
        return ""

    return f"""\
.. warning::

   {tool.friendly_name()} did not include a container. You can provide one through the command line by including
   the following instruction:

   .. code-block:: bash

      janis run --container-override '{tool.id()}=<organisation/container:version>' {tool.id()}
    """


def prepare_container_warning_for_workflow(tool: Workflow):
    def recursive_find_tools_without_container(wf: Workflow):
        tools = {}
        for s in wf.step_nodes.values():
            stool = s.tool
            if isinstance(stool, Workflow):
                tools.update(recursive_find_tools_without_container(stool))
            else:
                if not stool.container():
                    tools[stool.id()] = stool

        return tools

    tools_without_containers = recursive_find_tools_without_container(tool)
    if not tools_without_containers:
        return ""

    intro = (
        "A tool in" if len(tools_without_containers) == 1 else "Some of the tools in"
    )
    merged_command = ", ".join(
        f"{t.id()}=<organisation/container:version>"
        for t in tools_without_containers.values()
    )

    return f"""\
.. warning::

   {intro} {tool.friendly_name()} did not include a container. You could provide them through
    the command line by including the following instruction:

   .. code-block:: bash

      janis run --container-override '{merged_command}' {tool.id()}
    """


def get_tool_row(tools: Dict[str, Tool]):
    distincted = {t.version(): t for t in tools.values()}
    versions = sort_tool_versions(list(distincted.keys()))
    latestversion = versions[0]

    tool = distincted[latestversion]
    meta: Metadata = tool.bind_metadata() or tool.metadata
    sd = meta.short_documentation
    sdstr = f'<p style="color: black; margin-bottom: 10px">{sd}' if sd else ""

    href = tool.id().lower() + ".html"
    return f"""\
<a href="{href}">
  <p style="margin-bottom: 5px"><b>{tool.friendly_name()}</b> <span style="margin-left: 10px; color: darkgray">{tool.id()}</span></p>
  {sdstr}
  <p><span style="margin-right: 10px; color: darkgray">({len(versions)} versions)</span>{version_html(latestversion, href=href)}</p>
</a>
<hr />
    """


def sort_tool_versions(versions: List[str]) -> List[str]:
    try:
        return sorted(versions, key=StrictVersion, reverse=True)
    except:
        return sorted(versions, reverse=True)


def version_html(version, href=None):
    fhref = f'href="{href}"' if href else ""
    return f"""\
<a class="version-button" {fhref} style="margin-bottom: 10px">
  v<b>{version[1:] if (version and version.startswith("v")) else version}</b>
</a>"""


def prepare_quickstart(tool: Tool):
    required_python_input_map = "\n".join(
        " " * 15 + i.id() + "=None,"
        for i in tool.tool_inputs()
        if not i.intype.optional
    )

    python_step_name = tool.id().lower() + "_step"
    output_python_code = "\n".join(
        " " * 7 + f'wf.output("{o.id()}", source={python_step_name}.{o.id()})'
        for o in tool.tool_outputs()
    )
    python_codeblock = f"""\
    .. code-block:: python

       from {tool.__module__} import {tool.__class__.__name__}

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "{python_step_name}",
           {tool.__class__.__name__}(
{required_python_input_map}
           )
       )
{output_python_code}
    """

    return f"""\
Quickstart
-----------

{python_codeblock}

*OR*

{prepare_run_instructions(tool)}

"""


def prepare_source(source):
    if isinstance(source, list):
        return ", ".join(prepare_source(s) for s in source)
    elif isinstance(source, str):
        return source
    elif isinstance(source, dict):
        return "\n".join(f"* {k}: {prepare_source(v)}" for k, v in source.items())
    return str(source)


def prepare_run_instructions(tool: Tool):
    metadata = tool.bind_metadata() or tool.metadata
    has_array_of_arrays_inps = True
    bla = any(
        (isinstance(i.intype, Array) and isinstance(i.intype.subtype(), Array))
        for i in tool.tool_inputs()
    )

    static_input_tuples = [
        [i.id(), i.intype.id(), prepare_source(i.doc.source)]
        for i in tool.tool_inputs()
        if i.doc.quality == InputQualityType.static and i.doc.source is not None
    ]

    reference_information = ""

    if len(static_input_tuples) > 0:
        static_input_headers = ["Name", "Type", "Source"]
        reference_information = tabulate(
            static_input_tuples, headers=static_input_headers, tablefmt="rst"
        )

    # overrides = metadata.sample_input_overrides or {}
    user_inps = {}
    other_inps = {}

    for i in tool.tool_inputs():
        if i.intype.optional or i.default is not None:
            continue

        val = i.doc.example or prepare_default_for_type(i.id(), i.intype)
        if i.doc and i.doc.quality and i.doc.quality != InputQualityType.user:
            other_inps[i.id()] = val
        else:
            # catch None and InputQualityType.user
            user_inps[i.id()] = val

    if has_array_of_arrays_inps:
        return prepare_run_instructions_input_file(
            tool, user_inps, other_inps, reference_information
        )
    else:
        return prepare_run_instructions_cli(
            tool, user_inps, other_inps, reference_information
        )


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
    tool: Tool, user_inps: dict, other_inps: dict, reference_information: str
):
    yaml_user_inps = CwlTranslator.stringify_translated_inputs(user_inps)
    yaml_other_inps = CwlTranslator.stringify_translated_inputs(other_inps)
    indented_user = "".join(" " * 7 + s for s in yaml_user_inps.splitlines(True))
    indented_other = "".join(" " * 7 + s for s in yaml_other_inps.splitlines(True))

    not_localising_secondary_warning = ""
    if isinstance(tool, WorkflowBase):
        inputs_that_arent_localising_secondary_files = [
            t.id() for t in tool.tool_inputs() if t.doc.skip_sourcing_secondary_files
        ]
        if len(inputs_that_arent_localising_secondary_files) > 0:
            not_localising_secondary_warning = f"""\
.. warning::

   The secondary files for the inputs '{"', '".join(inputs_that_arent_localising_secondary_files)}' will not automatically \
   localise using janis prepare and are built just after download. Please note this can take a few hours to build \
   before the pipeline runs. 
"""

    has_static = len(other_inps) > 0

    tb = " " * 4
    run_args = ["janis run [...run options]", tb + "--inputs inputs.yaml"]

    static_generation = (
        ""
        if not has_static
        else f"""\
   # static inputs
   janis inputs --static {tool.id()} > static.yaml"""
    )
    static_yaml = (
        ""
        if not has_static
        else f"""\
**static.yaml**

.. code-block:: yaml

{indented_other}"""
    )
    if has_static:
        run_args.append(tb + "--inputs static.yaml")

    if isinstance(tool, CommandTool) and not tool.container():
        run_args.append(
            tb + f"--container-override '{tool.id()}=<organisation/container:version>'"
        )

    run_args.append(tb + tool.id())
    run_statement = " \\\n".join(" " * 3 + el for el in run_args)

    if reference_information:
        reference_information = f"The following inputs have a suggested source. Using janis prepare with the relevant \
        ``--source-hint`` will automatically download these files. See `below <#additional-configuration-inputs>`_ for \
        more information about inputs for {tool.id()}.\n{reference_information}"

    return f"""\
1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user {'and static ' if has_static else ''}input files for {tool.id()}:

.. code-block:: bash

   # user inputs
   janis inputs {"--user " if has_static else ""}{tool.id()} > inputs.yaml

{static_generation}

**inputs.yaml**

.. code-block:: yaml

{indented_user}

{static_yaml}

5. Run {tool.id()} with:

.. code-block:: bash

{run_statement}

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this {tool.type()}. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \\
          --inputs inputs.yaml \\
          --output-dir $OUTPUT_DIR \\
          {tool.id()}

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh

{indent(not_localising_secondary_warning, "   ")}


{indent(reference_information, '   ')}


"""


def prepare_run_instructions_cli(
    workflow: Tool, user_inps: dict, other_inps: dict, reference_information: str
):
    sp = 7 * " "
    cli = []
    iters = [*list(user_inps.items()), *list(other_inps.items())]
    for k, v in iters:
        vv = v
        if isinstance(v, bool):
            cli.append("--" + k)
        else:
            if isinstance(v, list):
                vv = " ".join(str(v))
            cli.append("--" + k + " " + str(vv))

    clis = (" \\\n" + sp).join(cli)

    return f"""\
1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.

{reference_information}

4. In your terminal, you can run the workflow with the following command:

.. code-block:: bash

   janis run [...workflow options] {workflow.id()} \\
       {clis}

See below for documentation for each input.
"""

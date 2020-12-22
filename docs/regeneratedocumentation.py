#!/usr/bin/env python3

"""
This class regenerates all the tool definitions

It's a bit of a random collection of things that should be refactored:

    - NestedDictionary (store values in nested structure (with root key))
    - RST Helpers
"""
from inspect import isfunction, ismodule, isabstract, isclass
from janis_assistant.templates import (
    get_all_templates,
    get_schema_for_template,
    EnvironmentTemplate,
)

from shutil import rmtree

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
    ToolType,
    WorkflowBase,
)

import janis_unix, janis_bioinformatics

# Import modules here so that the tool toolbox knows about them

# Output settings
from docs.generationhelpers.commandtool import prepare_commandtool_page
from docs.generationhelpers.codetool import prepare_code_tool_page
from docs.generationhelpers.datatype import prepare_data_type
from docs.generationhelpers.pipelines import (
    generate_pipeline_box,
    prepare_published_pipeline_page,
)
from docs.generationhelpers.template import prepare_template
from docs.generationhelpers.utils import (
    sort_tool_versions,
    get_tool_url,
    nested_keys_append_with_root,
    get_toc,
    get_tool_toc,
)
from docs.generationhelpers.workflow import prepare_workflow_page

docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"
dt_dir = docs_dir + "datatypes/"
templates_dir = docs_dir + "templates/"
pipelines_dir = docs_dir + "pipelines/"

modules = [janis_bioinformatics, janis_unix]


tool_module_information = {
    "bioinformatics": {
        "name": "Bioinformatics tools",
        "description": "A collection of tools for bioinformatics",
    },
    "unix": {
        "name": "Unix tools",
        "description": "A collection of tools available on unix systems",
    },
}

###### Shouldn't need to touch below this line #####

import os
import tabulate
from datetime import date, datetime
from typing import List, Set, Type, Tuple, Dict
import traceback


def prepare_tool(
    tool: Tool,
    toolversions: List[str],
    isorphan: bool,
    is_published_pipeline: bool = False,
):
    # Stuff to list on the documentation page:
    #   - Versions of tools
    #   - Generated command
    #   - Cool if it grouped the tools by vendor
    #   -

    if not tool:
        return None
    try:
        if is_published_pipeline:
            return ""
        if tool.type() == ToolType.CommandTool:
            return prepare_commandtool_page(tool, toolversions)
        elif tool.type() == ToolType.Workflow:
            return prepare_workflow_page(tool, toolversions)
        elif tool.type() == ToolType.CodeTool:
            return prepare_code_tool_page(tool, toolversions)
    except Exception as e:
        traceback.print_exc()
        Logger.critical(
            "Couldn't generate documentation for " + tool.id() + " " + str(e)
        )


def prepare_all_tools():
    JanisShed.hydrate(modules=[janis_unix, janis_bioinformatics])

    data_types = JanisShed.get_all_datatypes()
    tools = {
        ts[0].id(): {t.version(): t for t in ts} for ts in JanisShed.get_all_tools()
    }

    Logger.info(f"Preparing documentation for {len(tools)} tools")
    Logger.info(f"Preparing documentation for {len(data_types)} data_types")

    tool_module_index = {}
    dt_module_index = {}
    ROOT_KEY = "root"

    if os.path.exists(tools_dir):
        rmtree(tools_dir)

    for toolname, toolsbyversion in tools.items():
        # tool = tool_vs[0][0]()
        tool_versions = sort_tool_versions(list(toolsbyversion.keys()))
        default_version = tool_versions[0]
        Logger.log(
            f"Preparing {toolname}, found {len(tool_versions)} version[s] ({','.join(tool_versions)})"
        )

        defaulttool = toolsbyversion[default_version]
        if isclass(defaulttool):
            defaulttool = defaulttool()
        try:
            tool_path_components = list(
                filter(
                    lambda a: bool(a),
                    [defaulttool.tool_module(), defaulttool.tool_provider()],
                )
            )
        except Exception as e:
            Logger.critical(f"Failed to generate docs for {toolname}: {e}")
            continue

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
            if isinstance(tool, WorkflowBase):
                tool.get_dot_plot(output_directory=output_dir, log_to_stdout=False)
            if output_str is None:
                Logger.warn(f"Skipping {tool.id()}")
                continue
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
        indexed_submodules_tools = [m.lower() for m in submodule_keys]

        with open(module_filename, "w+") as module_file:
            module_file.write(
                get_tool_toc(
                    alltoolsmap=tools,
                    title=title,
                    intro_text=f"Automatically generated index page for {title}:",
                    subpages=indexed_submodules_tools,
                    tools=module_tools,
                    max_depth=max_depth,
                )
            )

        for submodule in submodule_keys:
            prepare_modules_in_index(
                contents=contents[submodule], title=submodule, dir=f"{dir}/{submodule}/"
            )

    def prepare_dtmodules_in_index(contents, title, dir, max_depth=1):
        module_filename = dir + "/index.rst"
        module_tools = sorted(set(contents[ROOT_KEY] if ROOT_KEY in contents else []))
        submodule_keys = sorted(m for m in contents.keys() if m != ROOT_KEY)
        indexed_submodules_tools = [m.lower() + "/index" for m in submodule_keys]

        with open(module_filename, "w+") as module_file:
            module_file.write(
                get_toc(
                    title=title,
                    intro_text=f"Automatically generated index page for {title}:",
                    subpages=indexed_submodules_tools + module_tools,
                    max_depth=max_depth,
                )
            )

        for submodule in submodule_keys:
            prepare_modules_in_index(
                contents=contents[submodule], title=submodule, dir=f"{dir}/{submodule}/"
            )

    prepare_modules_in_index(tool_module_index, title="Tools", dir=tools_dir)
    prepare_dtmodules_in_index(
        dt_module_index, title="Data Types", dir=dt_dir, max_depth=1
    )


def prepare_templates():
    """
    :return:
    """

    os.makedirs(templates_dir, exist_ok=True)

    templates = get_all_templates()

    for tkey, template in templates.items():
        with open(os.path.join(templates_dir, tkey + ".rst"), "w+") as f:
            f.write(prepare_template(tkey, template))

    introtext = """\

Templates
###########
    
These templates are used to configure Cromwell / CWLTool broadly. For more information, visit `Configuring Janis <https://janis.readthedocs.io/en/latest/references/configuration.html#cromwell>`__.
"""

    with open(os.path.join(templates_dir, "index.rst"), "w+") as f:
        f.write(
            introtext
            + get_toc(
                title="",
                intro_text="List of templates for ``janis-assistant``:",
                subpages=list(templates.keys()),
            )
        )


def generate_pipelines_page():
    import janis_pipelines

    print("Generating pipelines page...")

    modules = [janis_pipelines]
    workflows: List[Workflow] = []

    for module in modules:

        workflows.extend(
            cls()
            for n, cls in list(module.__dict__.items())
            if not n.startswith("__")
            and type(cls) != type
            and isclass(cls)
            and hasattr(cls, "type")
            and cls.type() == ToolType.Workflow
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
        toolstr = prepare_published_pipeline_page(w, [w.version()])
        w.get_dot_plot(output_directory=pipelines_dir, log_to_stdout=False)

        with open(os.path.join(pipelines_dir, w.id().lower() + ".rst"), "w+") as f:
            f.write(toolstr)


if __name__ == "__main__":
    prepare_all_tools()
    prepare_templates()
    generate_pipelines_page()

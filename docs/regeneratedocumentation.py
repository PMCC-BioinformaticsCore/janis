"""
This class regenerates all the tool definitions

It's a bit of a random collection of things that should be refactored:

    - NestedDictionary (store values in nested structure (with root key))
    - RST Helpers
"""
from inspect import isfunction, ismodule

from constants import PROJECT_ROOT_DIR

# Import modules here so that the tool registry knows about them
import bioinformatics
import janis.unix

# Output settings
docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"

###### Shouldn't need to touch below this line #####

import os
from typing import List

import tabulate
from datetime import date

from janis import CommandTool, Logger, Metadata
from janis.tool.tool import Tool
from janis.utils.metadata import Metadata


class NestedDictionaryTypeException(Exception):
    def __init__(self, key, error_type, keys=None):
        keychain = (" for keys: " + ".".join(keys)) if keys else ""
        super(NestedDictionaryTypeException, self).__init__(f"Incorrect type '{error_type}' for key '{key}'{keychain}")
        self.keys = keys
        self.key = key
        self.error_type = error_type


def format_rst_link(text, link):
    return f"`{text} <{link}>`_"


def prepare_tool(tool: Tool):
    # Stuff to list on the documentation page:
    #   - Versions of tools
    #   - Generated command
    #   - Cool if it grouped the tools by vendor
    #   -

    if not tool:
        return None

    # tool_modules = tool.__module__.split(".") # janis._bioinformatics_.tools.$toolproducer.$toolname.$version

    metadata = tool.metadata()
    if not tool.metadata():
        metadata = Metadata()

    fn = tool.friendly_name() if tool.friendly_name() else tool.id()
    en = f" ({tool.id()})" if fn != tool.id() else ""
    tn = fn + en

    formatted_url = format_rst_link(metadata.documentationUrl, metadata.documentationUrl) if metadata.documentationUrl \
        else "*No URL to the documentation was provided*"

    input_headers = ["name", "type", "prefix", "position", "documentation"]

    required_input_tuples = [[i.id(), i.input_type.id(), i.prefix, i.position, i.doc] for i in tool.inputs() if
                             not i.input_type.optional]
    optional_input_tuples = [[i.id(), i.input_type.id(), i.prefix, i.position, i.doc] for i in tool.inputs() if
                             i.input_type.optional]

    formatted_required_inputs = tabulate.tabulate(required_input_tuples, input_headers, tablefmt="rst")
    formatted_optional_inputs = tabulate.tabulate(optional_input_tuples, input_headers, tablefmt="rst")

    output_headers = ["name", "type", "documentation"]
    output_tuples = [[o.id(), o.output_type.id(), o.doc] for o in tool.outputs()]
    formatted_outputs = tabulate.tabulate(output_tuples, output_headers, tablefmt="rst")

    docker_tag = ""
    if isinstance(tool, CommandTool):
        docker_tag = "Docker\n******\n``" + tool.docker() + "``\n"

    return f"""
{fn}
{"=" * len(tn)}
Tool identifier: ``{tool.id()}``

Documentation
-------------

{docker_tag}
URL
******
{formatted_url}

Docstring
*********
{metadata.documentation if metadata.documentation else "*No documentation was provided: " + format_rst_link(
        "contribute one", "https://github.com/illusional") + "*"}

Outputs
-------
{formatted_outputs}

Inputs
------
Find the inputs below

Required inputs
***************

{formatted_required_inputs}

Optional inputs
***************

{formatted_optional_inputs}


*{fn} was last updated on {metadata.dateUpdated if metadata.dateUpdated else "**Unknown**"}*.
*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""


def nested_keys_add(d, keys: List[str], value, root_key):
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
                raise NestedDictionaryTypeException(key=key, error_type=type(d[key]), keys=keys)
            nested_keys_add(d[key], keys[1:], value, root_key=root_key)
            return d
        else:
            d[key] = nested_keys_add({}, keys[1:], value, root_key=root_key)
    except NestedDictionaryTypeException as de:
        raise NestedDictionaryTypeException(key=de.key, error_type=de.error_type, keys=keys)

    return d


def get_toc(title, intro_text, subpages, caption="Contents", max_depth=1):
    prepared_subpages = "\n".join("   " + m for m in subpages)
    return f"""
{title}
{"=" * len(title)}

{intro_text}

.. toctree::
   :maxdepth: {max_depth}
   :caption: {caption}:

{prepared_subpages}

*This page was auto-generated on {date.today().strftime(
        "%d/%m/%Y")}. Please do not directly alter the contents of this page.*
"""

def get_tools():
    import janis.bioinformatics, janis.unix
    modules = [janis.bioinformatics, janis.unix]

    tools = []
    for m in modules:
        try:
            tools_module = m.tools
            q = {n: cls for n, cls in list(tools_module.__dict__.items()) if not n.startswith("__") and type(cls) != type}
            for k in q:
                cls = q[k]
                if isfunction(cls): continue
                if ismodule(cls):
                    print("module: " + str(cls))
                if issubclass(type(cls), CommandTool):
                    print("Found tool")
                    tools.append(cls)
        except Exception as e:
            print(e)
            continue



    # Logger.info(f"Finding modules in {len(files)} files")
    # for file in files:
    #     if os.path.basename(file).startswith("__"):
    #         continue
    #     if os.path.basename(file) in ignore_files:
    #         continue
    #
    #     name = os.path.splitext(file.replace(PROJECT_ROOT_DIR + "/", ""))[0].replace("/", ".")
    #     try:
    #         module = importlib.import_module(name)
    #
    #         for cc in q:
    #             try_register_type(q[cc])
    #
    #     except Exception as e:
    #         Logger.log_ex(e)

    return tools


def prepare_all_tools():
    tools = get_tools()

    Logger.info(f"Preparing documentation for {len(tools)} tools")
    tool_module_index = {}
    ROOT_KEY = "root"

    for tool_vs in tools:
        tool = tool_vs[0][0]()
        Logger.log("Preparing " + tool.id())
        output_str = prepare_tool(tool)

        tool_path_components = list(filter(
            lambda a: bool(a),
            [tool.tool_module(), tool.tool_provider()]
        ))

        path_components = "/".join(tool_path_components)
        output_dir = f"{tools_dir}/{path_components}/"
        output_filename = output_dir + tool.id() + ".rst"

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        nested_keys_add(tool_module_index, tool_path_components, tool.id(), root_key=ROOT_KEY)

        with open(output_filename, "w+") as tool_file:
            tool_file.write(output_str)

        Logger.log("Prepared " + tool.id())

    def prepare_modules_in_index(contents, title="Tools", output_dir=tools_dir):
        module_filename = output_dir + "/index.rst"
        module_tools = sorted(contents[ROOT_KEY] if ROOT_KEY in contents else [])
        submodule_keys = sorted(m for m in contents.keys() if m != ROOT_KEY)
        indexed_submodules_tools = [m + "/index" for m in submodule_keys]

        with open(module_filename, "w+") as module_file:
            module_file.write(get_toc(
                title=title,
                intro_text="Automatically generated index page for {module} tools",
                subpages=indexed_submodules_tools + module_tools,
                max_depth=2
            ))

        for submodule in submodule_keys:
            prepare_modules_in_index(
                contents=contents[submodule],
                title=submodule,
                output_dir=f"{output_dir}/{submodule}/"
            )

    prepare_modules_in_index(tool_module_index)


prepare_all_tools()

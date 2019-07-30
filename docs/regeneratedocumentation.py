"""
This class regenerates all the tool definitions

It's a bit of a random collection of things that should be refactored:

    - NestedDictionary (store values in nested structure (with root key))
    - RST Helpers
"""
from requests.utils import requote_uri
from distutils.version import StrictVersion
from inspect import isfunction, ismodule, isabstract, isclass

import janis
from constants import PROJECT_ROOT_DIR
from janis_core import Array, DataType, Workflow, CommandTool, Tool, Metadata, Logger
from janis_core.types.common_data_types import all_types

# Import modules here so that the tool registry knows about them

# Output settings
docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"
dt_dir = docs_dir + "datatypes/"

modules = [janis.bioinformatics, janis.unix]

###### Shouldn't need to touch below this line #####

import os
from typing import List, Set, Type, Tuple, Dict

import tabulate
from datetime import date


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

    metadata = tool.metadata()
    if not tool.metadata():
        metadata = Metadata()

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


Metadata
********

Author: {metadata.creator if metadata.creator else "**Unknown**"}


*{fn} was last updated on {metadata.dateUpdated if metadata.dateUpdated else "**Unknown**"}*.
*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""


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


def get_tools_and_datatypes():

    tools: Dict[str, Dict[str, Type[Tool]]] = {}
    data_types: Set[Type[DataType]] = set(all_types)

    for m in modules:
        # noinspection PyTypeChecker
        tt, dt = get_tool_from_module(m.tools)
        # noinspection PyTypeChecker
        tdt, ddt = get_tool_from_module(m.data_types)

        tools = merge_versioned_tools_dicts(tools, merge_versioned_tools_dicts(tt, tdt))
        data_types = data_types.union(dt).union(ddt)

    return tools, list(data_types)


def merge_versioned_tools_dicts(d1: dict, d2: Dict[str, Dict[str, Type[Tool]]]):
    d = {**d1}
    for k, v in d2.items():
        if k in d:
            d[k] = {**d[k], **v}
        else:
            d[k] = v
    return d


def get_tool_from_module(
    module, seen_modules=None
) -> Tuple[Dict[str, Dict[str, Type[Tool]]], Set[Type[DataType]]]:
    q = {
        n: cls
        for n, cls in list(module.__dict__.items())
        if not n.startswith("__") and type(cls) != type
    }

    # name: version: Type[Tool]
    tools: Dict[str, Dict[str, Type[Tool]]] = {}
    data_types: Set[Type[DataType]] = set()

    if seen_modules is None:
        seen_modules = set()

    for k in q:
        cls = q[k]
        try:
            if hasattr(cls, "__name__"):
                if cls.__name__ in seen_modules:
                    continue
                seen_modules.add(cls.__name__)

            if isfunction(cls):
                continue
            if ismodule(cls):
                t, d = get_tool_from_module(cls, seen_modules)
                tools = merge_versioned_tools_dicts(tools, t)
                data_types = data_types.union(d)
            elif isabstract(cls):
                continue
            elif not isclass(cls):
                continue
            elif issubclass(cls, CommandTool):
                print("Found commandtool: " + cls.tool())
                v = cls.version() if cls.version() else "None"
                nested_keys_add(tools, [cls.tool(), v], cls)
            elif issubclass(cls, Workflow):
                c = cls()
                print("Found workflow: " + c.id())
                v = c.version() if c.version() else "None"
                nested_keys_add(tools, [cls().id(), v], cls)

            elif issubclass(cls, DataType):
                print("Found datatype: " + cls().id())
                data_types.add(cls)
        except Exception as e:
            print(f"{str(e)} for type {type(cls)}")
            # print(traceback.format_exc())
            continue

    return tools, data_types


def sort_tool_versions(versions: List[str]) -> List[str]:
    try:
        return sorted(versions, key=StrictVersion, reverse=True)
    except:
        return sorted(versions, reverse=True)


def prepare_all_tools():
    tools, data_types = get_tools_and_datatypes()

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
            output_str = prepare_tool(tool(), tool_versions, not isprimary)
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

        dt = d()
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


prepare_all_tools()

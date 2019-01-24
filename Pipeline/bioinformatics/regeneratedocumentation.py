from datetime import date

from Pipeline import CommandTool, Logger
import tabulate
# from Pipeline.bioinformatics.tools.bwa.mem.latest import BwaMemLatest
# from Pipeline.bioinformatics.tools.gatk4.haplotypecaller.latest import Gatk4HaplotypeCallerLatest

from Pipeline.tool.registry import get_tools
from Pipeline.tool.tool import Tool
from Pipeline.utils.metadata import Metadata
from constants import PROJECT_ROOT_DIR

docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"


def format_rst_link(text, link):
    return f"`{text} <{link}>`_"


def prepare_tool(tool: Tool):

    # Stuff to list on the documentation page:
    #   - Versions of tools
    #   - Generated command
    #   - Cool if it grouped the tools by vendor
    #   -

    if not tool:
        return None, None, None

    tool_module = tool.__module__.split(".")[1]     # Pipeline._bioinformatics_.tools.$toolproducer.$toolname.$version

    tool_dir = tools_dir + tool_module + "/" + tool.id() + ".rst"

    metadata = tool.metadata()
    if not metadata:
        metadata = Metadata()

    fn = tool.friendly_name() if tool.friendly_name() else tool.id()
    en = f" ({tool.id()})" if fn != tool.id() else ""
    tn = fn + en

    formatted_url = format_rst_link(metadata.documentationUrl, metadata.documentationUrl) if metadata.documentationUrl \
        else "*No URL to the documentation was provided*: " + format_rst_link("contribute one", "https://github.com/illusional")

    input_headers = ["name", "type", "prefix", "position", "documentation"]
    input_tuples = [[i.id(), i.input_type.id(), i.prefix, i.position, i.doc] for i in tool.inputs()]
    formatted_inputs = tabulate.tabulate(input_tuples, input_headers, tablefmt="rst")

    output_headers = ["name", "type", "documentation"]
    output_tuples = [[o.id(), o.output_type.id(), o.doc] for o in tool.outputs()]
    formatted_outputs = tabulate.tabulate(output_tuples, output_headers, tablefmt="rst")

    docker_tag = ""
    if isinstance(tool, CommandTool):
        docker_tag = "Docker\n******\n``" + tool.docker() + "``\n"

    return tool_module, tool_dir, f"""
{fn}
{"=" * len(tn)}
*{tool_module}*{en}

Documentation
-------------

{docker_tag}
URL
******
{formatted_url}

Docstring
*********
{metadata.documentation if metadata.documentation else tool.doc()}

Outputs
-------
{formatted_outputs}

Inputs
------
{formatted_inputs}

*{fn} was last updated on {metadata.dateUpdated if metadata.dateUpdated else "**Unknown**"}*.

*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""


def prepare_all_tools():
    import Pipeline.bioinformatics
    # tools = [[Gatk4HaplotypeCallerLatest], [BwaMemLatest]]
    tools = get_tools()

    Logger.info(f"Preparing documentation for {len(tools)} tools")
    tool_module_index = {}

    for tool_vs in tools:
        tool = tool_vs[0][0]()
        Logger.log("Preparing " + tool.id())
        module, output_filename, output_str, = prepare_tool(tool)

        if module in tool_module_index: tool_module_index[module].append(tool.id())
        else: tool_module_index[module] = [tool.id()]

        with open(output_filename, "w+") as tool_file:
            tool_file.write(output_str)
        Logger.log("Prepared " + tool.id())

    for module in tool_module_index:
        module_filename = tools_dir + module + "/index.rst"
        module_list = "\n".join("   " + m for m in sorted(tool_module_index[module]))
        with open(module_filename, "w+") as module_file:
            module_file.write(f"""
{module}
{"=" * len(module)}

Automatically generated index page for {module} tools.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

{module_list}
""")

    modules = "\n".join("   " + k + "/index" for k in sorted(tool_module_index.keys()))
    tool_index_page = f"""
Tools
======

.. toctree::
   :maxdepth: 2
   :caption: Contents:

{modules}

*This page was auto-generated. Please do not directly alter the contents of this page.*
"""
    with open(tools_dir + "index.rst", "w+") as tool_index_file:
        tool_index_file.write(tool_index_page)


prepare_all_tools()

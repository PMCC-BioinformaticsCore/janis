from Pipeline import CommandTool
from Pipeline.bioinformatics.tools.bwa.mem.latest import BwaMemLatest
from Pipeline.bioinformatics.tools.gatk4.haplotypecaller.latest import Gatk4HaplotypeCallerLatest

from Pipeline.tool.registry import get_tools
from constants import PROJECT_ROOT_DIR

docs_dir = PROJECT_ROOT_DIR + "/docs/"
tools_dir = docs_dir + "tools/"

def format_rst_link(text, link):
    return f"`{text} <{link}/>`_"

def prepare_tool(tool: CommandTool):

    tool_module = tool.__module__.split(".")[1]     # Pipeline._bioinformatics_.tools.$toolproducer.$toolname.$version

    tool_dir = tools_dir + tool_module + "/" + tool.id() + ".rst"

    formatted_url = format_rst_link(tool.docurl(), tool.docurl()) if tool.docurl() \
        else "*No URL to the documentation was provided, " + format_rst_link("contribute one", "github.com/illusional")

    return tool_module, tool_dir, f"""
{tool.id()}
{"=" * len(tool.id())}
*{tool_module}*

Documentation
-------------

URL
******
{formatted_url}

Docstring
*********
{tool.doc()}

    """


def prepare_all_tools():
    import Pipeline.bioinformatics

    tools = [[Gatk4HaplotypeCallerLatest], [BwaMemLatest]]
    tool_module_index = {}

    for tool_vs in tools:
        tool = tool_vs[0]()
        module, output_filename, output_str, = prepare_tool(tool)

        if module in tool_module_index: tool_module_index[module].append(tool.id())
        else: tool_module_index[module] = [tool.id()]

        with open(output_filename, "w+") as tool_file:
            tool_file.write(output_str)

    for module in tool_module_index:
        module_filename = tools_dir + module + "/index.rst"
        module_list = "\n".join("   " + m for m in tool_module_index[module])
        with open(module_filename, "w+") as module_file:
            module_file.write(f"""
{module}
{"=" * len(module)}

Automatically generated index page for {module} related tools.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

{module_list}
""")

    modules = "\n".join("   " + k + "/index" for k in tool_module_index.keys())
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




if __name__ == "__main___":
    prepare_all_tools()

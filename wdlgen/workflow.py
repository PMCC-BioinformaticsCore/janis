from typing import List, Any
import re

from wdlgen.common import Input, Output
from wdlgen.util import WdlBase
from wdlgen.workflowcall import WorkflowCallBase


class Workflow(WdlBase):

    def __init__(self, name, inputs: List[Input]=None, outputs: List[str]=None, calls: List[WorkflowCallBase]=None,
                 imports: List[Any]=None):
        """

        :param name:
        :param inputs:
        :param outputs:
        :param calls: List[WorkflowXCallBase]
        :param imports:
        :type imports: List[WorkflowImport]
        """

        # validate

        self.name = name.replace("-", "_")
        self.inputs = inputs if inputs else []
        self.outputs = outputs if outputs else []
        self.calls = calls if calls else []
        self.imports = imports if imports else []

        self.format = """
{imports_block}

workflow {name} {{
{inputs_block}
{call_block}
{output_block}
}}""".strip()

    def get_string(self):
        tb = "  "

        name = self.name
        inputs_block, output_block, call_block, imports_block = "", "", "", ""

        if self.inputs:
            ins = []
            for i in self.inputs:
                wd = i.get_string()
                if isinstance(wd, list):
                    ins.extend(tb + ii for ii in wd)
                else:
                    ins.append(tb + wd)
            inputs_block = "\n".join(ins)

        if self.outputs:
            outs = []
            # either str | Output | list[str | Output]
            for o in self.outputs:
                if isinstance(o, Output):
                    wd = o.get_string()
                    if isinstance(wd, list):
                        outs.extend((2 * tb) + w for w in wd)
                    else:
                        outs.append((2 * tb) + wd)
                else:
                    outs.append(str(o))
            output_block = "{tb}output {{\n{outs}\n{tb}}}".format(
                tb=tb,
                outs="\n".join(outs)
            )

        if self.calls:
            call_block = "\n".join(c.get_string(indent=1) for c in self.calls)
            # map calls

            # do the imports as well

        if self.imports:
            imports_block = "\n".join(i.get_string() for i in self.imports)

        return self.format.format(
            name=name,
            inputs_block=inputs_block,
            output_block=output_block,
            call_block=call_block,
            imports_block=imports_block
        )

    class WorkflowImport(WdlBase):
        def __init__(self, name: str, alias: str, tools_dir="tools/"):
            self.name = name
            self.alias = alias
            self.tools_dir = tools_dir

            if tools_dir and not self.tools_dir.endswith("/"):
                tools_dir += "/"

        def get_string(self):
            as_alias = " as " + self.alias if self.alias else ""
            return "import \"{tools_dir}{tool}.wdl\"{as_alias}".format(
                tools_dir=self.tools_dir if self.tools_dir else "",
                tool=self.name,
                as_alias=as_alias
            )

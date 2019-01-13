from typing import List

from wdlgen.common import Input, Output
from wdlgen.util import WdlBase
from wdlgen.workflowcall import WorkflowCallBase


class Workflow(WdlBase):

    def __init__(self, name, inputs: List[Input]=None, outputs: List[str]=None, calls: List[WorkflowCallBase]=None):

        self.name = name
        self.inputs = inputs if inputs else []
        self.outputs = outputs if outputs else []
        self.calls = calls if calls else []

        self.format = """
{imports_block}

workflow {name} {{
{inputs_block}
{call_block}
{output_block}
}}""".strip()

    def wdl(self, tools_dir="tools/"):
        tb = "  "

        name = self.name
        inputs_block, output_block, call_block, imports_block = "", "", "", ""

        if self.inputs:
            ins = []
            for i in self.inputs:
                wd = i.wdl()
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
                    wd = o.wdl()
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
            call_block = "\n".join(c.wdl(indent=1) for c in self.calls)
            # map calls

            # do the imports as well

        return self.format.format(
            name=name,
            inputs_block=inputs_block,
            output_block=output_block,
            call_block=call_block,
            imports_block=imports_block
        )


"""
        get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])

        tools: List[Tool] = [s.step.tool() for s in self._steps]
        tool_name_to_tool: Dict[str, Tool] = {t.id().lower(): t for t in tools}
        tool_name_to_alias = {}
        steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in self._steps}

        aliases = set()

        for tool in tool_name_to_tool:
            a = get_alias(tool).upper()
            s = a
            idx = 2
            while s in aliases:
                s = a + str(idx)
                idx += 1
            aliases.add(s)
            tool_name_to_alias[tool] = s

        tab_char = '  '
        nline_char = '\n'

        import_str = "import \"tools/{tool_file}.wdl\" as {alias}"
        input_str = "{tb}{data_type} {identifier}"
        step_str = "{tb}call {tool_file}.{tool} as {alias} {{ input: {tool_mapping} }}"
        output_str = "{tb2}{data_type} {identifier} = {alias}.{outp}"

        imports = [import_str.format(
            tb=tab_char,
            tool_file=t,
            alias=tool_name_to_alias[t.lower()].upper()
        ) for t in tool_name_to_tool]

        inputs = [input_str.format(
            tb=tab_char,
            data_type=i.input.data_type.wdl(),
            identifier=i.id()
        ) for i in self._inputs]

        steps = [step_str.format(
            tb=tab_char,
            tool_file=tool_name_to_alias[s.step.tool().id().lower()].upper(),
            tool=s.step.tool().wdl_name(),
            alias=s.id(),
            tool_mapping=', '.join(s.wdl_map())  # [2 * tab_char + w for w in s.wdl_map()])
        ) for s in self._steps]

        outputs = [output_str.format(
            tb2=2 * tab_char,
            data_type=o.output.data_type.wdl(),
            identifier=o.id(),
            alias=first_value(first_value(o.connection_map).source_map).start.id(),
            outp=first_value(first_value(o.connection_map).source_map).stag
        ) for o in self._outputs]

        # imports = '\n'.join([f"import \"tools/{t}.wdl\" as {tool_name_to_alias[t.lower()].upper()}" for t in tool_name_to_tool])
        # inputs = '\n'.join([f"{tab_char}{i.input.data_type.wdl()} {i.id()}" for i in self._inputs])
        # steps = '\n'.join([f"{tab_char}call {tool_name_to_alias[s.step.get_tool().tool().lower()].upper()}"
        #                    f".{s.step.get_tool().wdl_name()} as {s.id()} {{input: \n{(',' + nline_char).join([2 * tab_char + w for w in s.wdl_map()])}\n{tab_char}}}" for s in self._steps])
        # outputs = '\n'.join([f"{2*tab_char}{o.output.data_type.wdl()} {o.id()} = {steps_to_alias[next(iter(o.connection_map.values()))[0].split('/')[0].lower()].lower()}.{o.id()}" for o in self._outputs])

        workflow = f\"""
{nline_char.join(imports)}

workflow {self.identifier} {{
{nline_char.join(inputs)}

{nline_char.join(steps)}

{tab_char}output {{
{nline_char.join(outputs)}
{tab_char}}}
}}"\"
        tools = {t.id(): t.wdl() for t in tools}

        inp = {f"{self.identifier}.{i.id()}": i.input.wdl_input() for i in self._inputs}

        return workflow, inp, tools
"""

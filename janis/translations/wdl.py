import os
from typing import List, Dict, Optional, Any

import wdlgen as wdl
from janis.workflow.step import StepNode
from janis.graph.stepinput import Edge, StepInput
from janis.types import InputSelector, WildcardSelector, CpuSelector, MemorySelector
from janis.types.common_data_types import Stdout, Array, Boolean, Filename, File
from janis.utils.validators import Validators
from janis.tool.commandtool import CommandTool
from janis.utils import first_value, convert_expression_to_wdl
from janis.tool.tool import Tool, ToolInput

from janis.utils.logger import Logger


def value_or_default(ar, default):
    return ar if ar else default


def dump_wdl(workflow, to_console=True, to_disk=False, with_docker=False, with_hints=False,
             with_resource_overrides=False, write_inputs_file=False, should_validate=False, should_zip=True):
    import json
    wf_wdl, inp_dict, tool_dicts = translate_workflow(workflow, with_docker=with_docker,
                                                      with_resource_overrides=with_resource_overrides)

    wf_str = wf_wdl.get_string()
    inp_str = json.dumps(inp_dict, sort_keys=True, indent=4, separators=(',', ': '))
    tls_strs = [("tools/" + t + ".wdl", tool_dicts[t].get_string()) for t in tool_dicts]

    if to_console:
        print("=== WORKFLOW ===")
        print(wf_str)
        print("\n=== INPUTS ===")
        print(inp_str)

        print("\n=== TOOLS ===")
        [print(t[1]) for t in tls_strs]

    d = os.path.expanduser("~") + f"/Desktop/{workflow.id()}/wdl/"

    if write_inputs_file:
        with open(d + workflow.id() + "-job.json", "w+") as wdlfile:
            Logger.log(f"Writing {workflow.id()}-job.yml to disk")
            wdlfile.write(inp_str)
            # ruamel.yaml.dump(inp_dict, cwl, default_flow_style=False)
            Logger.log(f"Written {workflow.id()}-job.yml to disk")
    else:
        Logger.log("Skipping writing input (yaml) job file")

    if to_disk:
        d_tools = d + "tools/"

        if not os.path.isdir(d):
            os.makedirs(d)
        if not os.path.isdir(d_tools):
            os.makedirs(d_tools)

        os.chdir(d)
        wf_filename = d + workflow.id() + ".wdl"
        with open(wf_filename, "w+") as wdlfile:
            Logger.log(f"Writing {workflow.id()}.wdl to disk")
            wdlfile.write(wf_str)
            # ruamel.yaml.dump(wf_dict, cwl, default_flow_style=False)
            Logger.log(f"Written {workflow.id()}.wdl to disk")

        # z = zipfile.ZipFile(d + "tools.zip", "w")
        for (tool_filename, tool) in tls_strs:
            with open(d + tool_filename, "w+") as wdlfile:
                Logger.log(f"Writing {tool_filename} to disk")
                wdlfile.write(tool)
                Logger.log(f"Written {tool_filename} to disk")

        import subprocess

        Logger.info("Validing outputted WDL")

        import subprocess
        if should_validate:
            Logger.info("Validing outputted CWL")

            womtool_result = subprocess.run(["java", "-jar", "$womtool", "--validate", wf_filename])
            if womtool_result.returncode == 0:
                Logger.info("Exported workflow is valid WDL.")
            else:
                Logger.critical(womtool_result.stderr)

        if should_zip:
            Logger.info("Zipping tools")
            os.chdir(d)

            zip_result = subprocess.run(["zip", "-r", "tools.zip", "tools/"])
            if zip_result.returncode == 0:
                Logger.info("Zipped tools")
            else:
                Logger.critical(zip_result.stderr)

    return wf_str, inp_str, tls_strs


def build_aliases(steps):
    """
    :param steps:
    :return:
    """

    get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])
    aliases = set()

    tools: List[Tool] = [s.step.tool() for s in steps]
    tool_name_to_tool: Dict[str, Tool] = {t.id().lower(): t for t in tools}
    tool_name_to_alias = {}
    steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in steps}

    for tool in tool_name_to_tool:
        a = get_alias(tool).upper()
        s = a
        idx = 2
        while s in aliases:
            s = a + str(idx)
            idx += 1
        aliases.add(s)
        tool_name_to_alias[tool] = s

    return tool_name_to_alias, steps_to_alias


def get_secondary_tag_from_original_tag(original, secondary):
    secondary_without_punctuation = secondary.replace(".", "").replace("^", "")
    return original + "_" + secondary_without_punctuation


def apply_secondary_file_format_to_filename(filename: Optional[str], secondary_file: str):
    if not filename: return None

    fixed_sec = secondary_file.lstrip("^")
    leading = len(secondary_file) - len(fixed_sec)
    if leading <= 0:
        return filename + fixed_sec

    split = filename.split(".")
    return ".".join(split[:-min(leading, len(split) - 1)]) + fixed_sec


def translate_workflow(wf, with_docker=True, is_nested_tool=False, with_resource_overrides=False):
    """

    :param with_docker:
    :param is_nested_tool:
    :return:
    """
    from janis.workflow.workflow import Workflow

    # Notes:
    #       All wdlgen classes have a .get_string(**kwargs) function
    #       The wdlgen Workflow class requires a

    w = wdl.Workflow(wf.identifier, version="development")
    tools: List[Tool] = [s.step.tool() for s in wf._steps]

    wtools = {}  # Store all the tools by their name in this dictionary
    tool_aliases, step_aliases = build_aliases(wf._steps)  # Generate call and import aliases

    # Convert self._inputs -> wdl.Input
    for i in wf._inputs:
        wd = i.input.data_type.wdl()
        w.inputs.append(wdl.Input(wd, i.id(), i.input.data_type.default()))
        is_array = isinstance(i.input.data_type, Array)
        if i.input.data_type.secondary_files() or \
                (is_array and i.input.data_type.subtype().secondary_files()):
            secs = i.input.data_type.secondary_files() if not is_array else i.input.data_type.subtype().secondary_files()
            w.inputs.extend(
                wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s)) for s in secs)

    resource_inputs = []
    if with_resource_overrides:
        resource_inputs = build_wdl_resource_inputs(wf)
        w.inputs.extend(resource_inputs)

    # Convert self._outputs -> wdl.Output
    for o in wf._outputs:
        w.outputs.append(wdl.Output(
            o.output.data_type.wdl(),
            o.id(),
            "{a}.{b}".format(  # Generate fully qualified stepid.tag identifier (MUST be step node)
                a=first_value(first_value(o.connection_map).source_map).start.id(),
                b=first_value(first_value(o.connection_map).source_map).stag
            )))
        if o.output.data_type.secondary_files():
            w.outputs.extend(wdl.Output(
                o.output.data_type.wdl(),
                get_secondary_tag_from_original_tag(o.id(), s),
                "{a}.{b}".format(  # Generate fully qualified stepid.tag identifier (MUST be step node)
                    a=first_value(first_value(o.connection_map).source_map).start.id(),
                    b=get_secondary_tag_from_original_tag(first_value(first_value(o.connection_map).source_map).stag, s)
                )
            ) for s in o.output.data_type.secondary_files())

    # Generate import statements (relative tool dir is later?)
    w.imports = [
        wdl.Workflow.WorkflowImport(
            t.id(),
            tool_aliases[t.id().lower()].upper(),
            None if is_nested_tool else "tools/")
        for t in tools]

    # Step[] -> (wdl.Task | wdl.Workflow)[]
    for s in wf._steps:
        t = s.step.tool()

        if isinstance(t, Workflow):
            wf_wdl, _, wf_tools = translate_workflow(t, with_docker=with_docker, is_nested_tool=True,
                                                     with_resource_overrides=with_resource_overrides)
            wtools[t.id()] = wf_wdl
            wtools.update(wf_tools)

        elif isinstance(t, CommandTool):
            wtools[t.id()] = translate_tool(t, with_docker=with_docker, with_resource_overrides=with_resource_overrides)

        resource_overrides = {}
        for r in resource_inputs:
            if not r.name.startswith(s.id()): continue

            resource_overrides[r.name[(len(s.id()) + 1):]] = r.name
        call = translate_step_node(s, tool_aliases[t.id().lower()].upper() + "." + t.id(), s.id(), resource_overrides)

        w.calls.append(call)

    inp = {}
    for i in wf._inputs:
        inp_key = f"{wf.id()}.{i.id()}"
        inp_val = i.input.wdl_input()
        inp[inp_key] = inp_val
        if i.input.data_type.secondary_files():
            for sec in i.input.data_type.secondary_files():
                inp[get_secondary_tag_from_original_tag(inp_key, sec)] = \
                    apply_secondary_file_format_to_filename(inp_val, sec)

    return w, inp, wtools




def get_input_value_from_potential_selector(value, tool_id, string_environment=True):
    if not value:
        return None
    if isinstance(value, str):
        return value if string_environment else f'"{value}"'
    elif isinstance(value, int) or isinstance(value, float):
        return value
    elif isinstance(value, InputSelector):
        return translate_input_selector(selector=value, string_environment=string_environment)
    elif isinstance(value, WildcardSelector):
        raise Exception(f"A wildcard selector cannot be used as an argument value for '{tool_id}'")
    elif isinstance(value, CpuSelector):
        return translate_cpu_selector(value, string_environment=string_environment)
    elif isinstance(value, MemorySelector):
        return translate_mem_selector(value, string_environment=string_environment)
    elif callable(getattr(value, "wdl", None)):
        return value.wdl()

    raise Exception("Could not detect type %s to convert to input value" % type(value))


def translate_input_selector(selector: InputSelector, string_environment=True):
    if not selector.input_to_select: raise Exception("No input was selected for input selector: " + str(selector))

    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""

    if string_environment:
        return f"{pre}${{{selector.input_to_select}}}{suf}"
    else:
        pref = ('"%s" + ' % pre) if pre else ""
        suff = (' + "%s"' % suf) if suf else ""
        return pref + selector.input_to_select + suff


def translate_cpu_selector(selector: CpuSelector, string_environment=True):
    if string_environment:
        return "${runtime_cpu}"
    return "runtime_cpu"


def translate_mem_selector(selector: MemorySelector, string_environment=True):
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""

    val = "floor(runtime_memory)"

    if string_environment:
        return f"{pre}${{{val}}}{suf}"
    else:
        pref = ('"%s" + ' % pre) if pre else ""
        suff = (' + "%s"' % suf) if suf else ""
        return pref + val + suff


def translate_wildcard_selector(selector: WildcardSelector):
    if not selector.wildcard: raise Exception("No wildcard was selected for wildcard selector: " + str(selector))
    return f"glob(\"{selector.wildcard}\")"

def translate_tool_str(tool, with_docker):
    return translate_tool(tool, with_docker=with_docker).get_string()

def translate_tool(tool, with_docker, with_resource_overrides=False):
    if not Validators.validate_identifier(tool.id()):
        raise Exception(f"The identifier '{tool.id()}' for class '{tool.__class__.__name__}' was not validated by "
                        f"'{Validators.identifier_regex}' (must start with letters, and then only contain letters, "
                        f"numbers or an underscore)")

    ins, outs = [], []
    for i in tool.inputs():
        wd = i.input_type.wdl()
        if isinstance(wd, list):
            ins.extend(wdl.Input(w, i.id()) for w in wd)
        else:
            default = i.default if i.default else i.input_type.default()
            default = get_input_value_from_potential_selector(default, tool.id(), string_environment=False)

            ins.append(wdl.Input(wd, i.id(), default, requires_quotes=False))

            sec = value_or_default(i.input_type.subtype().secondary_files() if isinstance(i.input_type, Array)
                                   else i.input_type.secondary_files(), default=[])
            ins.extend(wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s))
                       for s in sec)

    for o in tool.outputs():
        outs.extend(translate_output_node(o, tool))

    command_ins = []
    if tool.inputs():
        for i in tool.inputs():
            cmd = translate_command_input(i)
            if cmd:
                command_ins.append(cmd)

    command_args = None
    if tool.arguments():
        command_args = []
        for a in tool.arguments():
            if a.value is None:
                val = None
            val = get_input_value_from_potential_selector(a.value, tool.id())
            command_args.append(wdl.Task.Command.CommandArgument(a.prefix, val, a.position))

    commands = [prepare_move_statement_for_input_to_localise(ti) for ti in tool.inputs() if ti.localise_file]

    bc = " ".join(tool.base_command()) if isinstance(tool.base_command(), list) else tool.base_command()

    commands.append(wdl.Task.Command(bc, command_ins, command_args))

    r = wdl.Task.Runtime()
    if with_docker:
        r.add_docker(tool.docker())

    if with_resource_overrides:
        # generate resource inputs, for memory, cpu and disk at the moment
        ins.extend([
            wdl.Input(wdl.WdlType.parse_type("Int?"), "runtime_cpu"),
            wdl.Input(wdl.WdlType.parse_type("String?"), "runtime_memory"),
            wdl.Input(wdl.WdlType.parse_type("String?"), "runtime_disks"),
        ])

        r.add_cpus("runtime_cpu")
        r.add_memory("${runtime_memory}")
        r.kwargs["disks"] = "runtime_disks"

    return wdl.Task(tool.id(), ins, outs, commands, r, version="development")


def prepare_move_statement_for_input_to_localise(ti: ToolInput):
    it = ti.input_type

    if issubclass(type(it), File):
        return wdl.Task.Command(f"mv ${{{ti.id()}}} -t .")
    if isinstance(it, Array) and issubclass(type(it.subtype()), File):
        return wdl.Task.Command("mv ${{sep=' ' {s}}} -t .".format(s=ti.id()))

    raise Exception(f"WDL is unable to localise type '{type(it)}'")


def translate_command_input(tool_input: ToolInput):

    # make sure it has some essence of a command line binding, else we'll skip it
    # TODO: make a property on ToolInput (.bind_to_commandline) and set default to true
    if not (tool_input.position is not None or tool_input.prefix): return None

    name = tool_input.id()
    optional = tool_input.input_type.optional
    position = tool_input.position
    separate_value_from_prefix = tool_input.separate_value_from_prefix
    prefix = tool_input.prefix
    true = None
    sep = tool_input.separator

    if tool_input.localise_file:
        name = "basename(%s)" % name

    separate_arrays = tool_input.nest_input_binding_on_array
    if separate_arrays is None and sep is None and isinstance(tool_input.input_type, Array):
        separate_arrays = True

    if isinstance(tool_input.input_type, Boolean):
        true = tool_input.prefix
        prefix = None

    return wdl.Task.Command.CommandInput(
        name=name,
        optional=optional,
        prefix=prefix,
        position=position,
        separate_value_from_prefix=separate_value_from_prefix if separate_value_from_prefix is not None else True,
        # default=default,
        true=true,
        separator=sep,
        separate_arrays=separate_arrays,
    )


def translate_output_node(o, tool) -> List[wdl.Output]:
    if isinstance(o.output_type, Stdout):
        base_expression = "stdout()"
        return [wdl.Output(o.output_type.wdl(), o.id(), base_expression)]

    elif isinstance(o.glob, InputSelector):
        base_expression = translate_input_selector(o.glob, string_environment=False)

        tool_in = tool.inputs_map().get(o.glob.input_to_select)
        if not tool_in:
            raise Exception(f"The InputSelector for tool '{tool.id()}.{o.id()}' did not select an input (tried: '{o.glob.input_to_select}')")

        expression = base_expression if not tool_in.localise_file else f'basename({base_expression})'

        outputs = [wdl.Output(o.output_type.wdl(), o.id(), expression)]
        for s in value_or_default(o.output_type.secondary_files(), []):
            sec_expression = None
            if "^" not in s:
                # do stuff here
                sec_expression = f'{expression} + "{s.replace("^", "")}"'

            elif isinstance(tool_in.input_type, Filename) and tool_in.input_type.extension:
                # use the wdl function: sub
                sec_expression = 'sub({inp}, "\\\\{old_ext}$", "{new_ext}")'\
                    .format(inp=expression, old_ext=tool_in.input_type.extension, new_ext=s.replace("^", ""))

            elif File().can_receive_from(tool_in.input_type):
                # use basename
                sec_expression = f'basename({expression}, \"{tool_in.input_type.extension}\") + "{s.replace("^", "")}"'

            outputs.append(wdl.Output(
                o.output_type.wdl(),
                get_secondary_tag_from_original_tag(o.id(), s),
                sec_expression
            ))
        return outputs

    elif isinstance(o.glob, WildcardSelector):
        base_expression = translate_wildcard_selector(o.glob)
        if not isinstance(o.output_type, Array):
            Logger.warn(f"The command tool '{tool.id()}.{o.tag}' used a star-bind (*) glob to find the output, "
                        f"but the return type was not an array. For WDL, the first element will be used, "
                        f"ie: '{base_expression}[0]'")
            base_expression += "[0]"
        wdl_type = wdl.WdlType.parse_type(o.output_type.wdl())
        outputs = [wdl.Output(wdl_type, o.id(), base_expression)]

        secondary = o.output_type.secondary_files()
        if secondary:
            outputs.extend(wdl.Output(wdl_type, get_secondary_tag_from_original_tag(o.id(), s), base_expression)
                           for s in o.output_type.secondary_files())
        return outputs

    else:
        Logger.warn(f"Tool '{tool.id()}' has the non-selector glob: '{o.glob}', this is deprecated. "
                    f"Please use the WildcardSelector to build output for '{o.id()}'")
        glob = convert_expression_to_wdl(o.glob)
        if glob is not None and "*" in glob:
            glob = f'glob({glob})'
            if not isinstance(o.output_type, Array):
                Logger.warn(f"The command tool '{tool.id()}.{o.tag}' used a star-bind (*) glob to find the output, "
                            f"but the return type was not an array. For WDL, the first element will be used, "
                            f"ie: '{glob}[0]'")
                glob = glob + "[0]"
        base_expression = glob
        wdl_type = wdl.WdlType.parse_type(o.output_type.wdl())
        outputs = [wdl.Output(wdl_type, o.id(), base_expression)]

        secondary = o.output_type.secondary_files()
        if secondary:
            outputs.extend(wdl.Output(wdl_type, get_secondary_tag_from_original_tag(o.id(), s), base_expression)
                           for s in o.output_type.secondary_files())
        return outputs


def translate_step_node(node, step_identifier: str, step_alias: str, resource_overrides: Dict[str, str])\
        -> wdl.WorkflowCallBase:

    ins = node.inputs()

    # One step => One WorkflowCall. We need to traverse the edge list to see if there's a scatter
    # then we can build up the WorkflowCall / ScatterCall
    scatterable: List[StepInput] = []

    for k in node.inputs():
        if not (k in node.connection_map and node.connection_map[k].has_scatter()): continue
        step_input: StepInput = node.connection_map[k]
        if step_input.multiple_inputs or isinstance(step_input.source(), list):
            raise NotImplementedError(f"The edge '{step_input.dotted_source()}' on node '{node.id()}' has multiple "
                                      f"inputs, and I don't know how this can be implemented in WDL")
        scatterable.append(step_input)

    # We need to replace the scatterable key(s) with some random variable, eg: for i in iterable:
    ordered_variable_identifiers = ["i", "j", "k", "x", "y", "z", "a", "b", "c", "ii", "jj", "kk", "xx", "yy", "zz"]
    old_to_new_identifier = {k.dotted_source(): k.dotted_source() for k in scatterable
                             if not isinstance(k.dotted_source(), list)}
    current_identifiers = set(old_to_new_identifier.values())

    # We'll wrap everything in the scatter block later, but let's replace the fields we need to scatter
    # with the new scatter variable (we'll try to guess one based on the fieldname. We might need to eventually
    # pass the workflow inputs to make sure now conflict will arise.
    # Todo: Pass Workflow input tags to wdl scatter generation to ensure scatter var doesn't conflict with inputs

    for s in scatterable:
        current_identifiers.remove(s.dotted_source())
        e: Edge = first_value(s.source_map)
        new_var = e.stag[0] if e.stag else e.source()[0]

        while new_var in current_identifiers:
            new_var = ordered_variable_identifiers.pop(0)
        old_to_new_identifier[s.dotted_source()] = new_var
        current_identifiers.add(new_var)

    # Let's map the inputs, to the source.
    # We're using a dictionary for the map atm, but WDL requires the format:
    #       fieldName: sourceCall.Output

    # So, the current way of mixing accessory files is not really supported, but a little complicated
    # basically, if our scatterable edge contains secondary files, they'll all be arrays of separate files, eg:
    #
    # File[] bams = [...]
    # File[] bais = [...]
    #
    # We can handle this by zipping (for two), or actually transposing the array of both items, eg:
    # transpose([bam1, bam2, ..., bamn], [bai1, bai2, ..., bai3]) => [[bam1, bai1], [bam2, bai2], ..., [bamn, bain]]
    #
    # Providing we're delicate, we can zip these in a particlar way, and then unwrap them using their indices
    # and hoefully have everything line up:
    #
    # Source: https://software.broadinstitute.org/wdl/documentation/spec#arrayarrayx-transposearrayarrayx
    #

    inputs_map = {}
    for k in ins:
        inp = ins[k]
        if k not in node.connection_map:
            if inp.input_type.optional:
                continue
            else:
                raise Exception(f"Error when building connections for step '{node.id()}', "
                                f"could not find required connection: '{k}'")

        edge: StepInput = node.connection_map[k]
        source: Edge = edge.source()  # potentially single item or array

        if isinstance(source, list):
            if len(source) == 1:
                source = source[0]
            elif len(source) > 1:
                input_name_maps = ", ".join(edge.dotted_source())

                Logger.warn(f"Conversion to WDL for field '{node.id()}.{k}' does not fully support multiple sources."
                            f" This will only work if all of the inputs ({input_name_maps}) have the same secondaries "
                            f"AND this field ('{k}') is not scattered")

                unique_types = set((first_value(x.start.outputs()) if not x.stag else x.start.outputs()[x.stag]).output_type.name() for x in edge.source())
                if len(unique_types) > 1:
                    Logger.warn(f"There is more than one type of unique types mapped to the field '{node.id()}.{k}', "
                                f"per previous warning this might cause issues at runtime")

                inputs_map[k] = "[" + ", ".join(edge.dotted_source()) + "]"
                f = edge.finish.inputs()[edge.ftag]
                secs = f.input_type.subtype().secondary_files() if isinstance(f.input_type, Array) \
                    else f.input_type.secondary_files()
                if secs:
                    for sec in secs:
                        inputs_map[get_secondary_tag_from_original_tag(k, sec)] = \
                            "[" + ", ".join(get_secondary_tag_from_original_tag(kk, sec) for kk in edge.dotted_source()) + "]"
                continue
                # source = source[0]
                # raise Exception("Conversion to WDL does not currently support multiple sources")
                # ds = f'[{", ".join(ds)}]'

        secondary = None
        if source and isinstance(source.finish, StepNode):

            it = source.finish.inputs()[source.ftag].input_type
            secondary = it.subtype().secondary_files() if isinstance(it, Array) else it.secondary_files()
            if secondary and isinstance(source.start, StepNode):

                ot = source.start.outputs()[source.stag].output_type

                sec_out = set(value_or_default(ot.subtype().secondary_files() if isinstance(ot, Array)
                                               else ot.secondary_files(), default=[]))
                sec_in = set(secondary)
                if not sec_out.issubset(sec_in):
                    raise Exception(f"An error occurred when connecting '{source.dotted_source()}' to "
                                    f"'{source.finish.id()}.{source.ftag}', there were secondary files in the final node "
                                    f"that weren't present in the source: {', '.join(sec_out.difference(sec_in))}")
        if not source:
            # edge but no source, probably a default
            if not edge.default:
                Logger.critical(f"Skipping connection to '{edge.finish}.{edge.ftag}' had no source or default, "
                                f"please raise an issue as investigation may be required")
                continue

            if isinstance(edge.default, bool):
                inputs_map[k] = "true" if edge.default else "false"
            elif isinstance(edge.default, str):
                inputs_map[k] = f'"{edge.default}"'
            else:
                inputs_map[k] = edge.default

        elif edge in scatterable and secondary:
            # We're ensured through inheritance and .receiveBy that secondary files will match.
            ds = source.dotted_source()
            Logger.log(f"Oh boii, we're gonna have some complicated scattering here with {len(secondary)} secondary file(s)")

            identifier = old_to_new_identifier[ds]
            inputs_map[k] = identifier + "[0]"
            for idx in range(len(secondary)):
                sec = secondary[idx]
                inputs_map[get_secondary_tag_from_original_tag(k, sec)] = f"{identifier}[{idx + 1}]"
        else:
            ds = source.dotted_source()
            if ds in old_to_new_identifier and old_to_new_identifier[ds]:
                # can't get here with secondary
                inputs_map[k] = old_to_new_identifier[ds]
            else:
                inputs_map[k] = ds
                if secondary:
                    for idx in range(len(secondary)):
                        sec = secondary[idx]
                        inputs_map[get_secondary_tag_from_original_tag(k, sec)] = get_secondary_tag_from_original_tag(
                            ds, sec)

    inputs_map.update(resource_overrides)

    call = wdl.WorkflowCall(step_identifier, step_alias, inputs_map)

    for s in scatterable:
        if not isinstance(s.finish, StepNode):
            raise Exception("An internal error has occured when generating scatterable input map")
        secondary = s.finish.step.tool().inputs_map()[s.ftag].input_type.secondary_files()
        if secondary:
            print("There are secondary files here")
            ds = s.dotted_source()
            transformed = f"transform([{ds}, {', '.join(get_secondary_tag_from_original_tag(ds, sec) for sec in secondary)}])"
            call = wdl.WorkflowScatter(old_to_new_identifier[s.dotted_source()], transformed, [call])

        else:
            call = wdl.WorkflowScatter(old_to_new_identifier[s.dotted_source()], s.dotted_source(), [call])

    return call


def build_wdl_resource_inputs_dict(wf, hints, prefix=None) -> Dict[str, Any]:
    from janis.workflow.workflow import Workflow

    # returns a list of key, value pairs
    steps = {}
    if not prefix:
        prefix = wf.id() + "."
    else:
        prefix += "_"

    for s in wf._steps:
        tool: Tool = s.step.tool()


        if isinstance(tool, CommandTool):
            tool_pre = prefix + tool.id() + "_"
            steps.update([
                (tool_pre + "runtime_memory", tool.memory(hints)),
                (tool_pre + "runtime_cpu", tool.cpus(hints)),
                (tool_pre + "runtime_disks", None)
            ])
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            steps.update(build_wdl_resource_inputs_dict(tool, hints, tool_pre))

    return steps


def build_wdl_resource_inputs(wf, prefix=None) -> List[wdl.Input]:
    from janis.workflow.workflow import Workflow

    # returns a list of key, value pairs
    inputs = []
    if not prefix:
        prefix = "" # wf.id() + "."
    else:
        prefix += "_"

    for s in wf._steps:
        tool: Tool = s.step.tool()

        if isinstance(tool, CommandTool):
            tool_pre = prefix + tool.id() + "_"
            inputs.extend([
                wdl.Input(wdl.WdlType.parse_type("String?"), tool_pre + "runtime_memory"),
                wdl.Input(wdl.WdlType.parse_type("Int?"), tool_pre + "runtime_cpu"),
                wdl.Input(wdl.WdlType.parse_type("String?"), tool_pre + "runtime_disks")
            ])
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            inputs.extend(build_wdl_resource_inputs(tool, tool_pre))

    return inputs

import os
from typing import List, Dict, Optional
import itertools

import wdlgen as wdl
from janis.workflow.step import StepNode

from janis.graph.stepinput import Edge, StepInput

from janis.types import InputSelector, WildcardSelector

from janis.types.common_data_types import Stdout, Array

from janis.utils.validators import Validators

from janis.tool.commandtool import CommandTool

from janis.utils import first_value, convert_expression_to_wdl

from janis.tool.tool import Tool

from janis.utils.logger import Logger


def dump_wdl(workflow, to_console=True, to_disk=False, with_docker=False, with_hints=False,
             with_resource_overrides=False, write_inputs_file=False):
    wf_wdl, inp_str, tool_dicts = translate_workflow(workflow, with_docker=with_docker)

    wf_str = wf_wdl.get_string()
    tls_strs = [(t + ".wdl", tool_dicts[t].get_string()) for t in tool_dicts]

    if to_console:
        print("=== WORKFLOW ===")
        print(wf_str)
        print("\n=== INPUTS ===")
        print(inp_str)
        print("\n=== TOOLS ===")
        [print(t[1]) for t in tls_strs]

    if to_disk:
        d = os.path.expanduser("~") + f"/Desktop/{workflow.id()}/cwl/"
        d_tools = d + "tools/"

        if not os.path.isdir(d):
            os.makedirs(d)
        if not os.path.isdir(d_tools):
            os.makedirs(d_tools)

        os.chdir(d)
        wf_filename = d + workflow.id() + ".cwl"
        with open(wf_filename, "w+") as cwl:
            Logger.log(f"Writing {workflow.id()}.cwl to disk")
            cwl.write(wf_str)
            # ruamel.yaml.dump(wf_dict, cwl, default_flow_style=False)
            Logger.log(f"Written {workflow.id()}.cwl to disk")

        if write_inputs_file:
            with open(d + workflow.id() + "-job.yml", "w+") as cwl:
                Logger.log(f"Writing {workflow.id()}-job.yml to disk")
                cwl.write(inp_str)
                # ruamel.yaml.dump(inp_dict, cwl, default_flow_style=False)
                Logger.log(f"Written {workflow.id()}-job.yml to disk")
        else:
            Logger.log("Skipping writing input (yaml) job file")

        # z = zipfile.ZipFile(d + "tools.zip", "w")
        for (tool_filename, tool) in tls_strs:
            with open(d_tools + tool_filename, "w+") as cwl:
                Logger.log(f"Writing {tool_filename} to disk")
                cwl.write(tool)
                Logger.log(f"Written {tool_filename} to disk")

        import subprocess

        Logger.info("Validing outputted CWL")

        cwltool_result = subprocess.run(["cwltool", "--validate", wf_filename])
        if cwltool_result.returncode == 0:
            Logger.info("Exported workflow is valid CWL.")
        else:
            Logger.critical(cwltool_result.stderr)

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
    return ".".join(split[:-min(leading, len(split)-1)]) + fixed_sec


def translate_workflow(wf, with_docker=True, is_nested_tool=False):
    """

    :param with_docker:
    :param is_nested_tool:
    :return:
    """
    from janis.workflow.workflow import Workflow

    # Notes:
    #       All wdlgen classes have a .get_string(**kwargs) function
    #       The wdlgen Workflow class requires a

    w = wdl.Workflow(wf.identifier)
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
            wf_wdl, _, wf_tools = translate_workflow(t, with_docker=with_docker, is_nested_tool=True)
            wtools[s.id()] = wf_wdl
            wtools.update(wf_tools)
        elif isinstance(t, CommandTool):
            wtools[t.id()] = translate_tool(t, with_docker=with_docker)

        w.calls.append(
            translate_step_node(s, tool_aliases[t.id().lower()].upper() + "." + t.id(), s.id())
        )

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


def translate_input_selector(selector: InputSelector):
    if not selector.input_to_select: raise Exception("No input was selected for input selector: " + str(selector))
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""
    return f"{pre}${{{selector.input_to_select}}}{suf}"


def translate_wildcard_selector(selector: WildcardSelector):
    if not selector.wildcard: raise Exception("No wildcard was selected for wildcard selector: " + str(selector))
    return f"glob({selector.wildcard})"


def translate_tool(tool, with_docker):
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
            ins.append(wdl.Input(wd, i.id()))
            ins.extend(
                wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s)) for s in i.input_type.secondary_files())

    for o in tool.outputs():
        outs.extend(translate_output_node(o, tool))

    command_ins = [

        wdl.Task.Command.CommandInput(
            name=i.id(),
            optional=i.input_type.optional,
            prefix=i.prefix,
            position=i.position,
            separate_value_from_prefix=i.separate_value_from_prefix,
            default=i.default if i.default else i.input_type.default()
        ) for i in tool.inputs()] if tool.inputs() else None

    command_args = None
    if tool.arguments():
        command_args = []
        for a in tool.arguments():
            if a.value is None:
                val = None
            elif callable(getattr(a.value, "wdl", None)):
                val = a.value.wdl()
            else:
                val = a.value
            command_args.append(wdl.Task.Command.CommandArgument(a.prefix, val, a.position))

    command = wdl.Task.Command(tool.base_command(), command_ins, command_args)

    r = wdl.Task.Runtime()
    if with_docker:
        r.add_docker(tool.docker())

    return wdl.Task(tool.id(), ins, outs, command, r)


def translate_output_node(o, tool) -> List[wdl.Output]:
    if isinstance(o.output_type, Stdout):
        expression = "stdout()"
        return [wdl.Output(o.output_type.wdl(), o.id(), expression)]

    elif isinstance(o.glob, InputSelector):
        expression = translate_input_selector(o.glob)
        return [wdl.Output(o.output_type.wdl(), o.id(), f'"{expression}"')] + \
               [wdl.Output(o.output_type.wdl(), get_secondary_tag_from_original_tag(o.id(), s), f'"{expression}{s}"')
                for s in o.output_type.secondary_files()]

    elif isinstance(o.glob, WildcardSelector):
        expression = translate_wildcard_selector(o.glob)
        if not isinstance(o.output_type, Array):
            Logger.warn(f"The command tool '{tool.id()}.{o.tag}' used a star-bind (*) glob to find the output, "
                        f"but the return type was not an array. For WDL, the first element will be used, "
                        f"ie: '{expression}[0]'")
        expression += "[0]"

    else:
        print("NON SELECTOR GLOB: " + o.glob)
        glob = convert_expression_to_wdl(o.glob)
        if glob is not None and "*" in glob:
            glob = f'glob({glob})'
            if not isinstance(o.output_type, Array):
                Logger.warn(f"The command tool '{tool.id()}.{o.tag}' used a star-bind (*) glob to find the output, "
                            f"but the return type was not an array. For WDL, the first element will be used, "
                            f"ie: '{glob}[0]'")
                glob = glob + "[0]"
        expression = glob
        wdl_type = wdl.WdlType.parse_type(o.output_type.wdl())
        return [wdl.Output(wdl_type, o.id(), expression)] + \
               [wdl.Output(wdl_type, get_secondary_tag_from_original_tag(o.id(), s), expression)
                for s in o.output_type.secondary_files()]


def translate_step_node(node, step_identifier: str, step_alias: str):
    import wdlgen as wdl

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
        source: Edge = edge.source() # potentially single item or array

        if isinstance(source, list):
            if len(source) == 1:
                source = source[0]
            elif len(source) > 1:
                raise Exception("Conversion to WDL does not currently support multiple sources")
                # ds = f'[{", ".join(ds)}]'

        ds = source.dotted_source()
        secondary = None
        if isinstance(source.finish, StepNode):
            secondary = source.finish.inputs()[source.ftag].input_type.secondary_files()
            if isinstance(source.start, StepNode):
                sec_out = set(source.start.outputs()[source.stag].output_type.secondary_files())
                sec_in = set(secondary)
                if sec_out.issubset(sec_in) > 0:
                    raise Exception(f"An error occurred when connecting '{source.dotted_source()}' to "
                                    f"'{source.finish.id()}.{source.ftag}', there were secondary files in the final node "
                                    f"that weren't present in the source: {', '.join(sec_out.difference(sec_in))}")

        if edge in scatterable and secondary:
            # We're ensured through inheritance and .receiveBy that secondary files will match.
            print(f"Oh boii, we're gonna have some complicated scattering here with {len(secondary)} secondary file(s)")

            identifier = old_to_new_identifier[ds]
            inputs_map[k] = identifier + "[0]"
            for idx in range(len(secondary)):
                sec = secondary[idx]
                inputs_map[get_secondary_tag_from_original_tag(k, sec)] = f"{identifier}[{idx+1}]"

        else:
            if ds in old_to_new_identifier and old_to_new_identifier[ds]:
                # can't get here with secondary
                inputs_map[k] = old_to_new_identifier[ds]
            elif edge.default is not None:
                if isinstance(edge.default, bool):
                    inputs_map[k] = "true" if edge.default else "false"
                elif isinstance(edge.default, str):
                    inputs_map[k] = f'"{edge.default}"'
                else:
                    inputs_map[k] = edge.default
            else:
                inputs_map[k] = ds
                for idx in range(len(secondary)):
                    sec = secondary[idx]
                    inputs_map[get_secondary_tag_from_original_tag(k, sec)] = get_secondary_tag_from_original_tag(ds, sec)

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

# def map_inputs_for_non_secondary_scatterable_field():


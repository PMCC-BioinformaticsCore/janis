import os
from typing import List, Dict

import wdlgen as wdl
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
        w.inputs.extend(wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s)) for s in i.input.data_type.secondary_files())

    # Convert self._outputs -> wdl.Output
    for o in wf._outputs:
        w.outputs.append(wdl.Output(
            o.output.data_type.wdl(),
            o.id(),
            "{a}.{b}".format(  # Generate fully qualified stepid.tag identifier (MUST be step node)
                a=first_value(first_value(o.connection_map).source_map).start.id(),
                b=first_value(first_value(o.connection_map).source_map).stag
            )))
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

    inp = {f"{wf.id()}.{i.id()}": i.input.wdl_input() for i in wf._inputs}

    return w, inp, wtools


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
            ins.extend(wdl.Input(wd, get_secondary_tag_from_original_tag(i.id(), s)) for s in i.input_type.secondary_files())

    for o in tool.outputs():

        if isinstance(o.output_type, Stdout):
            expression = "stdout()"
        else:

            glob = convert_expression_to_wdl(o.glob)
            if glob is not None and "*" in glob:
                glob = f'glob({glob})'
                if not isinstance(o.output_type, Array):
                    Logger.warn(f"The command tool '{tool.id()}.{o.tag}' used a star-bind (*) glob to find the output, "
                                f"but the return type was not an array. For WDL, the first element will be used, "
                                f"ie: '{glob}[0]'")
                    glob = glob + "[0]"
            expression = glob

        outs.append(wdl.Output(o.output_type.wdl(), o.id(), expression))
        # TODO: Come up with better expression, maybe an InputSelector wrapper as it's super common, and may
        #       avoid having to do the expression conversions (between JS and else), and may be easier to validate
        outs.extend(wdl.Input(o.output_type.wdl(), get_secondary_tag_from_original_tag(o.id(), s), expression)
                    for s in o.output_type.secondary_files())

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


def translate_step_node(node, step_identifier: str, step_alias: str):
    import wdlgen as wdl

    ins = node.inputs()

    # One step => One WorkflowCall. We need to traverse the edge list to see if there's a scatter
    # then we can build up the WorkflowCall / ScatterCall
    scatterable = [node.connection_map[k].dotted_source()
                   for k in node.inputs() if k in node.connection_map and node.connection_map[k].has_scatter()]

    # We need to replace the scatterable key(s) with some random variable, eg: for i in iterable:
    ordered_variable_identifiers = ["i", "j", "k", "x", "y", "z", "a", "b", "c", "ii", "jj", "kk", "xx", "yy", "zz"]
    new_to_old_identifier = {k.dotted_source(): k.dotted_source() for k in node.connection_map.values()
                             if not isinstance(k.dotted_source(), list)}

    # We'll wrap everything in the scatter block later, but let's replace the fields we need to scatter
    # with the new scatter variable (we'll try to guess one based on the fieldname. We might need to eventually
    # pass the workflow inputs to make sure now conflict will arise.
    # Todo: Pass Workflow input tags to wdl scatter generation to ensure scatter var doesn't conflict with inputs

    for s in scatterable:
        new_var = s.split(".")[-1][0].lower()
        while new_var in new_to_old_identifier:
            new_var = ordered_variable_identifiers.pop(0)
        new_to_old_identifier[s] = new_var

    # Let's map the inputs, to the source.
    # We're using a dictionary for the map atm, but WDL requires the format:
    #       fieldName: sourceCall.Output

    inputs_map = {}
    for k in ins:
        inp = ins[k]
        if k not in node.connection_map:
            if inp.input_type.optional:
                continue
            else:
                raise Exception(f"Error when building connections for step '{node.id()}', "
                                f"could not find required connection: '{k}'")

        edge = node.connection_map[k]
        ds = edge.dotted_source()

        if isinstance(ds, list):
            if len(ds) == 1:
                ds = ds[0]
            elif len(ds) > 1:
                Logger.critical("Conversion to WDL does not currently support multiple sources")
                ds = f'[{", ".join(ds)}]'

        if ds in new_to_old_identifier and new_to_old_identifier[ds]:
            inputs_map[k] = new_to_old_identifier[ds]
        elif edge.default is not None:
            if isinstance(edge.default, bool):
                inputs_map[k] = "true" if edge.default else "false"
            elif isinstance(edge.default, str):
                inputs_map[k] = f'"{edge.default}"'
            else:
                inputs_map[k] = edge.default
        else:
            inputs_map[k] = ds

    call = wdl.WorkflowCall(step_identifier, step_alias, inputs_map)

    for s in scatterable:
        call = wdl.WorkflowScatter(new_to_old_identifier[s], s, [call])

    return call



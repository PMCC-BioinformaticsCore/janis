import os
import ruamel.yaml
import cwlgen
from typing import List, Dict, Type

from janis.hints import Hint, HintEnum, HintArray, HINTS
from janis.tool.tool import Tool, ToolInput
from janis.tool.commandtool import CommandTool
from janis.types import InputSelector, Selector, WildcardSelector
from janis.types.common_data_types import Int, Stdout, Array
from janis.utils.logger import Logger
from janis.utils.metadata import WorkflowMetadata, ToolMetadata
from janis.utils.janisconstants import RESOURCE_OVERRIDE_KEY, HINTS_KEY
from janis.workflow.step import StepNode


def dump_cwl(workflow, to_console=True, to_disk=False, with_docker=False, with_hints=False,
             with_resource_overrides=False, write_inputs_file=False):
    wf_cwl, inp_dict, tools_cwl = translate_workflow(workflow,
                                                     with_docker=with_docker,
                                                     with_hints=with_hints,
                                                     with_resource_overrides=with_resource_overrides)

    wf_dict = wf_cwl.get_dict()
    tool_dicts = [("tools/" + t.id + ".cwl", t.get_dict()) for t in tools_cwl]

    ruamel.yaml.add_representer(cwlgen.utils.literal, cwlgen.utils.literal_presenter)

    wf_str = ruamel.yaml.dump(wf_dict, default_flow_style=False)
    inp_str = ruamel.yaml.dump(inp_dict, default_flow_style=False)
    tls_strs = [(t[0], ruamel.yaml.dump(t[1], default_flow_style=False)) for t in tool_dicts]

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


def translate_workflow(wf, is_nested_tool=False, with_docker=False, with_hints=False,
                       with_resource_overrides=False):
    from janis.workflow.workflow import Workflow

    metadata = wf.metadata() if wf.metadata() else WorkflowMetadata()
    w = cwlgen.Workflow(wf.identifier, wf.friendly_name(), metadata.documentation)

    w.inputs: List[cwlgen.InputParameter] = [translate_input(i.input) for i in wf._inputs]

    if with_resource_overrides:
        rOverride = cwlgen.InputParameter(
            RESOURCE_OVERRIDE_KEY,
            param_type=["null", generate_cwl_resource_override_schema_for_steps(wf)]
        )
        w.inputs.append(rOverride)

    w.steps: List[cwlgen.WorkflowStep] = [translate_step(s, is_nested_tool=is_nested_tool) for s in wf._steps]

    w.outputs = [translate_output_node(o) for o in wf._outputs]

    #
    keys = ["coresMin", "coresMax", "ramMin", "ramMax"]
    sins = [translate_tool_input(ToolInput(k, Int(optional=True))) for k in keys]
    if with_resource_overrides:
        for s in w.steps:
            # work out whether (the tool of) s is a workflow or tool
            resource_override_step_inputs = [cwlgen.WorkflowStepInput(
                input_id=k,
                source=RESOURCE_OVERRIDE_KEY,
                value_from=f"${{var k = \"{k}\";var stepId = \"{s.id}\";if(!self) return null;if (!(stepId in self)) "
                "return null;return self[stepId][k]}}"
            ) for k in keys]
            s.inputs.extend(resource_override_step_inputs)

    if with_hints:
        resource_schema = get_cwl_schema_for_recognised_hints()
        nullable_resource_schema = ["null", resource_schema]
        w.inputs.append(cwlgen.InputParameter(HINTS_KEY, param_type=nullable_resource_schema))
        for s in w.steps:
            s.inputs.append(cwlgen.WorkflowStepInput(HINTS_KEY, HINTS_KEY))

    w.requirements.append(cwlgen.InlineJavascriptReq())
    w.requirements.append(cwlgen.StepInputExpressionRequirement())

    if wf.has_scatter:
        w.requirements.append(cwlgen.ScatterFeatureRequirement())
    if wf.has_subworkflow:
        w.requirements.append(cwlgen.SubworkflowFeatureRequirement())
    if wf.has_multiple_inputs:
        w.requirements.append(cwlgen.MultipleInputFeatureRequirement())

    tools = []
    tools_to_build: Dict[str, Tool] = {s.step.tool().id(): s.step.tool() for s in wf._steps}
    for t in tools_to_build:
        tool: Tool = tools_to_build[t]
        if isinstance(tool, Workflow):
            wf_cwl, _, subtools = translate_workflow(tool,
                                                     is_nested_tool=True,
                                                     with_docker=with_docker,
                                                     with_hints=with_hints,
                                                     with_resource_overrides=with_resource_overrides)
            tools.append(wf_cwl)
            tools.extend(subtools)
        elif isinstance(tool, CommandTool):
            tool_cwl = translate_tool(tool, with_docker=with_docker)  # tool.cwl(with_docker=with_docker)
            if with_hints:
                tool_cwl.inputs.append(
                    cwlgen.InputParameter(HINTS_KEY, param_type=["null", get_cwl_schema_for_recognised_hints()]))
                hm = tool.hint_map()
                if hm:
                    tool_cwl.requirements.append(cwlgen.ResourceRequirement(
                        cores_min=generate_hint_selectors_for_hint_map("coresMin", hm),
                        cores_max=generate_hint_selectors_for_hint_map("coresMax", hm),
                        ram_min=generate_hint_selectors_for_hint_map("ramMin", hm),
                        ram_max=generate_hint_selectors_for_hint_map("ramMax", hm)
                    ))

            if with_resource_overrides:
                tool_cwl.inputs.extend(sins)

            tools.append(tool_cwl)
        else:
            raise Exception(f"Unknown tool type: '{type(tool)}'")

    inp = {i.id(): i.input.cwl_input() for i in wf._inputs}

    return w, inp, tools


def translate_tool_str(tool, with_docker):
    ruamel.yaml.add_representer(cwlgen.utils.literal, cwlgen.utils.literal_presenter)
    return ruamel.yaml.dump(translate_tool(tool, with_docker=with_docker).get_dict(), default_flow_style=False)


def translate_tool(tool, with_docker):
    metadata = tool.metadata() if tool.metadata() else ToolMetadata()
    stdouts = [o.output_type for o in tool.outputs() if isinstance(o.output_type, Stdout) and o.output_type.stdoutname]
    stdout = stdouts[0].stdoutname if len(stdouts) > 0 else None

    tool_cwl = cwlgen.CommandLineTool(
        tool_id=tool.id(),
        base_command=tool.base_command(),
        label=tool.id(),
        doc=metadata.documentation,
        # cwl_version=Cwl.kCUR_VERSION,
        stdin=None,
        stderr=None,
        stdout=stdout
    )

    tool_cwl.requirements.extend([
        cwlgen.InlineJavascriptReq()
    ])

    if tool.requirements():
        tool_cwl.requirements.extend(tool.requirements())

    if with_docker:
        tool_cwl.requirements.append(cwlgen.DockerRequirement(
            docker_pull=tool.docker(),
            # docker_load=None,
            # docker_file=None,
            # docker_import=None,
            # docker_image_id=None,
            # docker_output_dir=None
        ))

    tool_cwl.inputs.extend(translate_tool_input(i) for i in tool.inputs())
    tool_cwl.outputs.extend(translate_tool_output(o, tool=tool.id()) for o in tool.outputs())
    args = tool.arguments()
    if args:
        tool_cwl.arguments.extend(translate_tool_argument(a) for a in tool.arguments())

    return tool_cwl


def translate_input(inp):
    return cwlgen.InputParameter(
        param_id=inp.id(),
        label=inp.label,
        secondary_files=inp.data_type.secondary_files(),
        param_format=None,
        streamable=None,
        doc=inp.doc,
        input_binding=None,
        param_type=inp.data_type.cwl_type()
    )


def translate_output_node(node):
    return translate_output(node.output, next(iter(node.connection_map.values())).source())


def translate_output(outp, source):
    return cwlgen.WorkflowOutputParameter(
        param_id=outp.id(),
        output_source=source,
        label=outp.label,
        secondary_files=outp.data_type.secondary_files(),
        param_format=None,
        streamable=None,
        doc=outp.doc,
        param_type=outp.data_type.cwl_type(),
        output_binding=None,
        linkMerge=None
    )


def translate_tool_input(toolinput):
    default = toolinput.default if toolinput.default else toolinput.input_type.default()

    data_type = toolinput.input_type.cwl_type()
    input_binding = cwlgen.CommandLineBinding(
        # load_contents=toolinput.load_contents,
        position=toolinput.position,
        prefix=toolinput.prefix,
        separate=toolinput.separate_value_from_prefix,
        # item_separator=toolinput.item_separator,
        # value_from=toolinput.value_from,
        shell_quote=toolinput.shell_quote,
    )

    # Binding array inputs onto the console
    # https://www.commonwl.org/user_guide/09-array-inputs/
    if isinstance(toolinput.input_type, Array) and isinstance(data_type, cwlgen.CommandInputArraySchema):
        if toolinput.nest_input_binding_on_array:
            input_binding.prefix = None
            input_binding.separate = None
            nested_binding = cwlgen.CommandLineBinding(
                # load_contents=toolinput.load_contents,
                prefix=toolinput.prefix,
                separate=toolinput.separate_value_from_prefix,
                # item_separator=toolinput.item_separator,
                # value_from=toolinput.value_from,
                shell_quote=toolinput.shell_quote,
            )
            data_type.inputBinding = nested_binding
        else:
            input_binding.itemSeparator = ","

    return cwlgen.CommandInputParameter(
        param_id=toolinput.tag,
        label=toolinput.tag,
        secondary_files=toolinput.input_type.secondary_files(),
        # streamable=None,
        doc=toolinput.doc,
        input_binding=input_binding,
        default=default,
        param_type=data_type
    )


def translate_input_selector(selector: InputSelector):
    if not selector.input_to_select: raise Exception("No input was selected for input selector: " + str(selector))
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""
    return f"{pre}$(inputs.{selector.input_to_select}){suf}"


def translate_tool_argument(argument):
    if argument.value is None:
        val = None
    elif callable(getattr(argument.value, "cwl", None)):
        val = argument.value.cwl()
    elif isinstance(argument.value, InputSelector):
        val = translate_input_selector(argument.value)
    else:
        val = str(argument.value)
    return cwlgen.CommandLineBinding(
        # load_contents=False,
        position=argument.position,
        prefix=argument.prefix,
        separate=argument.separate_value_from_prefix,
        # item_separator=None,
        value_from=val,
        shell_quote=argument.shell_quote,
    )


def translate_tool_output(output, **debugkwargs):
    return cwlgen.CommandOutputParameter(
        param_id=output.tag,
        label=output.tag,
        secondary_files=output.output_type.secondary_files(),
        # param_format=None,
        # streamable=None,
        doc=output.doc,
        output_binding=cwlgen.CommandOutputBinding(
            glob=translate_to_cwl_glob(output.glob, outputtag=output.tag, **debugkwargs),
            # load_contents=False,
            # output_eval=None
        ),
        param_type=output.output_type.cwl_type()
    )


def translate_to_cwl_glob(glob, **debugkwargs):
    if not glob: return None

    if not isinstance(glob, Selector):
        Logger.warn("String globs are being phased out from tool output selections, please use the provided "
                    "Selector (InputSelector or WildcardSelector) classes. " + str(debugkwargs))
        return glob

    if isinstance(glob, InputSelector):
        return translate_input_selector(glob)

    elif isinstance(glob, WildcardSelector):
        return glob.wildcard

    raise Exception("Unimplemented selector type: " + glob.__class__.__name__)


def translate_step(step: StepNode, is_nested_tool=False):
    run_ref = ("{tool}.cwl" if is_nested_tool else "tools/{tool}.cwl").format(tool=step.step.tool().id())
    cwlstep = cwlgen.WorkflowStep(
        step_id=step.id(),
        run=run_ref,
        label=step.step.label,
        doc=step.step.doc,
        scatter=None,  # Filled by StepNode
        scatter_method=None  # Filled by StepNode
    )

    cwlstep.out = [cwlgen.WorkflowStepOutput(output_id=o.tag) for o in step.step.tool().outputs()]

    ins = step.inputs()
    scatterable = []

    for k in ins:
        inp = ins[k]
        if k not in step.connection_map:
            if inp.input_type.optional:
                continue
            else:
                raise Exception(f"Error when building connections for cwlstep '{step.id()}', "
                                f"could not find required connection: '{k}'")

        inp_t = step.inputs()[k].input_type
        edge = step.connection_map[k]
        default = edge.default if edge.default else inp_t.default()
        d = cwlgen.WorkflowStepInput(
            input_id=inp.tag,
            source=edge.source(),
            link_merge=None,  # this will need to change when edges have multiple source_map
            default=default,
            value_from=None
        )
        if edge.has_scatter():
            scatterable.append(k)

        cwlstep.inputs.append(d)

    if len(scatterable) > 0:
        if len(scatterable) > 1:
            Logger.info("Discovered more than one scatterable field on cwlstep '{step_id}', "
                        "deciding scatterMethod to be dot_product".format(step_id=step.id()))
            cwlstep.scatterMethod = "dot_product"
        cwlstep.scatter = scatterable

    return cwlstep


def generate_cwl_resource_override_schema():
    schema = cwlgen.CommandInputRecordSchema()
    schema.fields = [
        cwlgen.CommandInputRecordSchema.CommandInputRecordField("coresMin", ["long", "string", "null"]),
        cwlgen.CommandInputRecordSchema.CommandInputRecordField("coresMax", ["int", "string", "null"]),
        cwlgen.CommandInputRecordSchema.CommandInputRecordField("ramMin", ["long", "string", "null"]),
        cwlgen.CommandInputRecordSchema.CommandInputRecordField("ramMax", ["int", "string", "null"])
    ]
    return schema


def generate_cwl_resource_override_schema_for_steps(wf):
    from janis.workflow.workflow import Workflow
    schema = cwlgen.CommandInputRecordSchema()

    for step in wf._steps:
        tool = step.step.tool()
        if isinstance(tool, Workflow):
            key = tool.id() + "_resource_override"  # self.RESOURCE_OVERRIDE_KEY
            override_schema = generate_cwl_resource_override_schema_for_steps(tool)
        else:
            key = step.step.id()
            override_schema = generate_cwl_resource_override_schema()

        schema.fields.append(cwlgen.CommandInputRecordSchema.CommandInputRecordField(key, ["null", override_schema]))

    return schema


def generate_hint_selectors_for_hint_map(resource_key, hint_map):
    import json
    return f"""${{
var key = '{resource_key}';
if (inputs["{RESOURCE_OVERRIDE_KEY}"] && inputs["{RESOURCE_OVERRIDE_KEY}"][key])
    return inputs["{RESOURCE_OVERRIDE_KEY}"][key];
var hints = inputs.hints;
if (!hints) return null;
var hintMap = {json.dumps(hint_map)};
for (var hint in hintMap) {{
    var providedHintValue = hints[hint];
    if (!providedHintValue || !(providedHintValue in hintMap[hint])) continue;
    var hintValueToResourceMap = hintMap[hint][providedHintValue];
    if (hintValueToResourceMap[key]) return hintValueToResourceMap[key];
}}
return null;
}}"""  # .replace("    ", "").replace("\n", "")


def get_cwl_schema_for_recognised_hints():
    schema = cwlgen.CommandInputRecordSchema("hints")

    def prepare_hint(hint_class: Type[Hint]):
        name = hint_class.key() + "_schema"

        if issubclass(hint_class, HintEnum):
            # assume is an enum
            return cwlgen.CommandInputEnumSchema(label=hint_class.key(), name=name, symbols=hint_class.symbols())
        elif issubclass(hint_class, HintArray):
            return cwlgen.CommandInputArraySchema(items=hint_class.items(), label=hint_class.key())
        else:
            return "string?"

    schema.fields = [
        cwlgen.CommandInputRecordSchema.CommandInputRecordField(hint.key(), ["null", prepare_hint(hint)])
        for hint in HINTS]

    return schema

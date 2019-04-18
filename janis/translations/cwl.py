"""
CWL

This is one of the more complicated classes, it takes the janis in-memory representation of a workflow,
and converts it into the equivalent CWL objects. Janis was built alongside testing for CWL, so a lot of
the concepts directly or pretty closely match. There are a few extra things that Janis has that need to
be mapped back.

This file is logically structured similar to the WDL equiv:

- Imports
- dump_cwl
- translate_workflow
- translate_tool (command tool)
- other translate methods
- selector helpers (InputSelector, WildcardSelector, CpuSelector, MemorySelector)
- helper methods
"""

## IMPORTS

import os
from typing import List, Dict, Optional, Any

import cwlgen
import ruamel.yaml

from janis.tool.commandtool import CommandTool
from janis.tool.tool import Tool, ToolInput
from janis.types import InputSelector, Selector, WildcardSelector, MemorySelector, CpuSelector
from janis.types.common_data_types import Stdout, Array, File, Filename
from janis.utils.logger import Logger
from janis.utils.metadata import WorkflowMetadata, ToolMetadata
from janis.translations.exportpath import ExportPathKeywords
from janis.workflow.input import Input
from janis.workflow.step import StepNode

CWL_VERSION = "v1.0"


## TRANSLATION

def dump_cwl(workflow, to_console=True, with_docker=True, with_resource_overrides=False, to_disk=False,
             export_path=ExportPathKeywords.default, write_inputs_file=False, should_validate=False,
             should_zip=True):
    wf_cwl, inp_dict, tools_cwl = translate_workflow(workflow,
                                                     with_docker=with_docker,
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

    d = ExportPathKeywords.resolve(export_path, workflow_spec="cwl", workflow_name=workflow.id())

    if write_inputs_file:
        with open(d + workflow.id() + "-job.yml", "w+") as cwl:
            Logger.log(f"Writing {workflow.id()}-job.yml to disk")
            cwl.write(inp_str)
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
        wf_filename = d + workflow.id() + ".cwl"
        with open(wf_filename, "w+") as cwl:
            Logger.log(f"Writing {workflow.id()}.cwl to disk")
            cwl.write(wf_str)
            # ruamel.yaml.dump(wf_dict, cwl, default_flow_style=False)
            Logger.log(f"Written {workflow.id()}.cwl to disk")

        # z = zipfile.ZipFile(d + "tools.zip", "w")
        for (tool_filename, tool) in tls_strs:
            with open(d + tool_filename, "w+") as cwl:
                Logger.log(f"Writing {tool_filename} to disk")
                cwl.write(tool)
                Logger.log(f"Written {tool_filename} to disk")

        import subprocess
        if should_validate:
            Logger.info("Validing outputted CWL")

            cwltool_result = subprocess.run(["cwltool", "--validate", wf_filename])
            if cwltool_result.returncode == 0:
                Logger.info("Exported workflow is valid CWL.")
            else:
                Logger.critical(cwltool_result.stderr)

        if should_zip:
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
    w = cwlgen.Workflow(wf.identifier, wf.friendly_name(), metadata.documentation, cwl_version=CWL_VERSION)

    w.inputs: List[cwlgen.InputParameter] = [translate_input(i.input) for i in wf._inputs]

    resource_inputs = []
    if with_resource_overrides:
        resource_inputs = build_cwl_resource_inputs(wf)
        w.inputs.extend(resource_inputs)

    w.steps: List[cwlgen.WorkflowStep] = []

    for s in wf._steps:
        resource_overrides = {}
        for r in resource_inputs:
            if not r.id.startswith(s.id()): continue

            resource_overrides[r.id[(len(s.id()) + 1):]] = r.id
        w.steps.append(translate_step(s, is_nested_tool=is_nested_tool, resource_overrides=resource_overrides))

    w.outputs = [translate_output_node(o) for o in wf._outputs]

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
            tool_cwl = translate_tool(tool, with_docker=with_docker, with_resource_overrides=with_resource_overrides)

            tools.append(tool_cwl)
        else:
            raise Exception(f"Unknown tool type: '{type(tool)}'")

    inp = {i.id(): i.input.cwl_input() for i in wf._inputs}

    return w, inp, tools


def translate_tool_str(tool, with_docker, with_resource_overrides=False):
    ruamel.yaml.add_representer(cwlgen.utils.literal, cwlgen.utils.literal_presenter)
    return ruamel.yaml.dump(
        translate_tool(tool, with_docker=with_docker, with_resource_overrides=with_resource_overrides)
        .get_dict(), default_flow_style=False)


def translate_tool(tool, with_docker, with_resource_overrides=False):
    metadata = tool.metadata() if tool.metadata() else ToolMetadata()
    stdouts = [o.output_type for o in tool.outputs() if isinstance(o.output_type, Stdout) and o.output_type.stdoutname]
    stdout = stdouts[0].stdoutname if len(stdouts) > 0 else None

    if isinstance(stdout, InputSelector): stdout = translate_input_selector(stdout)

    tool_cwl = cwlgen.CommandLineTool(
        tool_id=tool.id(),
        base_command=tool.base_command(),
        label=tool.id(),
        doc=metadata.documentation,
        cwl_version=CWL_VERSION,
        stdin=None,
        stderr=None,
        stdout=stdout
    )

    tool_cwl.requirements.extend([
        cwlgen.InlineJavascriptReq()
    ])

    if tool.requirements():
        tool_cwl.requirements.extend(tool.requirements())

    inputs_that_require_localisation = [ti for ti in tool.inputs()
                                        if ti.localise_file and (issubclass(type(ti.input_type), File)
                                                                 or (issubclass(type(ti.input_type), Array))
                                                                 and issubclass(type(ti.input_type.subtype()), File))]
    if inputs_that_require_localisation:
        tool_cwl.requirements.append(cwlgen.InitialWorkDirRequirement([
            "$(inputs.%s)" % ti.id() for ti in inputs_that_require_localisation]))

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

    if with_resource_overrides:
        # work out whether (the tool of) s is a workflow or tool
        tool_cwl.inputs.extend([
            cwlgen.CommandInputParameter("runtime_memory", param_type="float?"),
            cwlgen.CommandInputParameter("runtime_cpu", param_type="int?"),
            # cwlgen.CommandInputParameter("runtime_disks", param_type="string?"),
        ])

        tool_cwl.requirements.append(cwlgen.ResourceRequirement(
            cores_min="$(inputs.runtime_cpu ? inputs.runtime_cpu : 1)",
            ram_min="$(inputs.runtime_memory ? Math.floor(1024 * inputs.runtime_memory) : 4096)",
        ))

    return tool_cwl


def translate_input(inp: Input):
    return cwlgen.InputParameter(
        param_id=inp.id(),
        default=inp.default,
        label=inp.label,
        secondary_files=inp.data_type.secondary_files(),
        param_format=None,
        streamable=None,
        doc=inp.doc,
        input_binding=None,
        param_type=inp.data_type.cwl_type()
    )


def translate_output_node(node):
    return translate_output(node.output, next(iter(node.connection_map.values())).slashed_source())


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


def translate_tool_input(toolinput: ToolInput) -> cwlgen.CommandInputParameter:
    data_type = toolinput.input_type.cwl_type()

    default, value_from = toolinput.default, None

    if isinstance(toolinput.input_type, Filename):
        default = toolinput.input_type.generated_filename()
    elif is_selector(default):
        default = None
        value_from = get_input_value_from_potential_selector_or_generator(toolinput.default, toolinput.id())

    input_binding = cwlgen.CommandLineBinding(
        # load_contents=toolinput.load_contents,
        position=toolinput.position,
        prefix=toolinput.prefix,
        separate=toolinput.separate_value_from_prefix,
        item_separator=toolinput.separator,
        value_from=value_from,
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


def translate_tool_argument(argument):
    return cwlgen.CommandLineBinding(
        # load_contents=False,
        position=argument.position,
        prefix=argument.prefix,
        separate=argument.separate_value_from_prefix,
        # item_separator=None,
        value_from=get_input_value_from_potential_selector_or_generator(argument.value, argument.value),
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


def translate_step(step: StepNode, is_nested_tool=False, resource_overrides=Dict[str, str]):
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

        edge = step.connection_map[k]
        d = cwlgen.WorkflowStepInput(
            input_id=inp.tag,
            source=edge.slashed_source(),
            link_merge=None,  # this will need to change when edges have multiple source_map
            value_from=None
        )
        if edge.has_scatter():
            scatterable.append(k)

        cwlstep.inputs.append(d)

    for r in resource_overrides:
        cwlstep.inputs.append(cwlgen.WorkflowStepInput(input_id=r, source=resource_overrides[r]))

    if len(scatterable) > 0:
        if len(scatterable) > 1:
            Logger.info("Discovered more than one scatterable field on cwlstep '{step_id}', "
                        "deciding scatterMethod to be dot_product".format(step_id=step.id()))
            cwlstep.scatterMethod = "dot_product"
        cwlstep.scatter = scatterable

    return cwlstep


## SELECTORS


def is_selector(selector):
    return issubclass(type(selector), Selector)


def get_input_value_from_potential_selector_or_generator(value, tool_id, string_environment=True):
    if value is None:
        return None
    if isinstance(value, str):
        return value if string_environment else f'"{value}"'
    elif isinstance(value, int) or isinstance(value, float):
        return value
    elif isinstance(value, Filename):
        return value.generated_filename()
    elif isinstance(value, InputSelector):
        return translate_input_selector(selector=value)
    elif isinstance(value, WildcardSelector):
        raise Exception(f"A wildcard selector cannot be used as an argument value for '{tool_id}'")
    elif isinstance(value, CpuSelector):
        return translate_cpu_selector(value)
    elif isinstance(value, MemorySelector):
        return translate_memory_selector(value)
    elif callable(getattr(value, "cwl", None)):
        return value.cwl()

    raise Exception("Could not detect type %s to convert to input value" % type(value))


def translate_input_selector(selector: InputSelector):
    if not selector.input_to_select: raise Exception("No input was selected for input selector: " + str(selector))
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""
    basename_extra = ".basename" if selector.use_basename else ""
    return f"{pre}$(inputs.{selector.input_to_select}{basename_extra}){suf}"


def translate_to_cwl_glob(glob, **debugkwargs):
    if not glob: return None

    if not isinstance(glob, Selector):
        Logger.critical("String globs are being phased out from tool output selections, please use the provided "
                        "Selector (InputSelector or WildcardSelector) classes. " + str(debugkwargs))
        return glob

    if isinstance(glob, InputSelector):
        return translate_input_selector(glob)

    elif isinstance(glob, WildcardSelector):
        return glob.wildcard

    raise Exception("Unimplemented selector type: " + glob.__class__.__name__)


def translate_cpu_selector(selector: CpuSelector):
    return "$(inputs.runtime_cpu)"


def translate_memory_selector(selector: MemorySelector):
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""

    val = "$(Math.floor(runtime_memory))"

    pref = ('"%s" + ' % pre) if pre else ""
    suff = (' + "%s"' % suf) if suf else ""
    return pref + val + suff


## OTHER HELPERS

def build_cwl_resource_inputs_dict(wf, hints, prefix=None) -> Dict[str, Any]:
    from janis.workflow.workflow import Workflow

    # returns a list of key, value pairs
    steps: Dict[str, Optional[Any]] = {}

    if not prefix:
        prefix = ""
    else:
        prefix += "_"

    for s in wf._steps:
        tool: Tool = s.step.tool()

        if isinstance(tool, CommandTool):
            tool_pre = prefix + s.id() + "_"
            steps.update({
                tool_pre + "runtime_memory": tool.memory(hints),
                tool_pre + "runtime_cpu": tool.cpus(hints),
                # tool_pre + "runtime_disks": None
            })
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            steps.update(build_cwl_resource_inputs_dict(tool, hints, tool_pre))

    return steps


def build_cwl_resource_inputs(wf, prefix=None) -> List[cwlgen.InputParameter]:
    from janis.workflow.workflow import Workflow

    # returns a list of key, value pairs
    inputs = []
    if not prefix:
        prefix = ""  # wf.id() + "."
    else:
        prefix += "_"

    for s in wf._steps:
        tool: Tool = s.step.tool()

        if isinstance(tool, CommandTool):
            tool_pre = prefix + s.id() + "_"
            inputs.extend([
                cwlgen.InputParameter(tool_pre + "runtime_memory", param_type="float?"),
                cwlgen.InputParameter(tool_pre + "runtime_cpu", param_type="int?"),
                # cwlgen.InputParameter(tool_pre + "runtime_disks", param_type="string?"),
            ])
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            inputs.extend(build_cwl_resource_inputs(tool, tool_pre))

    return inputs

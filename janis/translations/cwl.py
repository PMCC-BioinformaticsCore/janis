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
import re
from typing import List, Dict, Optional, Any, Tuple

import cwlgen
import ruamel.yaml

from janis.tool.commandtool import CommandTool
from janis.tool.tool import Tool, ToolInput
from janis.translations.translationbase import TranslatorBase
from janis.types import (
    InputSelector,
    Selector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
)
from janis.types.common_data_types import Stdout, Array, File, Filename
from janis.utils.logger import Logger
from janis.utils.metadata import WorkflowMetadata, ToolMetadata
from janis.translations.exportpath import ExportPathKeywords
from janis.workflow.input import Input
from janis.workflow.step import StepNode

CWL_VERSION = "v1.0"


## TRANSLATION


class CwlTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="cwl")
        ruamel.yaml.add_representer(
            cwlgen.utils.literal, cwlgen.utils.literal_presenter
        )

    @staticmethod
    def stringify_translated_workflow(wf):
        return ruamel.yaml.dump(wf.get_dict(), default_flow_style=False)

    @staticmethod
    def stringify_translated_tool(tool):
        return ruamel.yaml.dump(tool.get_dict(), default_flow_style=False)

    @staticmethod
    def stringify_translated_inputs(inputs):
        return ruamel.yaml.dump(inputs, default_flow_style=False)

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return ["cwltool", "--validate", wfpath]

    @classmethod
    def translate_workflow(
        cls, wf, with_docker=True, with_resource_overrides=False, is_nested_tool=False
    ) -> Tuple[any, Dict[str, any]]:
        from janis.workflow.workflow import Workflow

        metadata = wf.metadata() if wf.metadata() else WorkflowMetadata()
        w = cwlgen.Workflow(
            wf.identifier,
            wf.friendly_name(),
            metadata.documentation,
            cwl_version=CWL_VERSION,
        )

        w.inputs: List[cwlgen.InputParameter] = [
            translate_input(i.input) for i in wf._inputs
        ]

        resource_inputs = []
        if with_resource_overrides:
            resource_inputs = build_resource_override_maps_for_workflow(wf)
            w.inputs.extend(resource_inputs)

        w.steps: List[cwlgen.WorkflowStep] = []

        for s in wf._steps:
            resource_overrides = {}
            for r in resource_inputs:
                if not r.id.startswith(s.id()):
                    continue

                resource_overrides[r.id[(len(s.id()) + 1) :]] = r.id
            w.steps.append(
                translate_step(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                )
            )

        w.outputs = [translate_output_node(o) for o in wf._outputs]

        w.requirements.append(cwlgen.InlineJavascriptReq())
        w.requirements.append(cwlgen.StepInputExpressionRequirement())

        if wf.has_scatter:
            w.requirements.append(cwlgen.ScatterFeatureRequirement())
        if wf.has_subworkflow:
            w.requirements.append(cwlgen.SubworkflowFeatureRequirement())
        if wf.has_multiple_inputs:
            w.requirements.append(cwlgen.MultipleInputFeatureRequirement())

        tools = {}
        tools_to_build: Dict[str, Tool] = {
            s.step.tool().id(): s.step.tool() for s in wf._steps
        }
        for t in tools_to_build:
            tool: Tool = tools_to_build[t]
            if isinstance(tool, Workflow):
                wf_cwl, subtools = cls.translate_workflow(
                    tool,
                    is_nested_tool=True,
                    with_docker=with_docker,
                    with_resource_overrides=with_resource_overrides,
                )
                tools[tool.id()] = wf_cwl
                tools.update(subtools)
            elif isinstance(tool, CommandTool):
                tool_cwl = cls.translate_tool(
                    tool,
                    with_docker=with_docker,
                    with_resource_overrides=with_resource_overrides,
                )
                tools[tool.id()] = tool_cwl
            else:
                raise Exception(f"Unknown tool type: '{type(tool)}'")

        return w, tools

    @classmethod
    def build_inputs_file(
        cls,
        workflow,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
    ) -> Dict[str, any]:

        inp = {
            i.id(): i.input.cwl_input()
            for i in workflow._inputs
            if i.input.include_in_inputs_file_if_none or i.input.value
        }

        if merge_resources:
            inp.update(cls.build_resources_input(workflow, hints))

        if additional_inputs:
            for k in inp:
                if k in additional_inputs:
                    inp[k] = additional_inputs[k]

        return inp

    @classmethod
    def translate_workflow_to_all_in_one(
        cls, wf, with_resource_overrides=False, is_nested_tool=False
    ) -> cwlgen.Workflow:
        from janis.workflow.workflow import Workflow

        metadata = wf.metadata() if wf.metadata() else WorkflowMetadata()
        w = cwlgen.Workflow(
            wf.identifier,
            wf.friendly_name(),
            metadata.documentation,
            cwl_version=CWL_VERSION,
        )

        w.inputs: List[cwlgen.InputParameter] = [
            translate_input(i.input) for i in wf._inputs
        ]

        resource_inputs = []
        if with_resource_overrides:
            resource_inputs = build_resource_override_maps_for_workflow(wf)
            w.inputs.extend(resource_inputs)

        w.steps: List[cwlgen.WorkflowStep] = []

        for s in wf._steps:
            resource_overrides = {}
            for r in resource_inputs:
                if not r.id.startswith(s.id()):
                    continue

                resource_overrides[r.id[(len(s.id()) + 1) :]] = r.id

            w.steps.append(
                translate_step(
                    s,
                    is_nested_tool=is_nested_tool,
                    resource_overrides=resource_overrides,
                    use_run_ref=False,
                )
            )

        w.outputs = [translate_output_node(o) for o in wf._outputs]

        w.requirements.append(cwlgen.InlineJavascriptReq())
        w.requirements.append(cwlgen.StepInputExpressionRequirement())

        if wf.has_scatter:
            w.requirements.append(cwlgen.ScatterFeatureRequirement())
        if wf.has_subworkflow:
            w.requirements.append(cwlgen.SubworkflowFeatureRequirement())
        if wf.has_multiple_inputs:
            w.requirements.append(cwlgen.MultipleInputFeatureRequirement())

        return w

    @classmethod
    def translate_tool(cls, tool, with_docker=True, with_resource_overrides=False):
        metadata = tool.metadata() if tool.metadata() else ToolMetadata()
        stdouts = [
            o.output_type
            for o in tool.outputs()
            if isinstance(o.output_type, Stdout) and o.output_type.stdoutname
        ]
        stdout = stdouts[0].stdoutname if len(stdouts) > 0 else None

        if isinstance(stdout, InputSelector):
            stdout = translate_input_selector(stdout, code_environment=False)

        tool_cwl = cwlgen.CommandLineTool(
            tool_id=tool.id(),
            base_command=tool.base_command(),
            label=tool.id(),
            doc=metadata.documentation,
            cwl_version=CWL_VERSION,
            stdin=None,
            stderr=None,
            stdout=stdout,
        )

        if any(not i.shell_quote for i in tool.inputs()):
            tool_cwl.requirements.append(cwlgen.ShellCommandRequirement())

        tool_cwl.requirements.extend([cwlgen.InlineJavascriptReq()])

        if tool.requirements():
            tool_cwl.requirements.extend(tool.requirements())

        inputs_that_require_localisation = [
            ti
            for ti in tool.inputs()
            if ti.localise_file
            and (
                issubclass(type(ti.input_type), File)
                or (issubclass(type(ti.input_type), Array))
                and issubclass(type(ti.input_type.subtype()), File)
            )
        ]
        if inputs_that_require_localisation:
            tool_cwl.requirements.append(
                cwlgen.InitialWorkDirRequirement(
                    [
                        "$(inputs.%s)" % ti.id()
                        for ti in inputs_that_require_localisation
                    ]
                )
            )

        if with_docker:
            tool_cwl.requirements.append(
                cwlgen.DockerRequirement(
                    docker_pull=tool.docker(),
                    # docker_load=None,
                    # docker_file=None,
                    # docker_import=None,
                    # docker_image_id=None,
                    # docker_output_dir=None
                )
            )

        tool_cwl.inputs.extend(translate_tool_input(i) for i in tool.inputs())
        tool_cwl.outputs.extend(
            translate_tool_output(o, tool=tool.id()) for o in tool.outputs()
        )

        args = tool.arguments()
        if args:
            tool_cwl.arguments.extend(
                translate_tool_argument(a) for a in tool.arguments()
            )

        if with_resource_overrides:
            # work out whether (the tool of) s is a workflow or tool
            tool_cwl.inputs.extend(
                [
                    cwlgen.CommandInputParameter("runtime_memory", param_type="float?"),
                    cwlgen.CommandInputParameter("runtime_cpu", param_type="int?"),
                    # cwlgen.CommandInputParameter("runtime_disks", param_type="string?"),
                ]
            )

            tool_cwl.requirements.append(
                cwlgen.ResourceRequirement(
                    cores_min="$(inputs.runtime_cpu ? inputs.runtime_cpu : 1)",
                    ram_min="$(inputs.runtime_memory ? Math.floor(1024 * inputs.runtime_memory) : 4096)",
                )
            )

        return tool_cwl

    @classmethod
    def build_resources_input(cls, workflow, hints, prefix=None):
        from janis.workflow.workflow import Workflow

        # returns a list of key, value pairs
        steps: Dict[str, Optional[Any]] = {}

        if not prefix:
            prefix = ""
        else:
            prefix += "_"

        for s in workflow._steps:
            tool: Tool = s.step.tool()

            if isinstance(tool, CommandTool):
                tool_pre = prefix + s.id() + "_"
                steps.update(
                    {
                        tool_pre + "runtime_memory": tool.memory(hints),
                        tool_pre + "runtime_cpu": tool.cpus(hints),
                        # tool_pre + "runtime_disks": None
                    }
                )
            elif isinstance(tool, Workflow):
                tool_pre = prefix + s.id()
                steps.update(cls.build_resources_input(tool, hints, tool_pre))

        return steps

    @staticmethod
    def workflow_filename(workflow):
        return workflow.id() + ".cwl"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + "-inp.yml"

    @staticmethod
    def tool_filename(tool):
        return (tool.id() if isinstance(tool, Tool) else str(tool)) + ".cwl"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.yml"


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
        param_type=inp.data_type.cwl_type(inp.default is not None),
    )


def translate_output_node(node):
    return translate_output(
        node.output, next(iter(node.connection_map.values())).slashed_source()
    )


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
        linkMerge=None,
    )


def translate_tool_input(toolinput: ToolInput) -> cwlgen.CommandInputParameter:

    default, value_from = toolinput.default, None

    if isinstance(toolinput.input_type, Filename):
        default = toolinput.input_type.generated_filename()
        # value_from = get_input_value_from_potential_selector_or_generator(toolinput.input_type, code_environment=False, toolid=toolinput.id())
    elif is_selector(default):
        default = None
        value_from = get_input_value_from_potential_selector_or_generator(
            toolinput.default, code_environment=False, toolid=toolinput.id()
        )

    data_type = toolinput.input_type.cwl_type(default is not None)

    input_binding = cwlgen.CommandLineBinding(
        # load_contents=toolinput.load_contents,
        position=toolinput.position,
        prefix=toolinput.prefix,
        separate=toolinput.separate_value_from_prefix,
        item_separator=toolinput.separator,
        value_from=value_from,
        shell_quote=toolinput.shell_quote,
    )

    non_optional_dt_component = (
        [t for t in data_type if t != "null"][0]
        if isinstance(data_type, list)
        else data_type
    )

    # Binding array inputs onto the console
    # https://www.commonwl.org/user_guide/09-array-inputs/
    if isinstance(toolinput.input_type, Array) and isinstance(
        non_optional_dt_component, cwlgen.CommandInputArraySchema
    ):
        if toolinput.prefix_applies_to_all_elements:
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
            non_optional_dt_component.inputBinding = nested_binding

    return cwlgen.CommandInputParameter(
        param_id=toolinput.tag,
        label=toolinput.tag,
        secondary_files=toolinput.input_type.secondary_files(),
        # streamable=None,
        doc=toolinput.doc,
        input_binding=input_binding,
        default=default,
        param_type=data_type,
    )


def translate_tool_argument(argument):
    return cwlgen.CommandLineBinding(
        # load_contents=False,
        position=argument.position,
        prefix=argument.prefix,
        separate=argument.separate_value_from_prefix,
        # item_separator=None,
        value_from=get_input_value_from_potential_selector_or_generator(
            argument.value, code_environment=False
        ),
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
            glob=translate_to_cwl_glob(
                output.glob, outputtag=output.tag, **debugkwargs
            ),
            # load_contents=False,
            # output_eval=None
        ),
        param_type=output.output_type.cwl_type(),
    )


def translate_step(
    step: StepNode,
    is_nested_tool=False,
    resource_overrides=Dict[str, str],
    use_run_ref=True,
):

    tool = step.step.tool()
    if use_run_ref:
        run_ref = ("{tool}.cwl" if is_nested_tool else "tools/{tool}.cwl").format(
            tool=step.step.tool().id()
        )
    else:
        from janis.workflow.workflow import Workflow

        has_resources_overrides = len(resource_overrides) > 0
        if isinstance(tool, Workflow):
            run_ref = CwlTranslator.translate_workflow_to_all_in_one(
                tool, with_resource_overrides=has_resources_overrides
            )
        else:
            run_ref = CwlTranslator.translate_tool(
                tool, True, with_resource_overrides=has_resources_overrides
            )

    cwlstep = cwlgen.WorkflowStep(
        step_id=step.id(),
        run=run_ref,
        label=step.step.label,
        doc=step.step.doc,
        scatter=None,  # Filled by StepNode
        scatter_method=None,  # Filled by StepNode
    )

    cwlstep.out = [
        cwlgen.WorkflowStepOutput(output_id=o.tag) for o in step.step.tool().outputs()
    ]

    ins = step.inputs()
    scatterable = []

    for k in ins:
        inp = ins[k]
        if k not in step.connection_map:
            if inp.input_type.optional:
                continue
            else:
                raise Exception(
                    f"Error when building connections for cwlstep '{step.id()}', "
                    f"could not find required connection: '{k}'"
                )

        edge = step.connection_map[k]
        d = cwlgen.WorkflowStepInput(
            input_id=inp.tag,
            source=edge.slashed_source(),
            link_merge=None,  # this will need to change when edges have multiple source_map
            value_from=None,
        )
        if edge.has_scatter():
            scatterable.append(k)

        cwlstep.inputs.append(d)

    for r in resource_overrides:
        cwlstep.inputs.append(
            cwlgen.WorkflowStepInput(input_id=r, source=resource_overrides[r])
        )

    if len(scatterable) > 0:
        if len(scatterable) > 1:
            Logger.info(
                "Discovered more than one scatterable field on cwlstep '{step_id}', "
                "deciding scatterMethod to be dot_product".format(step_id=step.id())
            )
            cwlstep.scatterMethod = "dot_product"
        cwlstep.scatter = scatterable

    return cwlstep


## SELECTORS


def is_selector(selector):
    return issubclass(type(selector), Selector)


def get_input_value_from_potential_selector_or_generator(
    value, code_environment=True, **debugkwargs
):
    if value is None:
        return None
    if isinstance(value, str):
        return f'"{value}"' if code_environment else value
    elif isinstance(value, int) or isinstance(value, float):
        return value
    elif isinstance(value, Filename):
        # value.generated_filenamecwl() if code_environment else f"$({value.generated_filenamecwl()})"
        return (
            f'"{value.generated_filename()}"'
            if code_environment
            else value.generated_filename()
        )
    elif isinstance(value, StringFormatter):
        return translate_string_formatter(
            value, code_environment=code_environment, **debugkwargs
        )
    elif isinstance(value, InputSelector):
        return translate_input_selector(
            selector=value, code_environment=code_environment
        )
    elif isinstance(value, WildcardSelector):
        raise Exception(
            f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
        )
    elif isinstance(value, CpuSelector):
        return translate_cpu_selector(value)
    elif isinstance(value, MemorySelector):
        return translate_memory_selector(value)
    elif callable(getattr(value, "cwl", None)):
        return value.cwl()

    raise Exception("Could not detect type %s to convert to input value" % type(value))


def translate_input_selector(selector: InputSelector, code_environment):
    if not selector.input_to_select:
        raise Exception("No input was selected for input selector: " + str(selector))

    basename_extra = ".basename" if selector.use_basename else ""
    base = f"inputs.{selector.input_to_select}{basename_extra}"
    return base if code_environment else f"$({base})"


def translate_string_formatter(
    selector: StringFormatter, code_environment=True, **debugkwargs
):

    if len(selector.kwargs) == 0:
        return selector._format

    kwargreplacements = [
        f".replace(/{re.escape('{' +k + '}')}/g, {get_input_value_from_potential_selector_or_generator(v, code_environment=True, **debugkwargs)})"
        for k, v in selector.kwargs.items()
    ]
    return f'$("{selector._format}"' + "".join(kwargreplacements) + ")"


def translate_to_cwl_glob(glob, **debugkwargs):
    if not glob:
        return None

    if not isinstance(glob, Selector):
        Logger.critical(
            "String globs are being phased out from tool output selections, please use the provided "
            "Selector (InputSelector or WildcardSelector) classes. " + str(debugkwargs)
        )
        return glob

    if isinstance(glob, InputSelector):
        return translate_input_selector(glob, code_environment=False)

    elif isinstance(glob, StringFormatter):
        return translate_string_formatter(glob)

    elif isinstance(glob, WildcardSelector):
        return glob.wildcard

    raise Exception("Unimplemented selector type: " + glob.__class__.__name__)


def translate_cpu_selector(selector: CpuSelector):
    return "$(inputs.runtime_cpu)"


def translate_memory_selector(selector: MemorySelector):
    pre = selector.prefix if selector.prefix else ""
    suf = selector.suffix if selector.suffix else ""

    val = "$(Math.floor(inputs.runtime_memory))"

    pref = ('"%s" + ' % pre) if pre else ""
    suff = (' + "%s"' % suf) if suf else ""
    return pref + val + suff


## OTHER HELPERS


def build_resource_override_maps_for_workflow(
    wf, prefix=None
) -> List[cwlgen.InputParameter]:
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
            inputs.extend(
                [
                    cwlgen.InputParameter(
                        tool_pre + "runtime_memory", param_type="float?"
                    ),
                    cwlgen.InputParameter(tool_pre + "runtime_cpu", param_type="int?"),
                    # cwlgen.InputParameter(tool_pre + "runtime_disks", param_type="string?"),
                ]
            )
        elif isinstance(tool, Workflow):
            tool_pre = prefix + s.id()
            inputs.extend(build_resource_override_maps_for_workflow(tool, tool_pre))

    return inputs

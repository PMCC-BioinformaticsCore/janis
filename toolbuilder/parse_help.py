from typing import Optional, Union, List, Tuple

from janis_core.utils.validators import Validators

from .cltconvert import convert_command_tool_fragments
from janis_core import (
    ToolMetadata,
    String,
    Logger,
    JanisShed,
    DataType,
    Boolean,
    Filename,
)
from janis_core.tool.commandtool import ToolInput

from .templates import ToolTemplateType

container_exec = {
    "docker": ["docker", "run"],
    # "singularity": ["singularity", "exec"]
}

common_replacements = {
    "input": "inp",
    "output": "outputFilename",
}

option_markers = {
    "options:",
    "arguments:",
    "required arguments:",
    "optional arguments:",
}


def get_help_from_container(
    container: str,
    basecommand: Union[str, List[str]],
    help_param: Optional[str] = "--help",
    containersoftware="docker",
):
    import subprocess, os

    bc = basecommand if isinstance(basecommand, list) else [basecommand]
    cmd = [*container_exec[containersoftware], container, *bc]

    if help_param:
        cmd.append(help_param)

    print("Running command: " + " ".join(f"'{x}'" for x in cmd))
    try:
        help = (
            subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            .decode("utf-8")
            .rstrip()
        )
    except subprocess.CalledProcessError as e:
        help = str(e.output)
    return help


def get_version_from_container(
    container: str,
    basecommand: Union[str, List[str]],
    versionparam: Optional[str] = "--version",
    containersoftware="docker",
):
    import subprocess

    bc = basecommand if isinstance(basecommand, list) else [basecommand]
    cmd = [*container_exec[containersoftware], container, *bc, versionparam]

    print("Running command: " + " ".join(f"'{x}'" for x in cmd))
    try:
        help = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        return None

    return help.decode("utf-8").rstrip()


def first_or_default(iterable, default=None):
    filtered = [f for f in iterable if f is not None]
    if len(filtered) > 0:
        return filtered[0]
    return default


def parse_str(
    helpstr, option_marker: str = None, requires_prev_line_blank_or_param=False
):
    doc = ""
    args = []
    lines = helpstr.replace("\\n", "\n").split("\n")
    options_idx = None
    markers = option_markers
    if option_marker:
        markers = markers.union({option_marker.lower()})
    for il in range(len(lines)):
        line = lines[il]
        if not line.lstrip():
            continue

        ll = line.lower()

        if any(ll.startswith(m) for m in markers):
            options_idx = il
            break

        doc += line + "\n"

    if options_idx is None:
        raise Exception("Couldn't find the start of the inputs")

    prev_arg = None
    last_line_was_blank_or_param = True

    while options_idx < len(lines) - 1:
        options_idx += 1

        line = lines[options_idx]

        if not line.lstrip():
            # line is empty
            prev_arg = None
            last_line_was_blank_or_param = True
            continue

        line_args = [l.strip() for l in line.lstrip().split("  ") if l]
        largs = len(line_args)

        if largs == 0:
            raise Exception("No args when should have been filtered by previous step")

        tool_doc = ""

        if (
            not requires_prev_line_blank_or_param or last_line_was_blank_or_param
        ) and line_args[0].startswith("-"):
            # sometimes this section has two items
            processed_tags = [
                get_tag_and_cleanup_prefix(p) for p in line_args[0].split(",")
            ]
            tags = sorted(processed_tags, key=lambda l: len(l[1]), reverse=True,)
            potential_type = first_or_default([p[3] for p in processed_tags])

            if len(tags) > 1:
                tool_doc += "(" + ", ".join(t[0] for t in tags[1:]) + ") "

            if largs > 1:
                tool_doc += " ".join(line_args[1:])

            prefix, tag, has_equal, guessed_type = tags[0]
            eqifrequired = "=" if has_equal else ""

            if not potential_type:
                potential_type = Boolean

            if len(tag) == 1:
                while not Validators.validate_identifier(tag):
                    print(
                        f"The tag for '{prefix}' was invalid, we need you to come up with a new identifier for:"
                    )
                    print("\t" + tool_doc if tool_doc else line)
                    tag = str(input("New identifier: "))
            try:
                prev_arg = ToolInput(
                    tag,
                    potential_type(optional=True),
                    prefix=prefix + eqifrequired,
                    separate_value_from_prefix=not has_equal,
                    doc=tool_doc.replace('"', "'"),
                )
            except:
                print(f"Skipping '{tag}' as it wasn't validated correctly")
            args.append(prev_arg)
            # we'll get the longer one for the tag

        elif prev_arg:
            prev_arg.doc.doc += " " + line.lstrip()
        else:
            last_line_was_blank_or_param = False

    return doc, args


def guess_type(potential_type: str):
    if not potential_type:
        return None
    l = potential_type.lower()
    hopeful_type = JanisShed.get_datatype(l)

    if not hopeful_type:
        if "st" in potential_type:
            hopeful_type = String

    if hopeful_type:
        Logger.info(f"Found type {hopeful_type.__name__} from tag: {potential_type}")

    return hopeful_type


def get_tag_and_cleanup_prefix(prefix) -> Tuple[str, str, bool, Optional[DataType]]:
    """
    :param prefix:
    :return: (raw_element, potentialID, hasSeparator, potentialType)
    """
    # cases:
    # -a ADAPTER
    # --adapter=ADAPTER
    # --quality-cutoff=[5'CUTOFF,]3'CUTOFF
    el = prefix.lstrip()
    has_equals = False
    pretag = None
    potential_type = None

    if ":" in el:
        parts = el.split(":")
        if len(parts) > 2:
            Logger.warn(
                f"Unexpected number of components in the tag '{el}' to guess the type, using '{parts[0]}' and skipping type inference"
            )
        else:
            el, pt = parts[0], guess_type(parts[1])

            if not potential_type and pt:
                potential_type = pt

    if "=" in el:
        has_equals = True
        el = el.split("=")[0]
    elif " " in el:
        el = el.split(" ")[0]

    titleComponents = [l.strip().lower() for l in el.split("-") if l]
    if len(titleComponents) == 0:
        raise Exception(
            f"Title components for tag '{prefix}' does not have a component"
        )
    tag = "_".join(titleComponents)

    if tag.lower() in common_replacements:
        tag = common_replacements[tag.lower()]

    if tag.lower() == "outputfilename":
        potential_type = Filename

    return el, tag, has_equals, potential_type


def from_container(
    container: str,
    basecommand: List[str],
    helpcommand="-h",
    containersoftware="docker",
    optionsmarker: Optional[str] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
    type: ToolTemplateType = ToolTemplateType.base,
):
    helpstr = get_help_from_container(
        container=container,
        basecommand=basecommand,
        help_param=helpcommand,
        containersoftware=containersoftware,
    )
    tooldoc, args = parse_str(helpstr, option_marker=optionsmarker)

    if not version:
        comps = container.split(":")
        version = comps[-1] if len(comps) > 1 else "Latest"

    if not version:
        version = get_version_from_container(
            container,
            basecommand,
            versionparam="-v",
            containersoftware=containersoftware,
        )

    return (
        convert_command_tool_fragments(
            type=type,
            toolid=name or basecommand,
            basecommand=basecommand,
            friendly_name="".join(basecommand),
            toolprov="TOOLPROVIDER",
            ins=args,
            outs=[],
            metadata=ToolMetadata(documentation=tooldoc),
            container=container,
            version=version,
        ),
        helpstr,
    )

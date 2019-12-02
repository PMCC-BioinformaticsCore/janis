from typing import Optional, Union, List

from .cltconvert import convert_command_tool_fragments
from janis_core import ToolMetadata, String, Logger, JanisShed
from janis_core.tool.commandtool import ToolInput

container_exec = {
    "docker": ["docker", "run"],
    # "singularity": ["singularity", "exec"]
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


def parse_str(
    helpstr, option_marker: str = "Options:", requires_prev_line_blank_or_param=False
):
    doc = ""
    args = []
    lines = helpstr.replace("\\n", "\n").split("\n")
    options_idx = None
    for il in range(len(lines)):
        line = lines[il]
        if not line.lstrip():
            continue

        if line.startswith(option_marker):
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
            tags = sorted(
                [get_tag_and_cleanup_prefix(p) for p in line_args[0].split(",")],
                key=lambda l: len(l[1]),
                reverse=True,
            )

            if len(tags) > 1:
                tool_doc += "(" + ", ".join(t[0] for t in tags[1:]) + ") "

            if largs > 1:
                tool_doc += " ".join(line_args[1:])

            prefix, tag, has_equal = tags[0]
            eqifrequired = "=" if has_equal else ""

            datatype = String(optional=True)

            if ":" in tag:
                parts = tag.split(":")
                tag = parts[0]
                potentialtypestr = "".join(parts[1:])
                potentialtype = JanisShed.get_datatype(potentialtypestr)
                if potentialtype:
                    Logger.log(f"Found type {potentialtype.__name__} from tag: {tag}")
                    datatype = potentialtype(optional=True)

            if len(tag) == 1:
                print(
                    f"The tag for '{prefix}' was too short, we need you to come up with a new identifier for:"
                )
                print("\t" + tool_doc if tool_doc else line)
                tag = str(input("New identifier: "))
            try:
                prev_arg = ToolInput(
                    tag,
                    datatype,
                    prefix=prefix + eqifrequired,
                    separate_value_from_prefix=not has_equal,
                    doc=tool_doc.replace('"', "'"),
                )
            except:
                print(f"Skipping '{tag}' as it wasn't validated correctly")
            args.append(prev_arg)
            # we'll get the longer one for the tag

        elif prev_arg:
            prev_arg.doc += " " + line.lstrip()
        else:
            last_line_was_blank_or_param = False

    return doc, args


def get_tag_and_cleanup_prefix(prefix):
    # cases:
    # -a ADAPTER
    # --adapter=ADAPTER
    # --quality-cutoff=[5'CUTOFF,]3'CUTOFF
    el = prefix.lstrip()
    has_equals = False
    pretag = None

    if "=" in el:
        has_equals = True
        el = el.split("=")[0]
    elif " " in el:
        el = el.split(" ")[0]

    titleComponents = [l.strip().title() for l in el.split("-") if l]
    if len(titleComponents) == 0:
        raise Exception(
            f"Title components for tag '{prefix}' does not have a component"
        )
    titleComponents[0] = titleComponents[0].lower()
    tag = "".join(titleComponents)

    return el, tag, has_equals


def from_container(
    container: str,
    basecommand: List[str],
    helpcommand="-h",
    containersoftware="docker",
    optionsmarker: Optional[str] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
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

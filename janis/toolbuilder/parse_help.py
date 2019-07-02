from typing import Optional, Union, List

import janis.toolbuilder.cltconvert as clt
from janis import ToolMetadata
from janis.tool.tool import ToolInput
from janis.types.common_data_types import String
from janis.utils.logger import Logger


def from_docker(
    docker: str,
    basecommand: Union[str, List[str]],
    help_param: Optional[str] = "--help",
):
    import subprocess, os

    bc = base_command if isinstance(basecommand, list) else [base_command]
    cmd = ["docker", "run", docker, *bc]

    if help_param:
        cmd.append(help_param)

    print("Running command: " + " ".join(f"'{x}'" for x in cmd))

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, preexec_fn=os.setsid, stderr=subprocess.STDOUT
    )
    Logger.info("Starting docker process on pid=" + str(process.pid))
    help = process.communicate()[0].decode("utf-8").rstrip()
    print(help)
    return help


def parse_str(
    help, option_marker: str = "Options:", requires_prev_line_blank_or_param=False
):
    doc = ""
    args = []
    lines = help.split("\n")
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

            if len(tag) == 1:
                print(
                    f"The tag for '{prefix}' was too short, we need you to come up with a new identifier for:"
                )
                print("\t" + tool_doc if tool_doc else line)
                tag = str(input("New identifier: "))
            try:
                prev_arg = ToolInput(
                    tag,
                    String(),
                    prefix=prefix + eqifrequired,
                    separate_value_from_prefix=not has_equal,
                    doc=tool_doc,
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


if __name__ == "__main__":

    # name = str(input("Name: "))
    # friendly_name = str(input(f"Friendly Name (default='{name}'): "))
    # base_command = str(input("BaseCommand (case sensitive): "))
    # tool_provider = str(input("Tool Provider: "))
    # docker = str(input("Docker: "))

    name = "StrelkaSomatic"
    friendly_name = "StrelkaSomatic"
    base_command = "configureStrelkaSomaticWorkflow.py"
    tool_provider = "Illumina"
    docker = "michaelfranklin/strelka:2.9.10"

    help_str = from_docker(docker, base_command, help_param="--allHelp")
    doc, args = parse_str(help_str)

    st = clt.convert_command_tool_fragments(
        name,
        base_command,
        friendly_name,
        tool_provider,
        args,
        [],
        ToolMetadata(documentation=doc),
    )
    print(st)

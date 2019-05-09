from datetime import datetime
from typing import Optional

from janis.tool.tool import ToolInput
from janis.types.common_data_types import String

from janis.utils.logger import Logger


template = """
from datetime import datetime
from janis import CommandTool, ToolInput, ToolOutput, File, Boolean, String, Int, InputSelector, Filename

class {name}Base(CommandTool):

    def friendly_name(self) -> str:
        return "{friendly_name}"

    @staticmethod
    def tool_provider():
        return "{tool_provider}"

    @staticmethod
    def tool() -> str:
        return "{name}"

    @staticmethod
    def base_command():
        return "{base_command}"
    
    def inputs(self):
        return [
{inputs}
        ]
        
    def outputs(self):
        return [
{outputs}
        ]
        
    def metadata(self):
        return ToolMetadata(
            creator=None, 
            maintainer=None, maintainer_email=None,
            date_created=datetime({cy}, {cm}, {cd}), date_updated=datetime({cy}, {cm}, {cd}),
            institution=None, doi=None,
            citation=None,
            keywords=["{name}"],
            documentation_url="{url}",
            documentation=\"""{doc}""\")
"""


def from_docker(docker, basecommand, help_param:Optional[str]="--help"):
    import subprocess, os

    bc = base_command if isinstance(basecommand, list) else [base_command]
    cmd = ["docker", "run", docker, *bc]

    if help_param:
        cmd.append(help_param)

    print("Running command: " + " ".join(f"'{x}'" for x in cmd))

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, preexec_fn=os.setsid, stderr=subprocess.STDOUT)
    Logger.info("Starting docker process on pid=" + str(process.pid))
    help = process.communicate()[0].decode("utf-8").rstrip()
    print(help)
    print("end-help")
    return help


def parse_str(help):
    doc = ""
    args = []
    lines = help.split("\n")
    options_idx = None
    for il in range(len(lines)):
        line = lines[il]
        if not line.lstrip(): continue

        if line.startswith("Options:"):
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

        if last_line_was_blank_or_param and line_args[0].startswith("-"):
            # sometimes this section has two items
            tags = sorted(
                [get_tag_and_cleanup_prefix(p) for p in line_args[0].split(",")],
                key=lambda l: len(l[1]), reverse=True
            )

            if len(tags) > 1:
                tool_doc += "(" + ", ".join(t[0] for t in tags[1:]) + ") "

            if largs > 1:
                tool_doc += " ".join(line_args[1:])

            prefix, tag, has_equal = tags[0]
            eqifrequired = "=" if has_equal else ""

            if len(tag) == 1:
                print(f"The tag for '{prefix}' was too short, we need you to come up with a new identifier for:")
                print("\t" + tool_doc if tool_doc else line)
                tag = str(input("New identifier: "))
            try:
                prev_arg = ToolInput(
                    tag,
                    String(),
                    prefix=prefix + eqifrequired,
                    separate_value_from_prefix=not has_equal,
                    doc=tool_doc
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
        raise Exception(f"Title components for tag '{prefix}' does not have a component")
    titleComponents[0] = titleComponents[0].lower()
    tag = "".join(titleComponents)

    return el, tag, has_equals



if __name__ == "__main__":

    name = "SamToolsView" # str(input("Name: "))
    friendly_name = "samtools view" # str(input(f"Friendly Name (default='{name}'): "))
    base_command = ['samtools', 'view']   # str(input("BaseCommand (case sensitive): "))
    tool_provider = "Samtools" # str(input("BaseCommand (case sensitive): "))
    docker = "michaelfranklin/bwasamtools:0.7.17-1.9" # str(input("Docker: "))

    help_str = from_docker(docker, base_command)
    doc, args = parse_str(help_str)

    ins = ['\t\t\tToolInput("{tag}", String(), prefix="{prefix}", doc="{doc}"),'.format(tag=t.tag, prefix=t.prefix, doc=t.doc) for t in args if t and t.tag]

    print(template.format(
        name=name.replace(" ", ""),
        friendly_name=friendly_name if friendly_name else name,
        tool_provider=tool_provider,
        base_command=base_command,
        inputs="\n".join(ins),
        outputs="",
        cy=datetime.now().year,
        cm=datetime.now().month,
        cd=datetime.now().day,
        doc=doc,
        url=""
    ))

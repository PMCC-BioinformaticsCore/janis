import argparse, sys

from toolbuilder.parse_help import from_container


def process_args():
    cmds = {"fromcontainer": do_docker}

    parser = argparse.ArgumentParser(description="Execute a workflow")
    subparsers = parser.add_subparsers(help="subcommand help", dest="command")
    parser.add_argument("-d", "--debug", action="store_true")

    subparsers.add_parser("version")
    add_fromcontainer_args(subparsers.add_parser("fromcontainer"))
    # add_workflow_args(subparsers.add_parser("run-workflow"))

    args = parser.parse_args()
    return cmds[args.command](args)


def do_docker(args):
    tool, helpstr = from_container(
        container=args.container,
        basecommand=args.basecommand,
        helpcommand=args.help_str,
        containersoftware=args.container_tool,
        name=args.name,
        optionsmarker=args.options_marker,
        version=args.version,
    )

    if args.printhelp:
        print(helpstr, file=sys.stderr)

    if args.printtool:
        print(tool, file=sys.stderr)

    print(tool, file=sys.stdout)


def add_fromcontainer_args(parser):
    parser.description = (
        "Attempts to parse the help (-h) guide of a tool and convert it into a "
        "CommandTool representation of Janis. The output of the tool is returned to stdout."
    )

    parser.add_argument("container", help="container to run")
    parser.add_argument("basecommand", help="The command of your tool", nargs="+")
    parser.add_argument("--name", help="Name of tool, will default to name of command")
    parser.add_argument(
        "--printhelp", help="Print the output of -h to stderr", action="store_true"
    )
    parser.add_argument(
        "--printtool", help="Print the tool to stderr", action="store_true"
    )
    parser.add_argument(
        "--version",
        help="Version of this tool, will try to get it from the container tag if not present",
        action="store_true",
    )
    parser.add_argument("--container-tool", default="docker")
    parser.add_argument(
        "--options-marker",
        default="Options:",
        help="There's usually a header that separates the documentation from the parameters.",
    )
    parser.add_argument(
        "--help-str", help="String that your tool uses get the help guide", default="-h"
    )


if __name__ == "__main__":
    process_args()

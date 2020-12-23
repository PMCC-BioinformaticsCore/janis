import argparse, sys, os, subprocess

from toolbuilder.parse_help import from_container
from toolbuilder.templates import ToolTemplateType
from toolbuilder.runtest.runner import (
    run_test_case,
    find_test_cases,
    update_status,
    UpdateStatusOption,
    NotificationOption,
    cli_logging,
    send_slack_notification,
    execute,
    add_runtest_args,
)

from janis_assistant.management.configuration import JanisConfiguration


def process_args():
    cmds = {"container": do_container, "run-test": do_runtest}

    parser = argparse.ArgumentParser(description="Execute a workflow")
    subparsers = parser.add_subparsers(help="subcommand help", dest="command")
    parser.add_argument("-d", "--debug", action="store_true")

    subparsers.add_parser("version")
    add_container_args(subparsers.add_parser("container"))
    add_runtest_args(subparsers.add_parser("run-test"))
    # add_workflow_args(subparsers.add_parser("run-workflow"))

    args = parser.parse_args()
    return cmds[args.command](args)


def do_container(args):
    tooltype = ToolTemplateType.base
    if args.gatk4:
        tooltype = ToolTemplateType.gatk4

    outputdir = args.output

    (tool, toolversion), helpstr = from_container(
        container=args.container,
        basecommand=args.basecommand,
        helpcommand=args.help_str,
        containersoftware=args.container_tool,
        name=args.name,
        optionsmarker=args.options_marker,
        version=args.version,
        type=tooltype,
    )

    if args.printhelp:
        print(helpstr, file=sys.stderr)

    if args.printtool:
        print(tool, file=sys.stderr)

    if outputdir:
        from os import makedirs, path

        makedirs(outputdir, exist_ok=True)

        with open(path.join(outputdir, "base.py"), "w+") as f:
            f.write(tool)
        with open(path.join(outputdir, "versions.py"), "w+") as f:
            f.write(toolversion)
        with open(path.join(outputdir, "__init__.py"), "w+"):
            pass
    else:
        print(tool, file=sys.stdout)
        print(toolversion, file=sys.stdout)


def add_container_args(parser):
    parser.description = (
        "Attempts to parse the help (-h) guide of a tool and convert it into a "
        "CommandTool representation of Janis. The output of the tool is returned to stdout."
    )

    parser.add_argument("container", help="container to run")
    parser.add_argument("basecommand", help="The command of your tool", nargs="+")
    parser.add_argument(
        "-o",
        "--output",
        help="Directory to output base.py and versions.py to, otherwise these are both written to stdout",
    )

    annotations = parser.add_argument_group("Annotations")
    annotations.add_argument(
        "--version",
        help="Version of this tool, will try to get it from the container tag if not present",
        action="store_true",
    )
    annotations.add_argument(
        "--name", help="Name of tool, will default to name of command"
    )

    log_options = parser.add_argument_group("Logging options")

    log_options.add_argument(
        "--printhelp", help="Print the output of -h to stderr", action="store_true"
    )
    log_options.add_argument(
        "--printtool", help="Print the tool to stderr", action="store_true"
    )

    parser_info = parser.add_argument_group("Parsing options")
    parser_info.add_argument(
        "--options-marker",
        default="Options:",
        help="There's usually a header that separates the documentation from the parameters.",
    )
    parser_info.add_argument(
        "--help-str", help="String that your tool uses get the help guide", default="-h"
    )

    parser_info.add_argument("--container-tool", default="docker", choices=["docker"])

    extra = parser.add_argument_group("Extra options")
    extra.add_argument(
        "--gatk4", action="store_true", help="Use the GATK4 tool template"
    )


def do_runtest(args):
    config = None
    if args.config:
        config = JanisConfiguration.initial_configuration(path=args.config)

    from toolbuilder.runtest import runner

    runner_path = runner.__file__

    cli_args = sys.argv[2:]
    run_test_commands = ["python", runner_path] + cli_args

    precommands = []
    if config and config.template.id == "spartan":
        precommands = ["sbatch", "-p", "snowy", "--wrap"]
        commands = precommands + [" ".join(run_test_commands)]
    else:
        commands = run_test_commands

    subprocess.run(commands)


if __name__ == "__main__":
    process_args()

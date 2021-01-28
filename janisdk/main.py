import argparse, sys, os, subprocess

from janis_core import Logger

from janisdk.container import do_container, add_container_args
from janisdk.fromcwl import do_fromcwl, add_fromcwl_args
from janisdk.runtest import runner as test_runner

from janis_assistant.management.configuration import JanisConfiguration


def process_args():
    cmds = {"container": do_container, "run-test": do_runtest, "fromcwl": do_fromcwl}

    parser = argparse.ArgumentParser(description="Execute a workflow")
    subparsers = parser.add_subparsers(help="subcommand help", dest="command")
    parser.add_argument("-d", "--debug", action="store_true")

    subparsers.add_parser("version")
    add_container_args(subparsers.add_parser("container"))
    test_runner.add_runtest_args(subparsers.add_parser("run-test"))
    add_fromcwl_args(subparsers.add_parser("fromcwl"))

    args = parser.parse_args()
    return cmds[args.command](args)


def do_runtest(args):
    config = None
    if args.config:
        config = JanisConfiguration.initial_configuration(path=args.config)

    runner_path = test_runner.__file__

    cli_args = sys.argv[2:]
    run_test_commands = ["python", runner_path] + cli_args

    if config:
        commands = config.template.template.prepare_run_test_command(run_test_commands)
    else:
        commands = run_test_commands

    joined_command = "' '".join(commands)
    Logger.info(f"Deploying test with command: '{joined_command}'")
    subprocess.run(commands)


if __name__ == "__main__":
    process_args()

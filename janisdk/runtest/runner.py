import ast
import requests
from typing import List, Dict, Any, Optional
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool import test_helpers
from janis_core import Logger
from janis_assistant.engines.enginetypes import EngineType


class UpdateStatusOption:
    def __init__(self, url: str, token: str, method: Optional[str] = "patch"):
        self.url = url
        self.token = token
        self.method = method


class NotificationOption:
    def __init__(
        self, url: str, tool_name: str, test_case: str, test_id: Optional[str] = None
    ):
        self.url = url
        self.tool_name = tool_name
        self.test_case = test_case
        self.test_id = test_id


class TestCasesNotFound(Exception):
    pass


def run_test_case(
    tool_id: str,
    test_case: str,
    engine: EngineType,
    output: Optional[Dict] = None,
    config: str = None,
) -> Dict[str, Any]:
    tool = test_helpers.get_one_tool(tool_id)

    if not tool:
        raise Exception(f"Tool {tool_id} not found")

    runner = ToolTestSuiteRunner(tool, config=config)
    tests_to_run = [tc for tc in tool.tests() if tc.name.lower() == test_case.lower()]

    if not tests_to_run:
        raise Exception(f"Test case {test_case} not found")

    if len(tests_to_run) > 1:
        raise Exception(
            f"There is more than one test case with the same name {test_case}"
        )

    if output is not None:
        Logger.info(
            "Dryrun: validating test using provided output data without running the workflow"
        )

    failed = set()
    succeeded = set()
    execution_error = ""

    try:
        failed, succeeded, output = runner.run_one_test_case(
            t=tests_to_run[0], engine=engine, output=output
        )
    except Exception as e:
        execution_error = str(e)
    except SystemExit as e:
        execution_error = f"Workflow execution failed (exit code: {e.code})"

    return {
        "failed": list(failed),
        "succeeded": list(succeeded),
        "output": output,
        "execution_error": execution_error,
    }


def find_test_cases(tool_id: str):
    tool = test_helpers.get_one_tool(tool_id)

    if not tool:
        raise Exception(f"Tool {tool_id} not found")

    if tool.tests() is None:
        raise TestCasesNotFound(
            f"No test cases found for Tool {tool_id}. "
            f"You must implement the tests() function for this tool."
        )

    return [tc.name for tc in tool.tests()]


def update_status(result: Dict, option: UpdateStatusOption):
    Logger.info(f"Updating test status via {option.method} {option.url}")

    status = "test-failed"
    if not len(result["failed"]) and not result["execution_error"]:
        status = "test-succeeded"

    data = {"status": status, **result}

    headers = {"Authorization": f"Bearer {option.token}"}
    resp = requests.request(
        method=option.method, url=option.url, json=data, headers=headers
    )

    Logger.info("status updated")
    Logger.info(f"Response code {resp.status_code}")
    Logger.info(f"Response:\n{resp.text}")

    return resp.status_code, resp.text


def send_slack_notification(result: Dict, option: NotificationOption):
    Logger.info("sending notification to Slack")

    if len(result["failed"]) == 0 and not result["execution_error"]:
        failed = False
        status = "Test Succeeded"
        icon = ":white_check_mark:"
    else:
        failed = True
        status = "Test Failed"
        icon = ":x:"

    test_description = ""
    if option.test_id:
        test_description = f" *{option.test_id}*"

    summary_block = {
        "type": "section",
        "text": {
            "type": "mrkdwn",
            "text": f"{icon} {status}{test_description}: {option.tool_name} - {option.test_case}",
        },
    }

    blocks = [summary_block]

    if failed and result["failed"]:
        failed_expected_output = []

        for f in result["failed"]:
            failed_expected_output.append(f":black_small_square: {f}")

        failed_block = {
            "type": "section",
            "text": {"type": "mrkdwn", "text": "\n".join(failed_expected_output)},
        }

        blocks.append(failed_block)

    if result["execution_error"]:
        text = result["execution_error"].replace("\n", "<br />")
        execution_error_block = {
            "type": "section",
            "text": {"type": "mrkdwn", "text": f"{result['execution_error']}"},
        }

        blocks.append(execution_error_block)

    request = {"blocks": blocks}
    resp = requests.post(url=option.url, json=request)

    if resp.status_code == requests.codes.ok:
        Logger.info("Notification sent")
    else:
        Logger.warn("Failed to send slack notification")
        Logger.warn(f"{resp.status_code}: {resp.text}")

    return resp.status_code, resp.text


def cli_logging(result: Dict):
    name = result["test_case"]
    Logger.info(f"Test Case: {name}")
    Logger.info(f"Output: {result['output']}")

    if result["execution_error"]:
        Logger.critical(result["execution_error"])

    if len(result["succeeded"]) > 0:
        Logger.info(f"{len(result['succeeded'])} expected output PASSED")

        Logger.info("Succeeded expected output:")
        for s in result["succeeded"]:
            Logger.info(s)

    if len(result["failed"]) > 0:
        Logger.critical(f"{len(result['failed'])} expected output FAILED")

        Logger.critical("Failed expected output:")
        for f in result["failed"]:
            Logger.critical(f)

    if len(result["failed"]) == 0 and not result["execution_error"]:
        Logger.info(f"Test SUCCEEDED: {name}")
    else:
        Logger.critical(f"Test FAILED: {name}")


def add_runtest_args(parser):
    parser.add_argument("tool", help="Name of tool to test")

    parser.add_argument(
        "--test-case", help="Name of test case as listed in tool.tests()"
    )

    parser.add_argument("-e", "--engine", help="engine", default=EngineType.cromwell)

    parser.add_argument("-c", "--config", help="Path to janis config")

    parser.add_argument(
        "-o", "--output", help="Dry run test by providing a dictionary of output"
    )

    # For updating test-framework API endpoint
    parser.add_argument(
        "--test-manager-url", help="API endpoint to update run status in Test Manager"
    )

    parser.add_argument(
        "--test-manager-token", help="Authentication token for Test Manager API"
    )

    parser.add_argument(
        "--test-id",
        help="Test identification to be attached to the notification message",
    )

    parser.add_argument(
        "--slack-notification-url", help="Slack webhook to send notifications to"
    )


def execute(args):
    output = None
    if args.output:
        output = ast.literal_eval(args.output)

    try:
        available_test_cases = find_test_cases(args.tool)
        if args.test_case:
            if args.test_case not in available_test_cases:
                raise TestCasesNotFound(
                    f"Test case with name `{args.test_case}` NOT found."
                )
            test_cases = [args.test_case]
        else:
            test_cases = available_test_cases

    except Exception as e:
        Logger.critical("Unexpected error occurred when searching for test cases")
        Logger.critical(str(e))
        exit()

    for tc_name in test_cases:

        result = run_test_case(
            tool_id=args.tool,
            test_case=tc_name,
            engine=args.engine,
            output=output,
            config=args.config,
        )
        result["test_case"] = tc_name
        cli_logging(result)

        try:
            # send output to test framework API
            if args.test_manager_url and args.test_manager_token:
                option = UpdateStatusOption(
                    url=args.test_manager_url, token=args.test_manager_token
                )
                update_status(result, option)
        except Exception as e:
            Logger.warn(f"Failed to update test status to {args.test_manager_url}")

        try:
            # Send notification to Slack
            if args.slack_notification_url:
                option = NotificationOption(
                    url=args.slack_notification_url,
                    tool_name=args.tool,
                    test_case=tc_name,
                    test_id=args.test_id,
                )
                send_slack_notification(result=result, option=option)
        except Exception as e:
            Logger.warn(
                f"Failed to send notifications to Slack {args.slack_notification_url}"
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    add_runtest_args(parser)

    args = parser.parse_args()

    execute(args)

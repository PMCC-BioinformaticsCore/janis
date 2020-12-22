import json
import requests
from typing import List, Dict, Any, Optional
from janis_assistant.main import run_with_outputs
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool import test_helpers
from janis_core import Logger
from janis_core import JanisShed
from janis_assistant.engines.enginetypes import EngineType
import janis_bioinformatics


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

    failed, succeeded, output = runner.run_one_test_case(
        t=tests_to_run[0], engine=engine, output=output
    )

    return {"failed": list(failed), "succeeded": list(succeeded), "output": output}


def find_test_cases(tool_id: str):
    tool = test_helpers.get_one_tool(tool_id)

    if not tool:
        raise Exception(f"Tool {tool_id} not found")

    return [tc.name for tc in tool.tests()]


def update_status(result: Dict, option: UpdateStatusOption):
    Logger.info(f"Updating test status via {option.method} {option.url}")

    status = "test-failed"
    if not len(result["failed"]):
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

    if not len(result["failed"]):
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

    if failed:
        failed_expected_output = []

        for f in result["failed"]:
            failed_expected_output.append(f":black_small_square: {f}")

        failed_block = {
            "type": "section",
            "text": {"type": "mrkdwn", "text": "\n".join(failed_expected_output)},
        }

        blocks.append(failed_block)

    request = {"blocks": blocks}
    resp = requests.post(url=option.url, json=request)

    if resp.status_code == requests.codes.ok:
        Logger.info("Notification sent")
    else:
        Logger.warn("Failed to send slack notification")
        Logger.warn(f"{resp.status_code}: {resp.text}")

    return resp.status_code, resp.text


def cli_logging(name: str, result: Dict):
    Logger.info(f"Test Case: {name}")
    Logger.info(f"Output: {result['output']}")

    if len(result["succeeded"]) > 0:
        Logger.info(f"{len(result['succeeded'])} expected output PASSED")

        Logger.info("Succeeded expected output:")
        for s in result["succeeded"]:
            Logger.info(s)

    if len(result["failed"]) > 0:
        Logger.info(f"{len(result['failed'])} expected output FAILED")

        Logger.info("Failed expected output:")
        for f in result["failed"]:
            Logger.info(f)

    if len(result["failed"]) == 0:
        Logger.info(f"Test SUCCEEDED: {name}")
    else:
        Logger.critical(f"Test FAILED: {name}")

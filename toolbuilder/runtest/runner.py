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


def run_test_case(tool_id: str, test_case: str, engine: EngineType, output: Optional[Dict] = None) -> Dict[str, Any]:
    tool = test_helpers.get_one_tool(tool_id)

    if not tool:
        raise Exception(f"Tool {tool_id} not found")

    runner = ToolTestSuiteRunner(tool)
    tests_to_run = [tc for tc in tool.tests() if tc.name == test_case]

    if not tests_to_run:
        raise Exception(f"Test case {test_case} not found")

    if len(tests_to_run) > 1:
        raise Exception(f"There is more than one test case with the same name {test_case}")

    if output is not None:
        Logger.info("Dryrun: validating test using provided output data without running the workflow")

    failed, succeeded, output = runner.run_one_test_case(t=tests_to_run[0], engine=engine, output=output)

    return {
        "failed": list(failed),
        "succeeded": list(succeeded),
        "output": output
    }


def update_status(result: Dict, option: UpdateStatusOption):
    Logger.info(f"Updating test status via {option.method} {option.url}")

    status = "test-failed"
    if not len(result["failed"]):
        status = "test-succeeded"

    data = {"status": status, **result}

    headers = {"Authorization": f"Bearer {option.token}"}
    resp = requests.request(method=option.method, url=option.url, json=data, headers=headers)

    Logger.info("status updated")
    Logger.info(f"Response code {resp.status_code}")
    Logger.info(f"Response:\n{resp.text}")

    return resp.status_code, resp.text


def cli_logging(result: Dict):
    Logger.info(f"Output: {result['output']}")

    if len(result['succeeded']) > 0:
        Logger.info(f"{len(result['succeeded'])} expected output PASSED")

        Logger.info("Succeeded expected output:")
        for s in result["succeeded"]:
            Logger.info(s)

    if len(result['failed']) > 0:
        Logger.info(f"{len(result['failed'])} expected output FAILED")

        Logger.info("Failed expected output:")
        for f in result["failed"]:
            Logger.info(f)

    else:
        Logger.info("Test SUCCEEDED")

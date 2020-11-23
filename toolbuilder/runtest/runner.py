import requests
from typing import List, Dict, Any, Optional
from janis_assistant.main import run_with_outputs
from janis_core.tool.test_suite_runner import ToolTestSuiteRunner
from janis_core.tool.test_helpers import get_one_tool
from janis_core import Logger
from janis_core import JanisShed
from janis_assistant.engines.enginetypes import EngineType
import janis_bioinformatics


class UpdateStatusOption:

    def __init__(self, url: str, token: str, method: Optional[str] = "patch"):
        self.url = url
        self.token = token
        self.method = method


def run_test_case(tool_id: str, test_case: str, engine: EngineType) -> Dict[str, Any]:
    # shed = JanisShed
    # shed.hydrate(force=True)
    #
    # tool = shed.get_tool(tool=tool_id)()

    tool = get_one_tool(tool_id, modules=[janis_bioinformatics.tools])

    if not tool:
        raise Exception(f"Tool {tool_id} not found")

    runner = ToolTestSuiteRunner(tool)
    tests_to_run = [tc for tc in tool.tests() if tc.name == test_case]

    if not tests_to_run:
        raise Exception(f"Test case {test_case} not found")

    if len(tests_to_run) > 1:
        raise Exception(f"There is more than one test case with the same name {test_case}")

    failed, succeeded, output = runner.run_one_test_case(tests_to_run[0], engine)

    return {
        "failed": list(failed),
        "succeeded": list(succeeded),
        "output": output
    }


def update_status(result: Dict, option: UpdateStatusOption):
    Logger.info(f"Updating test status via {option.method} {option.url}")

    status = "failed"
    if not len(result["failed"]):
        status = "succeeded"

    headers = {"Authorization": f"Bearer {option.token}"}
    data = {"status": status}
    resp = requests.request(method=option.method, url=option.url, json=data, headers=headers)

    Logger.info(f"status updated: {resp.status_code} {resp.text}")

    return resp.status_code, resp.text


def cli_logging(result: Dict):
    Logger.info(f"Output: {result['output']}")

    Logger.info(f"{len(result['succeeded'])} expected output PASSED")
    if len(result['succeeded']) > 0:
        Logger.info("Succeeded expected output:")
        for s in result["succeeded"]:
            Logger.info(s)

    Logger.info(f"{len(result['failed'])} expected output FAILED")
    if len(result['failed']) > 0:
        Logger.info("Failed expected output:")
        for f in result["failed"]:
            Logger.info(f)
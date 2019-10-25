import janis_core as j
from janis.tools import Echo

w = j.WorkflowBuilder("workflowId")

w.input("input_id", j.String, default="Hello, World!")
w.step("step_id", Echo(inp=w.input_id))
w.output("output_id", source=w.step_id.out)

# Will print the CWL, input file and relevant tools to the console
if __name__ == "__main__":
    w.translate("cwl", to_disk=False)  # or "wdl"

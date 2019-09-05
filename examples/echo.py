import janis as j
from janis.unix.tools.echo import Echo

w = j.Workflow("workflowId")

w.input("inputIdentifier", j.String, default="Hello, World!")
w.step("stepIdentifier", Echo, inp=w.inputIdentifier)
w.output("outputIdentifier", source=w.stepIdentifier.out)

# Will print the CWL, input file and relevant tools to the console
if __name__ == "__main__":
    w.translate("cwl", to_disk=False)  # or "wdl"

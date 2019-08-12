import janis as j
from janis.unix.tools.echo import Echo

w = j.Workflow("workflowId")

inp = j.Input("inputIdentifier", data_type=j.String(), value="HelloWorld")
echo = j.Step("stepIdentifier", tool=Echo())
out = j.Output("outputIdentifier")

w.add_edges(
    [
        (inp, echo.inp),  # Connect 'inp' to 'echostep'
        (echo.out, out),  # Connect output of 'echostep' to 'out'
    ]
)

# Will print the CWL, input file and relevant tools to the console
if __name__ == "__main__":
    w.translate("cwl", to_disk=False)  # or "wdl"

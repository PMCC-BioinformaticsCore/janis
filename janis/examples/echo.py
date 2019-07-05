import janis as j
from janis.unix.tools.echo import Echo

w = j.Workflow("workflowId")

inp = j.Input("inputIdentifier", j.String(), value="HelloWorld")
echo = j.Step("stepIdentifier", Echo())
out = j.Output("outputIdentifier")

w.add_edges(
    [
        (inp, echo.inp),  # Connect 'inp' to 'echostep'
        (echo.out, out),  # Connect output of 'echostep' to 'out'
    ]
)

# Will print the CWL, input file and relevant tools to the console
w.translate("cwl", to_disk=False)  # or "wdl"

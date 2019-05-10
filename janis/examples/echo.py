import janis as j
from janis.unix.tools.echo import Echo

w = j.Workflow("echo")

inp = j.Input("input_identifier", j.String(), "Hello, World!")
step = j.Step("step_identifier", Echo())
outp = j.Output("output_identifier")

w.add_pipe(inp, step, outp)

# Will print the CWL, input file and relevant tools to the console
# w.translate("wdl", to_disk=True, write_inputs_file=True)
